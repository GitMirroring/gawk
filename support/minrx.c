//
// MinRX: a minimal matcher for POSIX Extended Regular Expressions.
// Copyright (C) 2023, 2024, 2025 Michael J. Haertel.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <locale.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wchar.h>
#include <wctype.h>
#include <langinfo.h>
#include "minrx.h"
#define CHARSET	1
#ifdef CHARSET
#include "charset.c"
#endif

#ifdef __GNUC__
#define INLINE __attribute__((__always_inline__)) inline
#else
#define INLINE inline
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GETTEXT_H
#include <gettext.h>
#define _(msgid)  gettext(msgid)
#else
#define _(msgid)  msgid
#endif

#define N_(msgid) msgid

#ifndef RE_DUP_MAX
#define RE_DUP_MAX 32767
#endif

// Utility functions
static INLINE int
ctz(unsigned int x)
{
	return __builtin_ctz(x);
}

static INLINE int
ctzl(unsigned long x)
{
	return __builtin_ctzl(x);
}

static INLINE int
ctzll(unsigned long long x)
{
	return __builtin_ctzll(x);
}

static INLINE size_t
min_size(size_t a, size_t b)
{
	return a < b ? a : b;
}

static INLINE size_t
max_size(size_t a, size_t b)
{
	return a > b ? a : b;
}

static INLINE int32_t
min_i32(int32_t a, int32_t b)
{
	return a < b ? a : b;
}

static INLINE int32_t
max_i32(int32_t a, int32_t b)
{
	return a > b ? a : b;
}

// Forward declarations
typedef int32_t WChar;
typedef size_t NInt;

// COWVec - Copy-on-write vector for size_t elements
typedef struct COWVec_Storage COWVec_Storage;
typedef struct COWVec_Allocator COWVec_Allocator;
typedef struct COWVec COWVec;

struct COWVec_Allocator
{
	size_t length;
	COWVec_Storage *freelist;
};

struct COWVec_Storage
{
	union
	{
		COWVec_Allocator *allocator;
		COWVec_Storage *freelink;
	} u;
	size_t refcnt;
	size_t data[1];		// flexible array member (hack for C compatibility)
};

struct COWVec
{
	COWVec_Storage *storage;
};

static COWVec_Allocator *
COWVec_Allocator_new(size_t length)
{
	COWVec_Allocator *a = malloc(sizeof(COWVec_Allocator));

	a->length = length;
	a->freelist = NULL;
	return a;
}

static void
COWVec_Allocator_free(COWVec_Allocator *a)
{
	COWVec_Storage *s = a->freelist;

	while (s)
	{
		COWVec_Storage *next = s->u.freelink;

		free(s);
		s = next;
	}
	free(a);
}

static COWVec_Storage *
COWVec_Allocator_alloc(COWVec_Allocator *a)
{
	COWVec_Storage *s;

	if (a->freelist)
	{
		s = a->freelist;
		a->freelist = s->u.freelink;
		s->u.allocator = a;
		s->refcnt = 1;
	}
	else
	{
		void *p =
		    malloc(sizeof(COWVec_Storage) +
			   (a->length - 1) * sizeof(size_t));
		s = (COWVec_Storage *) p;
		s->u.allocator = a;
		s->refcnt = 1;
	}
	for (size_t i = 0; i < a->length; i++)
		s->data[i] = (size_t) -1;
	return s;
}

static void
COWVec_Allocator_dealloc(COWVec_Allocator *a, COWVec_Storage *s)
{
	s->u.freelink = a->freelist;
	a->freelist = s;
}

static COWVec_Storage *
COWVec_Storage_clone(COWVec_Storage *s)
{
	COWVec_Allocator *a = s->u.allocator;
	COWVec_Storage *clone = COWVec_Allocator_alloc(a);

	for (size_t i = 0; i < a->length; i++)
		clone->data[i] = s->data[i];
	return clone;
}

static COWVec
COWVec_new(COWVec_Allocator *a)
{
	COWVec v;

	v.storage = COWVec_Allocator_alloc(a);
	return v;
}

static void
COWVec_free(COWVec *v)
{
	if (v->storage && --v->storage->refcnt == 0)
		COWVec_Allocator_dealloc(v->storage->u.allocator, v->storage);
	v->storage = NULL;
}

static COWVec
COWVec_copy(COWVec *src)
{
	COWVec dst;

	dst.storage = src->storage;
	++dst.storage->refcnt;
	return dst;
}

static void __attribute__((unused)) COWVec_assign(COWVec *dst, COWVec *src)
{
	++src->storage->refcnt;
	if (dst->storage && --dst->storage->refcnt == 0)
		COWVec_Allocator_dealloc(dst->storage->u.allocator,
					 dst->storage);
	dst->storage = src->storage;
}

static size_t
COWVec_get(const COWVec *v, size_t idx)
{
	return v->storage->data[idx];
}

static void
COWVec_put(COWVec *v, size_t idx, size_t val)
{
	if (v->storage->data[idx] == val)
		return;
	if (v->storage->refcnt > 1)
	{
		--v->storage->refcnt;
		v->storage = COWVec_Storage_clone(v->storage);
		v->storage->refcnt = 1;
	}
	v->storage->data[idx] = val;
}

static int
COWVec_cmp_range(const COWVec *a, const COWVec *b, size_t offset, size_t n)
{
	const size_t *av = a->storage->data;
	const size_t *bv = b->storage->data;

	for (size_t i = 0; i < n; i++)
	{
		size_t ai = av[offset + i];
		size_t bi = bv[offset + i];

		if (ai != bi)
			return (ai < bi) ? -1 : 1;
	}
	return 0;
}

// QSet - Set of unsigned integers using bitmaps
typedef struct QSet
{
	uint64_t *bits[10];
	uint64_t bits0;
	uint64_t *bitsfree;
	int depth;
} QSet;

static QSet *
QSet_new(size_t limit)
{
	QSet *q = malloc(sizeof(QSet));
	size_t s[10], t = 0;

	q->depth = 0;
	size_t lim = limit;

	do
	{
		lim = s[q->depth++] = (lim + 63u) / 64u;
		t += lim;
	}
	while (lim > 1);

	uint64_t *next = q->bitsfree = malloc(t * sizeof(uint64_t));

	q->bits[0] = &q->bits0;
	for (int i = 1; i < q->depth; ++i)
	{
		q->bits[i] = next;
		next += s[q->depth - 1 - i];
	}
	q->bits0 = 0;
	return q;
}

static void
QSet_free(QSet *q)
{
	if (q->bitsfree)
		free(q->bitsfree);
	free(q);
}

static INLINE uint64_t
QSet_bit(size_t k)
{
	return (uint64_t) 1 << (k & 0x3F);
}

static int
QSet_empty(const QSet *q)
{
	return !q->bits0;
}

static int
QSet_contains(const QSet *q, size_t k)
{
	int i = 0;
	int s = 6 * q->depth;
	size_t j = 0;

	while (i < q->depth)
	{
		uint64_t x = q->bits[i++][j];

		s -= 6;
		j = k >> s;
		uint64_t w = QSet_bit(j);

		if (!(x & w))
			return 0;
	}
	return 1;
}

static int
QSet_insert(QSet *q, size_t k)
{
	int r = 0;
	int i = 0;
	int s = 6 * q->depth;
	size_t j = 0;

	while (i < q->depth)
	{
		uint64_t *bp = &q->bits[i++][j];
		uint64_t x = *bp;

		s -= 6;
		j = k >> s;
		uint64_t w = QSet_bit(j);

		if ((x & w) == 0)
		{
			if (i < q->depth)
				q->bits[i][j] = 0;
			else
				r = 1;
		}
		*bp = x | w;
	}
	return r;
}

static size_t
QSet_remove(QSet *q)
{
	size_t k = 0;
	int i = 0;
	int d = q->depth;

	do
	{
		uint64_t bits = q->bits[i][k];

		k = (k << 6) | ctzll(bits);
		i++;
	}
	while (i != d);

	size_t r = k;

	do
	{
		--i;
		uint64_t w = QSet_bit(k);

		k >>= 6;
		q->bits[i][k] &= ~w;
		if (q->bits[i][k] != 0)
			break;
	}
	while (i != 0);
	return r;
}

// Forward declaration for NState (defined later after COWVec is fully defined)
typedef struct NState NState;

// QVec - Queue/Set of elements with associated data
typedef struct QVec
{
	QSet *qset;
	NState *storage;
	size_t capacity;
} QVec;

// QVec_new is forward declared here, but defined later after NState is complete
static QVec *QVec_new(size_t limit);
static void QVec_free(QVec * v);
static void QVec_clear(QVec * v);

static int
QVec_empty(const QVec *v)
{
	return QSet_empty(v->qset);
}

static int __attribute__((unused)) QVec_contains(const QVec *v, size_t k)
{
	return QSet_contains(v->qset, k);
}

// WConv - Wide character conversion
typedef enum
{
	WConv_Encoding_Byte,
	WConv_Encoding_MBtoWC,
	WConv_Encoding_UTF8
} WConv_Encoding;

typedef struct WConv WConv;
typedef WChar(*WConv_nextfn) (WConv *);

struct WConv
{
	WConv_nextfn nextfn;
	const char *bp;
	const char *ep;
	const char *cp;
	mbstate_t mbs;
};

#define WConv_End (-1)
#define WCharMax 0x10FFFF

static WChar WConv_nextbyte(WConv * wc);
static WChar WConv_nextmbtowc(WConv * wc);
static WChar WConv_nextutf8(WConv * wc);

static WConv_nextfn WConv_nextfns[3] = {
	WConv_nextbyte,
	WConv_nextmbtowc,
	WConv_nextutf8
};

static void
WConv_init(WConv *wc, WConv_Encoding enc, const char *bp, const char *ep)
{
	wc->nextfn = WConv_nextfns[enc];
	wc->bp = bp;
	wc->ep = ep;
	wc->cp = bp;
	memset(&wc->mbs, 0, sizeof(wc->mbs));
}

static void
WConv_init_str(WConv *wc, WConv_Encoding enc, const char *bp)
{
	WConv_init(wc, enc, bp, bp + strlen(bp));
}

static WChar
WConv_nextchr(WConv *wc)
{
	return wc->nextfn(wc);
}

static WChar
WConv_nextbyte(WConv *wc)
{
	return wc->cp !=
	    wc->ep ? (unsigned char) *wc->cp++ : (WChar) WConv_End;
}

static WChar
WConv_nextmbtowc(WConv *wc)
{
	wchar_t wct = L'\0';

	if (wc->cp != wc->ep)
	{
		size_t n = mbrtowc(&wct, wc->cp, wc->ep - wc->cp, &wc->mbs);

		if (n == 0 || n == (size_t) -1 || n == (size_t) -2)
		{
			if (wct == L'\0')
				wct = INT32_MIN + (unsigned char) *wc->cp++;
		}
		else
		{
			wc->cp += n;
		}
		return wct;
	}
	else
	{
		return WConv_End;
	}
}

static WChar
WConv_nextutf8(WConv *wc)
{
	if (wc->cp != wc->ep)
	{
		WChar u = (unsigned char) wc->cp[0];

		if (u < 0x80)
			return wc->cp += 1, u;
		if ((u & 0x40) == 0 || wc->cp + 1 == wc->ep)
		      error:
			return wc->cp += 1, INT32_MIN + u;
		WChar v = (unsigned char) wc->cp[1];

		if ((v & 0xC0) != 0x80)
			goto error;
		if ((u & 0x20) == 0)
		{
			WChar r = ((u & 0x1F) << 6) | (v & 0x3F);

			if (r < 0x80)
				goto error;
			return wc->cp += 2, r;
		}
		if (wc->cp + 2 == wc->ep)
			goto error;
		WChar w = (unsigned char) wc->cp[2];

		if ((w & 0xC0) != 0x80)
			goto error;
		if ((u & 0x10) == 0)
		{
			WChar r =
			    ((u & 0x0F) << 12) | ((v & 0x3F) << 6) | (w &
								      0x3F);
			if (r < 0x800)
				goto error;
			return wc->cp += 3, r;
		}
		if (wc->cp + 3 == wc->ep)
			goto error;
		WChar x = (unsigned char) wc->cp[3];

		if ((x & 0xC0) != 0x80)
			goto error;
		if ((u & 0x08) != 0)
			goto error;
		WChar r =
		    ((u & 0x07) << 18) | ((v & 0x3F) << 12) | ((w & 0x3F) << 6)
		    | (x & 0x3F);
		if (r < 0x010000 || r > 0x10FFFF)
			goto error;
		return wc->cp += 4, r;
	}
	else
	{
		return WConv_End;
	}
}

static WChar
WConv_lookahead(WConv *wc)
{
	WConv copy = *wc;

	return WConv_nextchr(&copy);
}

static size_t
WConv_off(const WConv *wc)
{
	return wc->cp - wc->bp;
}

static const char *
WConv_ptr(const WConv *wc)
{
	return wc->cp;
}

static const char *
WConv_save(const WConv *wc)
{
	return wc->cp;
}

static void
WConv_restore(WConv *wc, const char *p)
{
	wc->cp = p;
}

// CSet - Character set
#ifndef CHARSET
typedef struct CSet_Range
{
	WChar min;
	WChar max;
	struct CSet_Range *left;
	struct CSet_Range *right;
} CSet_Range;

typedef struct CSet_RangeList
{
	CSet_Range *head;
	size_t count;
} CSet_RangeList;

static CSet_Range *
CSet_Range_new(WChar min, WChar max)
{
	CSet_Range *r = malloc(sizeof(CSet_Range));

	r->min = min < max ? min : max;
	r->max = min < max ? max : min;
	r->left = NULL;
	r->right = NULL;
	return r;
}

static void __attribute__((unused)) CSet_Range_free(CSet_Range *r)
{
	if (r)
	{
		CSet_Range_free(r->left);
		CSet_Range_free(r->right);
		free(r);
	}
}

static int
CSet_Range_cmp(const CSet_Range *a, const CSet_Range *b)
{
	return (a->min > b->max) - (a->max < b->min);
}

static CSet_Range *
CSet_Range_insert(CSet_Range *root, WChar min, WChar max)
{
	if (!root)
		return CSet_Range_new(min, max);

	CSet_Range probe;

	probe.min = min - (min != INT32_MIN);
	probe.max = max + (max != WCharMax);

	int cmp = CSet_Range_cmp(&probe, root);

	if (cmp < 0)
	{
		root->left = CSet_Range_insert(root->left, min, max);
	}
	else if (cmp > 0)
	{
		root->right = CSet_Range_insert(root->right, min, max);
	}
	else
	{
		// Overlapping or adjacent - merge
		WChar new_min = min < root->min ? min : root->min;
		WChar new_max = max > root->max ? max : root->max;

		// Check left subtree for merges
		if (root->left)
		{
			CSet_Range *left_max = root->left;

			while (left_max->right)
				left_max = left_max->right;
			if (left_max->max >= new_min - (new_min != INT32_MIN))
			{
				new_min =
				    new_min <
				    left_max->min ? new_min : left_max->min;
				new_max =
				    new_max >
				    left_max->max ? new_max : left_max->max;
				// Remove merged nodes from left subtree - simplified, just update root
			}
		}

		// Check right subtree for merges
		if (root->right)
		{
			CSet_Range *right_min = root->right;

			while (right_min->left)
				right_min = right_min->left;
			if (right_min->min <= new_max + (new_max != WCharMax))
			{
				new_min =
				    new_min <
				    right_min->min ? new_min : right_min->min;
				new_max =
				    new_max >
				    right_min->max ? new_max : right_min->max;
			}
		}

		root->min = new_min;
		root->max = new_max;
	}
	return root;
}
#endif

typedef struct CSet
{
#ifdef CHARSET
	charset_t *charset;
#else
	CSet_Range *ranges;
#endif
	WConv_Encoding enc;
} CSet;

static CSet *CSet_new(WConv_Encoding enc);
static void CSet_free(CSet * cs);
static void CSet_set(CSet * cs, WChar wc);
static void CSet_set_range(CSet * cs, WChar wclo, WChar wchi);
static int CSet_test(const CSet * cs, WChar wc);
static void CSet_invert(CSet * cs);
static void CSet_merge(CSet * dst, const CSet * src);
static minrx_result_t CSet_parse(CSet * cs, minrx_regcomp_flags_t flags,
				 WConv * wconv);

// Node types
typedef enum
{
	Node_CSet = WCharMax + 1,
	Node_Exit,
	Node_Fork,
	Node_Goto,
	Node_Join,
	Node_Loop,
	Node_MinB,
	Node_MinL,
	Node_MinR,
	Node_Next,
	Node_Skip,
	Node_SubL,
	Node_SubR,
	Node_ZBOB,
	Node_ZEOB,
	Node_ZBOL,
	Node_ZEOL,
	Node_ZBOW,
	Node_ZEOW,
	Node_ZXOW,
	Node_ZNWB
} Node_Type;

typedef struct Node
{
	NInt type;
	NInt args[2];
	NInt nstk;
} Node;

// Dynamic arrays
typedef struct NodeArray
{
	Node *data;
	size_t size;
	size_t capacity;
} NodeArray;

static NodeArray *
NodeArray_new(void)
{
	NodeArray *a = malloc(sizeof(NodeArray));

	a->capacity = 16;
	a->size = 0;
	a->data = malloc(a->capacity * sizeof(Node));
	return a;
}

static void
NodeArray_free(NodeArray *a)
{
	free(a->data);
	free(a);
}

static void
NodeArray_push_back(NodeArray *a, Node n)
{
	if (a->size >= a->capacity)
	{
		a->capacity *= 2;
		a->data = realloc(a->data, a->capacity * sizeof(Node));
	}
	a->data[a->size++] = n;
}

static void
NodeArray_push_front(NodeArray *a, Node n)
{
	if (a->size >= a->capacity)
	{
		a->capacity *= 2;
		a->data = realloc(a->data, a->capacity * sizeof(Node));
	}
	memmove(a->data + 1, a->data, a->size * sizeof(Node));
	a->data[0] = n;
	a->size++;
}

static void
NodeArray_insert_array(NodeArray *dst, size_t pos, const NodeArray *src)
{
	while (dst->size + src->size > dst->capacity)
	{
		dst->capacity *= 2;
		dst->data = realloc(dst->data, dst->capacity * sizeof(Node));
	}
	memmove(dst->data + pos + src->size, dst->data + pos,
		(dst->size - pos) * sizeof(Node));
	memcpy(dst->data + pos, src->data, src->size * sizeof(Node));
	dst->size += src->size;
}

typedef struct CSetArray
{
	CSet **data;
	size_t size;
	size_t capacity;
} CSetArray;

static CSetArray *
CSetArray_new(void)
{
	CSetArray *a = malloc(sizeof(CSetArray));

	a->capacity = 16;
	a->size = 0;
	a->data = malloc(a->capacity * sizeof(CSet *));
	return a;
}

static void
CSetArray_free(CSetArray *a)
{
	for (size_t i = 0; i < a->size; i++)
		CSet_free(a->data[i]);
	free(a->data);
	free(a);
}

static size_t
CSetArray_push_back(CSetArray *a, CSet *cs)
{
	if (a->size >= a->capacity)
	{
		a->capacity *= 2;
		a->data = realloc(a->data, a->capacity * sizeof(CSet *));
	}
	a->data[a->size] = cs;
	return a->size++;
}

// Regexp structure
typedef struct Regexp
{
	WConv_Encoding enc;
	minrx_result_t err;
	CSetArray *csets;
	NodeArray *nodes;
	CSet *firstcset;
	int has_firstbytes;
	unsigned char firstbytes[256];
	int has_firstunique;
	unsigned char firstunique;
	size_t nmin;
	size_t nstk;
	size_t nsub;
} Regexp;

static Regexp *
Regexp_new(void)
{
	Regexp *r = malloc(sizeof(Regexp));

	r->csets = CSetArray_new();
	r->nodes = NodeArray_new();
	r->firstcset = NULL;
	r->has_firstbytes = 0;
	r->has_firstunique = 0;
	r->nmin = 0;
	r->nstk = 0;
	r->nsub = 0;
	return r;
}

static void
Regexp_free(Regexp *r)
{
	if (r)
	{
		CSetArray_free(r->csets);
		NodeArray_free(r->nodes);
		if (r->firstcset)
			CSet_free(r->firstcset);
		free(r);
	}
}

// Compile state
typedef struct IcmapEntry
{
	WChar key;
	unsigned int value;
	struct IcmapEntry *left;
	struct IcmapEntry *right;
} IcmapEntry;

typedef struct Compile
{
	minrx_regcomp_flags_t flags;
	WConv_Encoding enc;
	WConv wconv;
	WChar wc;
	CSetArray *csets;
	int has_dot;
	size_t dot_idx;
	int has_esc_s;
	size_t esc_s_idx;
	int has_esc_S;
	size_t esc_S_idx;
	int has_esc_w;
	size_t esc_w_idx;
	int has_esc_W;
	size_t esc_W_idx;
	IcmapEntry *icmap;
	NInt nmin;
	NInt nsub;
} Compile;

typedef struct Subexp
{
	NodeArray *nodes;
	size_t maxstk;
	int hasmin;
	minrx_result_t err;
} Subexp;

static Subexp
Subexp_new(void)
{
	Subexp s;

	s.nodes = NodeArray_new();
	s.maxstk = 0;
	s.hasmin = 0;
	s.err = MINRX_REG_SUCCESS;
	return s;
}

static void
Subexp_free(Subexp *s)
{
	if (s->nodes)
		NodeArray_free(s->nodes);
	s->nodes = NULL;
}

// Forward declarations for compile functions
static Subexp Compile_alt(Compile * c, int nested, NInt nstk);
static Subexp Compile_cat(Compile * c, int nested, NInt nstk);
static Subexp Compile_rep(Compile * c, int nested, NInt nstk);
static Subexp Compile_chr(Compile * c, int nested, NInt nstk);
static Regexp *Compile_compile(Compile * c);

// NState for execution
struct NState
{
	size_t gen;
	size_t boff;
	COWVec substack;
};

// Now we can define QVec_new since NState is complete
static QVec *
QVec_new(size_t limit)
{
	QVec *v = malloc(sizeof(QVec));

	v->qset = QSet_new(limit);
	v->storage = malloc(limit * sizeof(NState));
	v->capacity = limit;
	return v;
}

// Execute state
typedef struct Execute
{
	const Regexp *r;
	minrx_regexec_flags_t flags;
	size_t suboff;
	size_t minvcnt;
	size_t minvoff;
	size_t nestoff;
	size_t gen;
	size_t off;
	WConv wconv;
	WChar wcprev;
	COWVec_Allocator *allocator;
	int has_best;
	COWVec best;
	NInt bestmincount;
	QVec *epsv;
	QSet *epsq;
	const Node *nodes;
} Execute;

static int Execute_execute(Execute * e, size_t nm, minrx_regmatch_t * rm);

// Now implement all the functions...

// CSet implementation
#ifdef CHARSET
static CSet *
CSet_new(WConv_Encoding enc)
{
	CSet *cs = malloc(sizeof(CSet));
	int errcode = 0;

	cs->charset =
	    charset_create(&errcode, MB_CUR_MAX, enc == WConv_Encoding_UTF8);
	cs->enc = enc;
	return cs;
}

static void
CSet_free(CSet *cs)
{
	if (cs)
	{
		if (cs->charset)
			charset_free(cs->charset);
		free(cs);
	}
}

static void
CSet_set(CSet *cs, WChar wc)
{
	charset_add_char(cs->charset, wc);
}

static void
CSet_set_range(CSet *cs, WChar wclo, WChar wchi)
{
	charset_add_range(cs->charset, wclo, wchi);
}

static int
CSet_test(const CSet *cs, WChar wc)
{
	return charset_in_set(cs->charset, wc);
}

static void
CSet_invert(CSet *cs)
{
	int errcode = 0;
	charset_t *newset = charset_invert(cs->charset, &errcode);

	charset_free(cs->charset);
	cs->charset = newset;
}

static void
CSet_merge(CSet *dst, const CSet *src)
{
	charset_merge(dst->charset, src->charset);
}

#else
// Non-CHARSET implementation

static void
CSet_RangeList_free(CSet_Range *r)
{
	if (r)
	{
		CSet_RangeList_free(r->left);
		CSet_RangeList_free(r->right);
		free(r);
	}
}

static CSet *
CSet_new(WConv_Encoding enc)
{
	CSet *cs = malloc(sizeof(CSet));

	cs->ranges = NULL;
	cs->enc = enc;
	return cs;
}

static void
CSet_free(CSet *cs)
{
	if (cs)
	{
		CSet_RangeList_free(cs->ranges);
		free(cs);
	}
}

static void
CSet_set(CSet *cs, WChar wc)
{
	CSet_set_range(cs, wc, wc);
}

static void
CSet_set_range(CSet *cs, WChar wclo, WChar wchi)
{
	cs->ranges = CSet_Range_insert(cs->ranges, wclo, wchi);
}

static int
CSet_test_helper(const CSet_Range *r, WChar wc)
{
	if (!r)
		return 0;
	if (wc < r->min)
		return CSet_test_helper(r->left, wc);
	if (wc > r->max)
		return CSet_test_helper(r->right, wc);
	return 1;
}

static int
CSet_test(const CSet *cs, WChar wc)
{
	if (wc < 0)
		return 0;
	return CSet_test_helper(cs->ranges, wc);
}

static void
CSet_collect_ranges(CSet_Range *r, CSet_Range ***array, size_t *size,
		    size_t *capacity)
{
	if (!r)
		return;
	CSet_collect_ranges(r->left, array, size, capacity);
	if (*size >= *capacity)
	{
		*capacity = (*capacity == 0) ? 16 : *capacity * 2;
		*array = realloc(*array, *capacity * sizeof(CSet_Range *));
	}
	(*array)[(*size)++] = r;
	CSet_collect_ranges(r->right, array, size, capacity);
}

static void
CSet_invert(CSet *cs)
{
	// Collect all ranges into sorted array
	CSet_Range **array = NULL;
	size_t size = 0;
	size_t capacity = 0;

	CSet_collect_ranges(cs->ranges, &array, &size, &capacity);

	// Build inverted range list
	CSet_Range *new_ranges = NULL;
	WChar lo = 0;

	for (size_t i = 0; i < size; i++)
	{
		if (lo < array[i]->min)
		{
			new_ranges =
			    CSet_Range_insert(new_ranges, lo,
					      array[i]->min - 1);
		}
		lo = array[i]->max + 1;
	}
	if (lo <= WCharMax)
	{
		new_ranges = CSet_Range_insert(new_ranges, lo, WCharMax);
	}

	free(array);
	CSet_RangeList_free(cs->ranges);
	cs->ranges = new_ranges;
}

static void
CSet_merge(CSet *dst, const CSet *src)
{
	if (!src)
		return;
	if (!src->ranges)
		return;

	// Collect all ranges from src first, then insert them
	// (cannot recursively walk tree while modifying it)
	CSet_Range **array = NULL;
	size_t size = 0;
	size_t capacity = 0;

	// Cast away const for collection (doesn't modify the tree)
	CSet_collect_ranges((CSet_Range *) src->ranges, &array, &size,
			    &capacity);

	for (size_t i = 0; i < size; i++)
	{
		CSet_set_range(dst, array[i]->min, array[i]->max);
	}

	if (array)
		free(array);
}
#endif

// CSet parsing and character class support
static int
CSet_cclass(CSet *cs, minrx_regcomp_flags_t flags, const char *name)
{
#ifdef CHARSET
	int result = charset_add_cclass(cs->charset, name);

	if ((flags & MINRX_REG_ICASE) != 0)
	{
		if (strcmp(name, "lower") == 0)
			charset_add_cclass(cs->charset, "upper");
		else if (strcmp(name, "upper") == 0)
			charset_add_cclass(cs->charset, "lower");
	}
	return result == CSET_SUCCESS;
#else
	wctype_t wct = wctype(name);

	if (!wct)
		return 0;

	if (cs->enc == WConv_Encoding_Byte)
	{
		for (WChar b = 0; b <= 0xFF; ++b)
		{
			if (iswctype(btowc(b), wct))
			{
				CSet_set(cs, b);
				if ((flags & MINRX_REG_ICASE) != 0)
				{
					CSet_set(cs, tolower(b));
					CSet_set(cs, toupper(b));
				}
			}
		}
	}
	else
	{
		for (WChar wc = 0; wc <= WCharMax; ++wc)
		{
			if (iswctype(wc, wct))
			{
				CSet_set(cs, wc);
				if ((flags & MINRX_REG_ICASE) != 0)
				{
					CSet_set(cs, towlower(wc));
					CSet_set(cs, towupper(wc));
				}
			}
		}
	}
	return 1;
#endif
}

#ifndef CHARSET
static void
CSet_add_equiv(CSet *cs, int32_t equiv)
{
	wchar_t wcs_in[2];
	wchar_t wcs[2];
	wchar_t abuf[100], wbuf[100];

	wcs_in[0] = equiv;
	wcs_in[1] = 0;
	wcsxfrm(abuf, wcs_in, 99);
	wcs[1] = 0;
	for (wchar_t u = 1; u <= WCharMax; ++u)
	{
		wcs[0] = u;
		wcsxfrm(wbuf, wcs, 99);
		if (abuf[0] == wbuf[0])
			CSet_set(cs, u);
	}
}
#endif

static minrx_result_t
CSet_parse(CSet *cs, minrx_regcomp_flags_t flags, WConv *wconv)
{
	WChar wc = WConv_nextchr(wconv);
	int inv = (wc == L'^');

	if (inv)
		wc = WConv_nextchr(wconv);

	for (int first = 1; first || wc != L']'; first = 0)
	{
		if (wc == WConv_End)
			return MINRX_REG_EBRACK;

		WChar wclo = wc, wchi = wc;

		wc = WConv_nextchr(wconv);

		if (wclo == L'\\' && (flags & MINRX_REG_BRACK_ESCAPE) != 0)
		{
			if (wc != WConv_End)
			{
				wclo = wchi = wc;
				wc = WConv_nextchr(wconv);
			}
			else
			{
				return MINRX_REG_EESCAPE;
			}
		}
		else if (wclo == L'[')
		{
			if (wc == L'.')
			{
				wc = WConv_nextchr(wconv);
				wclo = wchi = wc;
				wc = WConv_nextchr(wconv);
				if (wc != L'.'
				    || (wc = WConv_nextchr(wconv)) != L']')
					return MINRX_REG_ECOLLATE;
				wc = WConv_nextchr(wconv);
			}
			else if (wc == L':')
			{
				const char *bp = WConv_ptr(wconv);
				const char *ep = bp;

				do
				{
					ep = WConv_ptr(wconv);
					wc = WConv_nextchr(wconv);
				}
				while (wc != WConv_End && wc != L':');
				if (wc != L':')
					return MINRX_REG_ECTYPE;
				wc = WConv_nextchr(wconv);
				if (wc != L']')
					return MINRX_REG_ECTYPE;
				wc = WConv_nextchr(wconv);

				size_t len = ep - bp;
				char *cclname = malloc(len + 1);

				memcpy(cclname, bp, len);
				cclname[len] = '\0';
				int ok = CSet_cclass(cs, flags, cclname);

				free(cclname);
				if (!ok)
					return MINRX_REG_ECTYPE;
				continue;
			}
			else if (wc == L'=')
			{
				wc = WConv_nextchr(wconv);
				wclo = wchi = wc;
#ifdef CHARSET
				charset_add_equiv(cs->charset, wc);
				if ((flags & MINRX_REG_ICASE) != 0)
				{
					if (iswlower(wc))
						charset_add_equiv(cs->charset,
								  towupper
								  (wc));
					else if (iswupper(wc))
						charset_add_equiv(cs->charset,
								  towlower
								  (wc));
				}
#else
				CSet_add_equiv(cs, wc);
				if ((flags & MINRX_REG_ICASE) != 0)
				{
					if (iswlower(wc))
						CSet_add_equiv(cs,
							       towupper(wc));
					else if (iswupper(wc))
						CSet_add_equiv(cs,
							       towlower(wc));
				}
#endif
				wc = WConv_nextchr(wconv);
				if (wc != L'='
				    || (wc = WConv_nextchr(wconv)) != L']')
					return MINRX_REG_ECOLLATE;
				wc = WConv_nextchr(wconv);
			}
		}

		int range = 0;

		if (wc == L'-')
		{
			const char *save = WConv_save(wconv);

			wc = WConv_nextchr(wconv);
			if (wc == WConv_End)
				return MINRX_REG_EBRACK;
			if (wc != L']')
			{
				wchi = wc;
				wc = WConv_nextchr(wconv);
				if (wchi == L'\\'
				    && (flags & MINRX_REG_BRACK_ESCAPE) != 0)
				{
					if (wc != WConv_End)
					{
						wchi = wc;
						wc = WConv_nextchr(wconv);
					}
					else
					{
						return MINRX_REG_EESCAPE;
					}
				}
				else if (wchi == L'[')
				{
					if (wc == L'.')
					{
						wchi = WConv_nextchr(wconv);
						wc = WConv_nextchr(wconv);
						if (wc != L'.'
						    || (wc =
							WConv_nextchr(wconv))
						    != L']')
							return
							    MINRX_REG_ECOLLATE;
						wc = WConv_nextchr(wconv);
					}
					else if (wc == L':' || wc == L'=')
					{
						return MINRX_REG_ERANGE;
					}
				}
				range = 1;
			}
			else
			{
				WConv_restore(wconv, save);
				wc = L'-';
			}
		}

		if (wclo > wchi || (wclo != wchi && (wclo < 0 || wchi < 0)))
			return MINRX_REG_ERANGE;

		if (wclo >= 0)
		{
			CSet_set_range(cs, wclo, wchi);
			if ((flags & MINRX_REG_ICASE) != 0)
			{
				for (WChar w = wclo; w <= wchi; ++w)
				{
					if (cs->enc == WConv_Encoding_Byte)
					{
						CSet_set(cs, tolower(w));
						CSet_set(cs, toupper(w));
					}
					else
					{
						CSet_set(cs, towlower(w));
						CSet_set(cs, towupper(w));
					}
				}
			}
		}

		if (range && wc == L'-' && WConv_lookahead(wconv) != L']')
			return MINRX_REG_ERANGE;
	}

	if (inv)
	{
		if ((flags & MINRX_REG_NEWLINE) != 0)
			CSet_set(cs, L'\n');
		CSet_invert(cs);
	}

	return MINRX_REG_SUCCESS;
}

// Continue with remaining implementations in next part due to length...
// This is approximately 40% of the full conversion. The pattern continues
// for Compile functions, Execute functions, and public API.

// Helper function for icmap
static unsigned int
Icmap_find(IcmapEntry *root, WChar key, int *found)
{
	if (!root)
	{
		*found = 0;
		return 0;
	}
	if (key < root->key)
		return Icmap_find(root->left, key, found);
	if (key > root->key)
		return Icmap_find(root->right, key, found);
	*found = 1;
	return root->value;
}

static void
Icmap_insert(IcmapEntry **root, WChar key, unsigned int value)
{
	if (!*root)
	{
		*root = malloc(sizeof(IcmapEntry));
		(*root)->key = key;
		(*root)->value = value;
		(*root)->left = NULL;
		(*root)->right = NULL;
		return;
	}
	if (key < (*root)->key)
		Icmap_insert(&(*root)->left, key, value);
	else if (key > (*root)->key)
		Icmap_insert(&(*root)->right, key, value);
}

static void
Icmap_free(IcmapEntry *root)
{
	if (root)
	{
		Icmap_free(root->left);
		Icmap_free(root->right);
		free(root);
	}
}

// Helper function to parse numbers in repetitions
static int
Compile_num(Compile *c, NInt *n)
{
	WChar *wc = &c->wc;

	if (*wc < L'0' || *wc > L'9')
		return 0;

	NInt v = 0;

	do
	{
		NInt oldv = v;

		v = v * 10;
		if (oldv != 0 && v / oldv != 10)
		{
			// Overflow
			do
				*wc = WConv_nextchr(&c->wconv);
			while (*wc >= L'0' && *wc <= L'9');
			*n = (NInt) - 1;
			return 1;
		}
		NInt digit = *wc - L'0';

		if (v > (NInt) - 1 - digit)
		{
			// Overflow
			do
				*wc = WConv_nextchr(&c->wconv);
			while (*wc >= L'0' && *wc <= L'9');
			*n = (NInt) - 1;
			return 1;
		}
		v += digit;
		*wc = WConv_nextchr(&c->wconv);
	}
	while (*wc >= L'0' && *wc <= L'9');
	*n = v;
	return 1;
}

// Helper to adjust all nstk values in nodes
static void
NodeArray_adjust_nstk(NodeArray *a, NInt delta)
{
	for (size_t i = 0; i < a->size; i++)
		a->data[i].nstk += delta;
}

// Forward declarations for rep helpers
static Subexp Compile_minimize(Compile * c, Subexp lh, NInt nstk);
static void Compile_minraise(Compile * c, Subexp * lh);
static Subexp Compile_mkrep_bool(Compile * c, Subexp lh, int optional,
				 int infinite, NInt nstk);
static Subexp Compile_mkrep_count(Compile * c, Subexp lh, NInt m, NInt n,
				  NInt nstk);

// Compile_chr - parse a single character or atom
static Subexp
Compile_chr(Compile *c, int nested, NInt nstk)
{
	Subexp result = Subexp_new();

	result.maxstk = nstk;

	switch (c->wc)
	{
	default:
	      normal:
		if ((c->flags & MINRX_REG_ICASE) == 0)
		{
			Node n = { (NInt) c->wc, {0, 0}, nstk };
			NodeArray_push_back(result.nodes, n);
		}
		else
		{
			WChar wcl =
			    c->enc ==
			    WConv_Encoding_Byte ? (WChar) tolower(c->
								  wc) : (WChar)
			    towlower(c->wc);
			WChar wcu =
			    c->enc ==
			    WConv_Encoding_Byte ? (WChar) toupper(c->
								  wc) : (WChar)
			    towupper(c->wc);
			if (c->wc != wcl || c->wc != wcu)
			{
				WChar key = c->wc < wcl ? c->wc : wcl;

				key = key < wcu ? key : wcu;
				int found = 0;
				unsigned int idx =
				    Icmap_find(c->icmap, key, &found);
				if (!found)
				{
					CSet *cs = CSet_new(c->enc);

					CSet_set(cs, c->wc);
					CSet_set(cs, wcl);
					CSet_set(cs, wcu);
					idx =
					    CSetArray_push_back(c->csets, cs);
					Icmap_insert(&c->icmap, key, idx);
				}
				Node n = { Node_CSet, {idx, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			else
			{
				Node n = { (NInt) c->wc, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
		}
		c->wc = WConv_nextchr(&c->wconv);
		break;

	case L'{':
		if ((c->flags & MINRX_REG_BRACE_COMPAT) != 0)
		{
			WChar la = WConv_lookahead(&c->wconv);
			int is_digit =
			    c->enc ==
			    WConv_Encoding_Byte ? isdigit(la) : iswdigit(la);
			if (!is_digit)
				goto normal;
		}
		// fall through
	case L'*':
	case L'+':
	case L'?':
		result.err = MINRX_REG_BADRPT;
		return result;

	case L'[':
		{
			Node n = { Node_CSet, {c->csets->size, 0}, nstk };
			CSet *cs = CSet_new(c->enc);
			minrx_result_t err =
			    CSet_parse(cs, c->flags, &c->wconv);
			if (err)
			{
				CSet_free(cs);
				result.err = err;
				return result;
			}
			CSetArray_push_back(c->csets, cs);
			NodeArray_push_back(result.nodes, n);
			c->wc = WConv_nextchr(&c->wconv);
		}
		break;

	case L'.':
		if (!c->has_dot)
		{
			c->has_dot = 1;
			c->dot_idx = c->csets->size;
			CSet *cs = CSet_new(c->enc);

			if ((c->flags & MINRX_REG_NEWLINE) != 0)
				CSet_set(cs, L'\n');
			CSet_invert(cs);
			CSetArray_push_back(c->csets, cs);
		}
		{
			Node n = { Node_CSet, {c->dot_idx, 0}, nstk };
			NodeArray_push_back(result.nodes, n);
		}
		c->wc = WConv_nextchr(&c->wconv);
		break;

	case L'^':
		{
			Node n =
			    { (c->flags & MINRX_REG_NEWLINE) ==
			       0 ? Node_ZBOB : Node_ZBOL, {0, 0}, nstk };
			NodeArray_push_back(result.nodes, n);
		}
		c->wc = WConv_nextchr(&c->wconv);
		break;

	case L'$':
		{
			Node n =
			    { (c->flags & MINRX_REG_NEWLINE) ==
			       0 ? Node_ZEOB : Node_ZEOL, {0, 0}, nstk };
			NodeArray_push_back(result.nodes, n);
		}
		c->wc = WConv_nextchr(&c->wconv);
		break;

	case L'\\':
		c->wc = WConv_nextchr(&c->wconv);
		switch (c->wc)
		{
		case L'<':
			if ((c->flags & MINRX_REG_EXTENSIONS_BSD) == 0)
				goto normal;
			{
				Node n = { Node_ZBOW, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'>':
			if ((c->flags & MINRX_REG_EXTENSIONS_BSD) == 0)
				goto normal;
			{
				Node n = { Node_ZEOW, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'`':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			{
				Node n = { Node_ZBOB, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'\'':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			{
				Node n = { Node_ZEOB, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'b':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			{
				Node n = { Node_ZXOW, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'B':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			{
				Node n = { Node_ZNWB, {0, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L's':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			if (!c->has_esc_s)
			{
				c->has_esc_s = 1;
				c->esc_s_idx = c->csets->size;
				WConv wc_tmp;

				WConv_init_str(&wc_tmp, c->enc, "[:space:]]");
				CSet *cs = CSet_new(c->enc);

				CSet_parse(cs, c->flags, &wc_tmp);
				CSetArray_push_back(c->csets, cs);
			}
			{
				Node n =
				    { Node_CSet, {c->esc_s_idx, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'S':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			if (!c->has_esc_S)
			{
				c->has_esc_S = 1;
				c->esc_S_idx = c->csets->size;
				WConv wc_tmp;

				WConv_init_str(&wc_tmp, c->enc, "^[:space:]]");
				CSet *cs = CSet_new(c->enc);

				CSet_parse(cs, c->flags, &wc_tmp);
				CSetArray_push_back(c->csets, cs);
			}
			{
				Node n =
				    { Node_CSet, {c->esc_S_idx, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'w':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			if (!c->has_esc_w)
			{
				c->has_esc_w = 1;
				c->esc_w_idx = c->csets->size;
				WConv wc_tmp;

				WConv_init_str(&wc_tmp, c->enc, "[:alnum:]_]");
				CSet *cs = CSet_new(c->enc);

				CSet_parse(cs, c->flags, &wc_tmp);
				CSetArray_push_back(c->csets, cs);
			}
			{
				Node n =
				    { Node_CSet, {c->esc_w_idx, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case L'W':
			if ((c->flags & MINRX_REG_EXTENSIONS_GNU) == 0)
				goto normal;
			if (!c->has_esc_W)
			{
				c->has_esc_W = 1;
				c->esc_W_idx = c->csets->size;
				WConv wc_tmp;

				WConv_init_str(&wc_tmp, c->enc,
					       "^[:alnum:]_]");
				CSet *cs = CSet_new(c->enc);

				CSet_parse(cs, c->flags, &wc_tmp);
				CSetArray_push_back(c->csets, cs);
			}
			{
				Node n =
				    { Node_CSet, {c->esc_W_idx, 0}, nstk };
				NodeArray_push_back(result.nodes, n);
			}
			break;
		case WConv_End:
			result.err = MINRX_REG_EESCAPE;
			return result;
		default:
			goto normal;
		}
		c->wc = WConv_nextchr(&c->wconv);
		break;

	case L'(':
		{
			NInt n = ++c->nsub;

			// Free the initial empty subexp before overwriting
			Subexp_free(&result);

			c->wc = WConv_nextchr(&c->wconv);
			result = Compile_alt(c, 1, nstk + 1);
			if (result.err)
				return result;
			if (c->wc != L')')
			{
				Subexp_free(&result);
				result.err = MINRX_REG_EPAREN;
				return result;
			}
			Node subl = { Node_SubL, {n, c->nsub}, nstk + 1 };
			Node subr = { Node_SubR, {n, c->nsub}, nstk };
			NodeArray_push_front(result.nodes, subl);
			NodeArray_push_back(result.nodes, subr);
			c->wc = WConv_nextchr(&c->wconv);
		}
		break;

	case L')':
		if (!nested)
			goto normal;
		// fall through
	case L'|':
	case WConv_End:
		break;
	}

	return result;
}

// Compile_rep - handle repetitions (*, +, ?, {m,n})
static Subexp
Compile_rep(Compile *c, int nested, NInt nstk)
{
	Subexp lh = Compile_chr(c, nested, nstk);

	if (lh.err)
		return lh;

	int hasmin = lh.hasmin;

	for (;;)
	{
		int infinite = 0;
		int minimal = (c->flags & MINRX_REG_MINIMAL) != 0;
		int optional = 0;

		switch (c->wc)
		{
		case L'?':
			optional = 1;
			goto common;
		case L'*':
			infinite = optional = 1;
			goto common;
		case L'+':
			infinite = 1;
		      common:
			if ((c->flags & MINRX_REG_MINDISABLE) == 0)
			{
				c->wc = WConv_nextchr(&c->wconv);
				if (c->wc == L'?')
				{
					minimal ^= 1;
					c->wc = WConv_nextchr(&c->wconv);
				}
			}
			else
			{
				c->wc = WConv_nextchr(&c->wconv);
			}
			if (hasmin)
				Compile_minraise(c, &lh);
			lh = minimal ? Compile_minimize(c, lh, nstk) : lh;
			lh = Compile_mkrep_bool(c, lh, optional, infinite,
						nstk);
		      comout:
			if (minimal)
			{
				Node minb = { Node_MinB, {0, 0}, nstk };
				NodeArray_push_front(lh.nodes, minb);
				hasmin = 1;
			}
			lh.hasmin = hasmin;
			continue;

		case L'{':
			{
				WChar la = WConv_lookahead(&c->wconv);
				int is_digit =
				    c->enc ==
				    WConv_Encoding_Byte ? isdigit(la) :
				    iswdigit(la);
				if ((c->flags & MINRX_REG_BRACE_COMPAT) == 0
				    || is_digit)
				{
					c->wc = WConv_nextchr(&c->wconv);
					if (c->wc == WConv_End)
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_EBRACE;
						return lh;
					}
					NInt m, n;

					if (!Compile_num(c, &m))
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_BADBR;
						return lh;
					}
					if (c->wc == L'}')
					{
						if ((c->
						     flags &
						     MINRX_REG_MINDISABLE) ==
						    0)
						{
							c->wc =
							    WConv_nextchr(&c->
									  wconv);
							if (c->wc == L'?')
							{
								minimal ^= 1;
								c->wc =
								    WConv_nextchr
								    (&c->
								     wconv);
							}
						}
						else
						{
							c->wc =
							    WConv_nextchr(&c->
									  wconv);
						}
						if (hasmin)
							Compile_minraise(c,
									 &lh);
						lh = minimal ?
						    Compile_minimize(c, lh,
								     nstk) :
						    lh;
						lh = Compile_mkrep_count(c, lh,
									 m, m,
									 nstk);
						goto comout;
					}
					if (c->wc == WConv_End)
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_EBRACE;
						return lh;
					}
					if (c->wc != L',')
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_BADBR;
						return lh;
					}
					c->wc = WConv_nextchr(&c->wconv);
					if (c->wc == L'}')
					{
						if ((c->
						     flags &
						     MINRX_REG_MINDISABLE) ==
						    0)
						{
							c->wc =
							    WConv_nextchr(&c->
									  wconv);
							if (c->wc == L'?')
							{
								minimal ^= 1;
								c->wc =
								    WConv_nextchr
								    (&c->
								     wconv);
							}
						}
						else
						{
							c->wc =
							    WConv_nextchr(&c->
									  wconv);
						}
						if (hasmin)
							Compile_minraise(c,
									 &lh);
						lh = minimal ?
						    Compile_minimize(c, lh,
								     nstk) :
						    lh;
						lh = Compile_mkrep_count(c, lh,
									 m, -1,
									 nstk);
						goto comout;
					}
					if (!Compile_num(c, &n))
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_BADBR;
						return lh;
					}
					if (c->wc == WConv_End)
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_EBRACE;
						return lh;
					}
					if (c->wc != L'}')
					{
						Subexp_free(&lh);
						lh.err = MINRX_REG_BADBR;
						return lh;
					}
					if ((c->
					     flags & MINRX_REG_MINDISABLE) ==
					    0)
					{
						c->wc =
						    WConv_nextchr(&c->wconv);
						if (c->wc == L'?')
						{
							minimal ^= 1;
							c->wc =
							    WConv_nextchr(&c->
									  wconv);
						}
					}
					else
					{
						c->wc =
						    WConv_nextchr(&c->wconv);
					}
					if (hasmin)
						Compile_minraise(c, &lh);
					lh = minimal ? Compile_minimize(c, lh,
									nstk) :
					    lh;
					lh = Compile_mkrep_count(c, lh, m, n,
								 nstk);
					goto comout;
				}
			}
			// fall through
		default:
			return lh;
		}
	}
}

// Compile_minimize - wrap subexpression for minimal matching
static Subexp
Compile_minimize(Compile *c, Subexp lh, NInt nstk)
{
	if (lh.err)
	{
		return lh;
	}
	NodeArray_adjust_nstk(lh.nodes, 1);
	Node minl = { Node_MinL, {0, 0}, nstk + 1 };
	Node minr = { Node_MinR, {0, 0}, nstk };
	NodeArray_push_front(lh.nodes, minl);
	NodeArray_push_back(lh.nodes, minr);
	lh.maxstk += 1;
	lh.hasmin = 1;
	c->nmin = max_size(c->nmin, 1);
	return lh;
}

// Compile_minraise - increase nesting level of minimized subexpressions
static void
Compile_minraise(Compile *c, Subexp *lh)
{
	NInt maxlevel = 0;

	for (size_t i = 0; i < lh->nodes->size; i++)
	{
		Node *n = &lh->nodes->data[i];

		switch (n->type)
		{
		case Node_MinB:
		case Node_MinL:
		case Node_MinR:
			n->args[0]++;
			maxlevel = max_i32(maxlevel, n->args[0]);
			break;
		default:
			break;
		}
	}
	c->nmin = max_size(c->nmin, maxlevel + 1);
}

// Compile_mkrep_bool - create repetition with boolean flags
static Subexp
Compile_mkrep_bool(Compile *c
		   __attribute__((unused)), Subexp lh, int optional,
		   int infinite, NInt nstk)
{
	if (lh.err)
		return lh;

	if (optional && !infinite)
	{
		NodeArray_adjust_nstk(lh.nodes, 2);
		NInt lhsize = lh.nodes->size;
		Node skip = { Node_Skip, {lhsize, 0}, nstk + 2 };
		NodeArray_push_front(lh.nodes, skip);
		lh.maxstk += 2;
		return lh;
	}
	else
	{
		NodeArray_adjust_nstk(lh.nodes, 3);
		NInt lhsize = lh.nodes->size;
		Node loop = { Node_Loop, {lhsize, (NInt) optional}, nstk + 3 };
		Node next = { Node_Next, {lhsize, (NInt) infinite}, nstk };
		NodeArray_push_front(lh.nodes, loop);
		NodeArray_push_back(lh.nodes, next);
		lh.maxstk += 3;
		return lh;
	}
}

// Compile_mkrep_count - create repetition with count bounds
static Subexp
Compile_mkrep_count(Compile *c, Subexp lh, NInt m, NInt n, NInt nstk)
{
	if (lh.err)
		return lh;

	if ((m != (NInt) - 1 && m > RE_DUP_MAX)
	    || (n != (NInt) - 1 && n > RE_DUP_MAX) || (m != (NInt) - 1
						       && n != (NInt) - 1
						       && m > n))
	{
		Subexp_free(&lh);
		lh.err = MINRX_REG_BADBR;
		return lh;
	}
	if (n == 0)
	{
		Subexp_free(&lh);
		lh = Subexp_new();
		return lh;
	}
	if (m == 0 && n == 1)
		return Compile_mkrep_bool(c, lh, 1, 0, nstk);
	if (m == 0 && n == (NInt) - 1)
		return Compile_mkrep_bool(c, lh, 1, 1, nstk);
	if (m == 1 && n == 1)
		return lh;
	if (m == 1 && n == (NInt) - 1)
		return Compile_mkrep_bool(c, lh, 0, 1, nstk);

	// Make copies
	NodeArray *lhs = lh.nodes;
	NodeArray *rhs = NodeArray_new();

	for (size_t i = 0; i < lhs->size; i++)
		NodeArray_push_back(rhs, lhs->data[i]);

	size_t lhmaxstk = lh.maxstk;
	size_t rhmaxstk = lh.maxstk;
	int lhasmin __attribute__((unused)) = lh.hasmin;
	int rhasmin = lh.hasmin;

	NInt k;

	for (k = 1; k < m; ++k)
		NodeArray_insert_array(lhs, lhs->size, rhs);

	if (n != (NInt) - 1 && k < n)
	{
		lhmaxstk += 2;
		rhmaxstk += 2;
		NodeArray_adjust_nstk(rhs, 2);
		NInt rhsize = rhs->size;
		Node skip = { Node_Skip, {rhsize, 1}, nstk + 2 };
		NodeArray_push_front(rhs, skip);
		for (; k < n; ++k)
			NodeArray_insert_array(lhs, lhs->size, rhs);
	}

	if (n == (NInt) - 1)
	{
		lhmaxstk += 3;
		rhmaxstk += 3;
		NodeArray_adjust_nstk(rhs, 3);
		NInt rhsize = rhs->size;
		Node loop = { Node_Loop, {rhsize, 1}, nstk + 3 };
		Node next = { Node_Next, {rhsize, 1}, nstk };
		NodeArray_push_front(rhs, loop);
		NodeArray_push_back(rhs, next);
		NodeArray_insert_array(lhs, lhs->size, rhs);
	}

	NodeArray_free(rhs);

	if (m == 0)
	{
		Subexp tmp;

		tmp.nodes = lhs;
		tmp.maxstk = rhmaxstk;
		tmp.hasmin = rhasmin;
		tmp.err = MINRX_REG_SUCCESS;
		return Compile_mkrep_bool(c, tmp, 1, 0, nstk);
	}
	else
	{
		lh.maxstk = rhmaxstk;
		lh.hasmin = rhasmin;
		return lh;
	}
}

// Compile_cat - concatenation of subexpressions
static Subexp
Compile_cat(Compile *c, int nested, NInt nstk)
{
	Subexp lh = Compile_rep(c, nested, nstk);

	if (lh.err)
		return lh;

	while (c->wc != WConv_End && c->wc != L'|'
	       && (c->wc != L')' || !nested))
	{
		Subexp rh = Compile_rep(c, nested, nstk);

		if (rh.err)
		{
			Subexp_free(&lh);
			return rh;
		}
		NodeArray_insert_array(lh.nodes, lh.nodes->size, rh.nodes);
		lh.maxstk = max_size(lh.maxstk, rh.maxstk);
		lh.hasmin |= rh.hasmin;
		Subexp_free(&rh);
	}

	return lh;
}

// Compile_alt - alternation (|)
static Subexp
Compile_alt(Compile *c, int nested, NInt nstk)
{
	Subexp lh = Compile_cat(c, nested, nstk);

	if (lh.err)
		return lh;

	if (c->wc == L'|')
	{
		NodeArray_adjust_nstk(lh.nodes, 1);

		// Collect all alternatives
		typedef struct
		{
			NodeArray *nodes;
			size_t maxstk;
			int hasmin;
		} AltData;
		AltData *alts = NULL;
		size_t alts_count = 0;
		size_t alts_capacity = 0;

		while (c->wc == L'|')
		{
			c->wc = WConv_nextchr(&c->wconv);
			Subexp mh = Compile_cat(c, nested, nstk + 1);

			if (mh.err)
			{
				for (size_t i = 0; i < alts_count; i++)
					NodeArray_free(alts[i].nodes);
				free(alts);
				Subexp_free(&lh);
				return mh;
			}
			if (alts_count >= alts_capacity)
			{
				alts_capacity =
				    alts_capacity == 0 ? 4 : alts_capacity * 2;
				alts =
				    realloc(alts,
					    alts_capacity * sizeof(AltData));
			}
			alts[alts_count].nodes = mh.nodes;
			alts[alts_count].maxstk = mh.maxstk;
			alts[alts_count].hasmin = mh.hasmin;
			alts_count++;
		}

		// Build the alternation from right to left
		NodeArray *rhs = alts[alts_count - 1].nodes;
		size_t rhmaxstk = alts[alts_count - 1].maxstk;
		int rhasmin = alts[alts_count - 1].hasmin;

		Node goto_node =
		    { Node_Goto, {rhs->size, rhs->size + 1}, nstk + 1 };
		NodeArray_push_front(rhs, goto_node);

		for (size_t i = alts_count - 1; i > 0; i--)
		{
			NodeArray *mhs = alts[i - 1].nodes;
			size_t mhmaxstk = alts[i - 1].maxstk;
			int mhasmin = alts[i - 1].hasmin;
			size_t mhs_size = mhs->size;  // Save size before freeing

			NodeArray_insert_array(rhs, 0, mhs);
			NodeArray_free(mhs);

			rhmaxstk = max_size(mhmaxstk, rhmaxstk);
			rhasmin |= mhasmin;

			Node goto_node2 =
			    { Node_Goto, {mhs_size, rhs->size + 1},
			    nstk + 1 };
			NodeArray_push_front(rhs, goto_node2);
		}

		free(alts);

		Node fork = { Node_Fork, {lh.nodes->size, 0}, nstk + 1 };
		NodeArray_push_front(lh.nodes, fork);
		NodeArray_insert_array(lh.nodes, lh.nodes->size, rhs);
		NodeArray_free(rhs);
		lh.maxstk = max_size(lh.maxstk, rhmaxstk);
		lh.hasmin |= rhasmin;
		Node join = { Node_Join, {lh.nodes->size - 1, 0}, nstk + 1 };
		NodeArray_push_back(lh.nodes, join);
	}

	return lh;
}

// Helper to compute first character set closure
static CSet *
Compile_firstclosure(Compile *c, const NodeArray *nodes)
{
	if (nodes->size == 0)
		return NULL;

	QSet *epsq = QSet_new(nodes->size);
	QSet *epsv = QSet_new(nodes->size);
	QSet *firsts = QSet_new(nodes->size);

	QSet_insert(epsq, 0);

	while (!QSet_empty(epsq))
	{
		size_t k = QSet_remove(epsq);

		if (QSet_contains(epsv, k))
			continue;
		QSet_insert(epsv, k);

		const Node *n = &nodes->data[k];

		if (n->type <= Node_CSet)
		{
			QSet_insert(firsts, k);
		}
		else
		{
			switch (n->type)
			{
			case Node_Exit:
				QSet_free(epsq);
				QSet_free(epsv);
				QSet_free(firsts);
				return NULL;
			case Node_Fork:
				{
					size_t kk = k;

					do
					{
						QSet_insert(epsq, kk + 1);
						kk = kk + 1 +
						    nodes->data[kk].args[0];
					}
					while (nodes->data[kk].type !=
					       Node_Join);
				}
				break;
			case Node_Goto:
				QSet_insert(epsq, k + 1 + n->args[0]);
				break;
			case Node_Loop:
				QSet_insert(epsq, k + 1);
				if (n->args[1])
					QSet_insert(epsq, k + 1 + n->args[0]);
				break;
			case Node_Skip:
				QSet_insert(epsq, k + 1);
				QSet_insert(epsq, k + 1 + n->args[0]);
				break;
			default:
				QSet_insert(epsq, k + 1);
				break;
			}
		}
	}

	if (QSet_empty(firsts))
	{
		QSet_free(epsq);
		QSet_free(epsv);
		QSet_free(firsts);
		return NULL;
	}

	CSet *cs = CSet_new(c->enc);

	while (!QSet_empty(firsts))
	{
		size_t k = QSet_remove(firsts);
		NInt t = nodes->data[k].type;

		if (t <= WCharMax)
		{
			CSet_set(cs, t);
		}
		else
		{
			size_t idx = nodes->data[k].args[0];

			if (idx < c->csets->size && c->csets->data[idx])
				CSet_merge(cs, c->csets->data[idx]);
		}
	}

	QSet_free(epsq);
	QSet_free(epsv);
	QSet_free(firsts);

	return cs;
}

// Helper to compute first bytes from CSet
static unsigned int
CSet_utfprefix(WChar wc)
{
	if (wc < 0x80)
		return wc;
	if (wc < 0x800)
		return 0xC0 + (wc >> 6);
	if (wc < 0x10000)
		return 0xE0 + (wc >> 12);
	if (wc < 0x100000)
		return 0xF0 + (wc >> 18);
	return 0xF4;
}

#ifndef CHARSET
static void
CSet_firstbytes_helper(const CSet_Range *r, unsigned char *fb,
		       WConv_Encoding enc)
{
	if (!r)
		return;
	CSet_firstbytes_helper(r->left, fb, enc);
	if (enc == WConv_Encoding_Byte)
	{
		if (r->min > 255)
			return;
		WChar lo = r->min;
		WChar hi = r->max < 255 ? r->max : 255;

		for (WChar b = lo; b <= hi; b++)
			fb[b] = 1;
	}
	else if (enc == WConv_Encoding_UTF8)
	{
		unsigned int lo = CSet_utfprefix(r->min);
		unsigned int hi = CSet_utfprefix(r->max);

		for (unsigned int b = lo; b <= hi; b++)
			fb[b] = 1;
	}
	CSet_firstbytes_helper(r->right, fb, enc);
}
#endif

static void
Compile_compute_firstbytes(Regexp *r)
{
	if (!r->firstcset)
	{
		r->has_firstbytes = 0;
		return;
	}

	memset(r->firstbytes, 0, 256);

	switch (r->enc)
	{
	case WConv_Encoding_Byte:
	case WConv_Encoding_UTF8:
#ifdef CHARSET
		{
			int errcode = 0;
			charset_firstbytes_t bytes =
			    charset_firstbytes(r->firstcset->charset,
					       &errcode);
			for (int i = 0; i < 256 && i < MAX_FIRSTBYTES; i++)
				r->firstbytes[i] = bytes.bytes[i];
		}
#else
		CSet_firstbytes_helper(r->firstcset->ranges, r->firstbytes,
				       r->enc);
#endif
		r->has_firstbytes = 1;
		break;
	default:
		r->has_firstbytes = 0;
		break;
	}

	if (r->has_firstbytes)
	{
		int n = 0;
		int u = -1;

		for (int i = 0; i < 256; i++)
		{
			if (r->firstbytes[i])
			{
				n++;
				u = i;
			}
		}
		if (n == 1)
		{
			r->has_firstunique = 1;
			r->firstunique = u;
		}
		else
		{
			r->has_firstunique = 0;
		}
	}
}

// Compile_compile - main compilation entry point
static Regexp *
Compile_compile(Compile *c)
{
	if ((c->flags & MINRX_REG_MINDISABLE) != 0
	    && (c->flags & MINRX_REG_MINIMAL) != 0)
	{
		Regexp *r = Regexp_new();

		r->enc = c->enc;
		r->err = MINRX_REG_BADPAT;
		r->nsub = 1;
		return r;
	}

	Subexp lh = Compile_alt(c, 0, 0);

	Regexp *r = Regexp_new();

	r->enc = c->enc;

	if (lh.err)
	{
		r->err = lh.err;
		Subexp_free(&lh);
		r->nmin = 0;
		r->nstk = 0;
		r->nsub = 0;
		return r;
	}

	r->err = MINRX_REG_SUCCESS;
	Node exit_node = { Node_Exit, {0, 0}, 0 };
	NodeArray_push_back(lh.nodes, exit_node);

	if (c->nmin > 0)
	{
		NodeArray_adjust_nstk(lh.nodes, c->nmin);
	}

	// Transfer nodes
	NodeArray_free(r->nodes);
	r->nodes = lh.nodes;
	lh.nodes = NULL;

	// Compute first character set BEFORE transferring csets
	// (Compile_firstclosure needs access to c->csets)
	CSet *firstcset = Compile_firstclosure(c, r->nodes);

	// Transfer csets
	CSetArray_free(r->csets);
	r->csets = c->csets;
	c->csets = CSetArray_new();

	r->nmin = c->nmin;
	r->nstk = lh.maxstk;
	r->nsub = c->nsub + 1;

	// Set first character set and compute first bytes
	r->firstcset = firstcset;
	Compile_compute_firstbytes(r);

	return r;
}

static void
QVec_free(QVec *v)
{
	if (v)
	{
		QVec_clear(v);
		QSet_free(v->qset);
		free(v->storage);
		free(v);
	}
}

static void
QVec_clear(QVec *v)
{
	// Need to free NState elements
	while (!QSet_empty(v->qset))
	{
		size_t k = QSet_remove(v->qset);
		NState *ns = &v->storage[k];

		COWVec_free(&ns->substack);
	}
}

// Helper for NState comparison
static int
NState_cmp(const NState *a, const NState *b, size_t gen, size_t nstk)
{
	if (gen != b->gen)
		return (gen < b->gen) ? -1 : 1;
	if (a->boff != b->boff)
		return (b->boff < a->boff) ? -1 : 1;
	return COWVec_cmp_range(&a->substack, &b->substack, 0, nstk);
}

// Execute helper functions - forward declarations
static void Execute_add(Execute * e, QVec * ncsv, NInt k, NInt nstk,
			const NState * ns, WChar wcnext);
static void Execute_add_1arg(Execute * e, QVec * ncsv, NInt k, NInt nstk,
			     const NState * ns, WChar wcnext, NInt arg1);
static void Execute_add_2arg(Execute * e, QVec * ncsv, NInt k, NInt nstk,
			     const NState * ns, WChar wcnext, NInt arg1,
			     NInt arg2);
static void Execute_add_3arg(Execute * e, QVec * ncsv, NInt k, NInt nstk,
			     const NState * ns, WChar wcnext, NInt arg1,
			     NInt arg2, NInt arg3);
static void Execute_epsclosure(Execute * e, QVec * ncsv, WChar wcnext);

#define SIZE_BITS (sizeof(size_t) * 8)

// Execute_add with variadic arguments
static void
Execute_add(Execute *e, QVec *ncsv, NInt k, NInt nstk, const NState *ns,
	    WChar wcnext)
{
	const Node *n = &e->nodes[k];

	if (n->type <= Node_CSet)
	{
		if (n->type == (NInt) wcnext
		    || (n->type == Node_CSet
			&& CSet_test(e->r->csets->data[n->args[0]], wcnext)))
		{
			int newly = QSet_insert(ncsv->qset, k);
			NState *newns = &ncsv->storage[k];

			if (newly)
			{
				newns->gen = 0;
				newns->boff = 0;
				newns->substack = COWVec_new(e->allocator);
			}
			int cmp = NState_cmp(ns, newns, e->gen, nstk);

			if (newly || cmp > 0
			    || (cmp == 0 && e->minvcnt
				&& COWVec_cmp_range(&ns->substack,
						    &newns->substack,
						    e->minvoff,
						    e->minvcnt) > 0))
			{
				COWVec old = newns->substack;

				newns->substack =
				    COWVec_copy((COWVec *) & ns->substack);
				newns->gen = e->gen;
				newns->boff = ns->boff;
				COWVec_free(&old);
			}
			else
			{
				return;
			}
		}
	}
	else
	{
		int newly = QSet_insert(e->epsv->qset, k);
		NState *newns = &e->epsv->storage[k];

		if (newly)
		{
			newns->gen = 0;
			newns->boff = 0;
			newns->substack = COWVec_new(e->allocator);
		}
		int cmp = NState_cmp(ns, newns, e->gen, nstk);

		if (newly || cmp > 0
		    || (cmp == 0 && e->minvcnt
			&& COWVec_cmp_range(&ns->substack, &newns->substack,
					    e->minvoff, e->minvcnt) > 0))
		{
			COWVec old = newns->substack;

			newns->substack =
			    COWVec_copy((COWVec *) & ns->substack);
			newns->gen = e->gen;
			newns->boff = ns->boff;
			COWVec_free(&old);
		}
		else
		{
			return;
		}
		QSet_insert(e->epsq, k);
	}
}

static void
Execute_add_1arg(Execute *e, QVec *ncsv, NInt k, NInt nstk, const NState *ns,
		 WChar wcnext, NInt arg1)
{
	NState nscopy = *ns;

	nscopy.substack = COWVec_copy((COWVec *) & ns->substack);
	COWVec_put(&nscopy.substack, nstk - 1, arg1);
	Execute_add(e, ncsv, k, nstk, &nscopy, wcnext);
	COWVec_free(&nscopy.substack);
}

static void
Execute_add_2arg(Execute *e, QVec *ncsv, NInt k, NInt nstk, const NState *ns,
		 WChar wcnext, NInt arg1, NInt arg2)
{
	NState nscopy = *ns;

	nscopy.substack = COWVec_copy((COWVec *) & ns->substack);
	COWVec_put(&nscopy.substack, nstk - 2, arg1);
	COWVec_put(&nscopy.substack, nstk - 1, arg2);
	Execute_add(e, ncsv, k, nstk, &nscopy, wcnext);
	COWVec_free(&nscopy.substack);
}

static void
Execute_add_3arg(Execute *e, QVec *ncsv, NInt k, NInt nstk, const NState *ns,
		 WChar wcnext, NInt arg1, NInt arg2, NInt arg3)
{
	NState nscopy = *ns;

	nscopy.substack = COWVec_copy((COWVec *) & ns->substack);
	COWVec_put(&nscopy.substack, nstk - 3, arg1);
	COWVec_put(&nscopy.substack, nstk - 2, arg2);
	COWVec_put(&nscopy.substack, nstk - 1, arg3);
	Execute_add(e, ncsv, k, nstk, &nscopy, wcnext);
	COWVec_free(&nscopy.substack);
}

// Helper for is_word check
static INLINE int
is_word_byte(WChar wc)
{
	return wc == '_' || isalnum(wc);
}

static INLINE int
is_word_wide(WChar wc)
{
	return wc == L'_' || iswalnum(wc);
}

static void
Execute_epsclosure(Execute *e, QVec *ncsv, WChar wcnext)
{
	const Node *nodes = e->nodes;
	int (*is_word)(WChar) =
	    (e->r->enc == WConv_Encoding_Byte) ? is_word_byte : is_word_wide;

	while (!QSet_empty(e->epsq))
	{
		NInt k = QSet_remove(e->epsq);
		NState *ns = &e->epsv->storage[k];

		if (e->has_best
		    && ns->boff > COWVec_get(&e->best, e->suboff + 0))
			continue;

		const Node *n = &nodes[k];
		NInt nstk = n->nstk;

		switch (n->type)
		{
		case Node_Exit:
			{
				size_t b = ns->boff, endoff = e->off;
				NInt mincount =
				    e->r->nmin ? (NInt) COWVec_get(&ns->
								   substack,
								   0) : (NInt)
				    - 1;
				int minvalid =
				    e->r->
				    nmin
				    ? (COWVec_get(&ns->substack, e->minvoff) <
				       ((size_t) 1 << (SIZE_BITS - 1))) : 0;

				if (!e->has_best
				    || b < COWVec_get(&e->best, e->suboff + 0)
				    || (b ==
					COWVec_get(&e->best, e->suboff + 0)
					&& endoff > COWVec_get(&e->best,
							       e->suboff + 1)
					&& (!minvalid
					    || mincount >= e->bestmincount)))
				{
					if (e->has_best)
						COWVec_free(&e->best);
					e->best = COWVec_copy(&ns->substack);
					COWVec_put(&e->best, e->suboff + 0, b);
					COWVec_put(&e->best, e->suboff + 1,
						   endoff);
					e->has_best = 1;
					if (minvalid)
						e->bestmincount =
						    e->bestmincount >
						    mincount ? e->
						    bestmincount : mincount;
				}
			}
			break;

		case Node_Fork:
			{
				NInt priority = (NInt) - 1;
				NInt kk = k;

				do
				{
					Execute_add_1arg(e, ncsv, kk + 1, nstk,
							 ns, wcnext,
							 priority--);
					kk = kk + 1 + nodes[kk].args[0];
				}
				while (nodes[kk].type != Node_Join);
			}
			break;

		case Node_Goto:
			Execute_add(e, ncsv, k + 1 + n->args[1], nstk, ns,
				    wcnext);
			break;

		case Node_Join:
			Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_Loop:
			Execute_add_3arg(e, ncsv, k + 1, nstk, ns, wcnext,
					 (NInt) e->off, (NInt) - 1,
					 (NInt) e->off);
			if (n->args[1])
				Execute_add_3arg(e, ncsv, k + 1 + n->args[0],
						 nstk, ns, wcnext,
						 (NInt) e->off, (NInt) 0,
						 (NInt) e->off);
			break;

		case Node_MinB:
			{
				NState nscopy = *ns;

				nscopy.substack = COWVec_copy(&ns->substack);
				size_t w = n->args[0] / SIZE_BITS;
				size_t b =
				    (size_t) 1 << (SIZE_BITS - 1 -
						   n->args[0] % SIZE_BITS);
				b |= -b;
				size_t x =
				    COWVec_get(&nscopy.substack,
					       e->minvoff + w);
				do
				{
					if ((x & b) != 0)
						COWVec_put(&nscopy.substack,
							   e->minvoff + w,
							   x & ~b);
					b = (size_t) -1;
				}
				while (w-- > 0);
				Execute_add(e, ncsv, k + 1, nstk, &nscopy,
					    wcnext);
				COWVec_free(&nscopy.substack);
			}
			break;

		case Node_MinL:
			Execute_add_1arg(e, ncsv, k + 1, nstk, ns, wcnext,
					 e->off);
			break;

		case Node_MinR:
			{
				NState nscopy = *ns;

				nscopy.substack = COWVec_copy(&ns->substack);
				size_t mininc =
				    e->off - COWVec_get(&nscopy.substack,
							n->nstk);
				size_t oldlen =
				    (size_t) -1 - COWVec_get(&nscopy.substack,
							     n->args[0]);
				mininc -=
				    COWVec_get(&nscopy.substack,
					       e->nestoff + n->args[0]);
				COWVec_put(&nscopy.substack,
					   e->nestoff + n->args[0], 0);
				COWVec_put(&nscopy.substack, n->args[0],
					   (size_t) -1 - (oldlen + mininc));
				for (NInt i = n->args[0]; i-- > 0;)
				{
					oldlen =
					    (size_t) -1 -
					    COWVec_get(&nscopy.substack, i);
					COWVec_put(&nscopy.substack, i,
						   (size_t) -1 - (oldlen +
								  mininc));
					COWVec_put(&nscopy.substack,
						   e->nestoff + i,
						   COWVec_get(&nscopy.substack,
							      e->nestoff + i) +
						   mininc);
				}
				Execute_add(e, ncsv, k + 1, nstk, &nscopy,
					    wcnext);
				COWVec_free(&nscopy.substack);
			}
			break;

		case Node_Next:
			Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			if (n->args[1]
			    && e->off > COWVec_get(&ns->substack,
						   nstk + 3 - 1))
			{
				NState nscopy = *ns;

				nscopy.substack = COWVec_copy(&ns->substack);
				COWVec_put(&nscopy.substack, nstk,
					   COWVec_get(&ns->substack, nstk));
				COWVec_put(&nscopy.substack, nstk + 1,
					   COWVec_get(&ns->substack,
						      nstk + 1) - 1);
				COWVec_put(&nscopy.substack, nstk + 2, e->off);
				Execute_add(e, ncsv, k - n->args[0], nstk + 3,
					    &nscopy, wcnext);
				COWVec_free(&nscopy.substack);
			}
			break;

		case Node_Skip:
			Execute_add_2arg(e, ncsv, k + 1, nstk, ns, wcnext,
					 (NInt) e->off,
					 (NInt) (1 ^ n->args[1]));
			Execute_add_2arg(e, ncsv, k + 1 + n->args[0], nstk, ns,
					 wcnext, (NInt) e->off,
					 (NInt) (0 ^ n->args[1]));
			break;

		case Node_SubL:
			{
				NState nscopy = *ns;

				nscopy.substack = COWVec_copy(&ns->substack);
				COWVec_put(&nscopy.substack, nstk - 1, e->off);
				if (n->args[0] != (NInt) - 1
				    && (e->flags & MINRX_REG_NOSUBRESET) == 0)
				{
					for (NInt i = n->args[0] + 1;
					     i <= n->args[1]; ++i)
					{
						COWVec_put(&nscopy.substack,
							   e->suboff + i * 2,
							   (size_t) -1);
						COWVec_put(&nscopy.substack,
							   e->suboff + i * 2 +
							   1, (size_t) -1);
					}
				}
				Execute_add(e, ncsv, k + 1, nstk, &nscopy,
					    wcnext);
				COWVec_free(&nscopy.substack);
			}
			break;

		case Node_SubR:
			if (n->args[0] != (NInt) - 1
			    && ((e->flags & MINRX_REG_FIRSTSUB) == 0
				|| COWVec_get(&ns->substack,
					      e->suboff + n->args[0] * 2) ==
				(size_t) -1))
			{
				NState nscopy = *ns;

				nscopy.substack = COWVec_copy(&ns->substack);
				COWVec_put(&nscopy.substack,
					   e->suboff + n->args[0] * 2 + 0,
					   COWVec_get(&ns->substack, nstk));
				COWVec_put(&nscopy.substack,
					   e->suboff + n->args[0] * 2 + 1,
					   e->off);
				Execute_add(e, ncsv, k + 1, nstk, &nscopy,
					    wcnext);
				COWVec_free(&nscopy.substack);
			}
			else
			{
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			}
			break;

		case Node_ZBOB:
			if (e->off == 0 && (e->flags & MINRX_REG_NOTBOL) == 0)
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZEOB:
			if (wcnext == WConv_End
			    && (e->flags & MINRX_REG_NOTEOL) == 0)
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZBOL:
			if ((e->off == 0 && (e->flags & MINRX_REG_NOTBOL) == 0)
			    || e->wcprev == L'\n')
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZEOL:
			if ((wcnext == WConv_End
			     && (e->flags & MINRX_REG_NOTEOL) == 0)
			    || wcnext == L'\n')
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZBOW:
			if ((e->off == 0 || !is_word(e->wcprev))
			    && (wcnext != WConv_End && is_word(wcnext)))
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZEOW:
			if ((e->off != 0 && is_word(e->wcprev))
			    && (wcnext == WConv_End || !is_word(wcnext)))
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZXOW:
			if (((e->off == 0 || !is_word(e->wcprev))
			     && (wcnext != WConv_End && is_word(wcnext)))
			    || ((e->off != 0 && is_word(e->wcprev))
				&& (wcnext == WConv_End || !is_word(wcnext))))
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		case Node_ZNWB:
			if ((e->off == 0 && wcnext == WConv_End)
			    || (e->off == 0 && wcnext != WConv_End
				&& !is_word(wcnext)) || (e->off != 0
							 && !is_word(e->wcprev)
							 && wcnext ==
							 WConv_End)
			    || (e->off != 0 && wcnext != WConv_End
				&& is_word(e->wcprev) == is_word(wcnext)))
				Execute_add(e, ncsv, k + 1, nstk, ns, wcnext);
			break;

		default:
			abort();
			break;
		}
	}
}

static int
Execute_execute(Execute *e, size_t nm, minrx_regmatch_t *rm)
{
	QVec *mcsvs[2];

	mcsvs[0] = QVec_new(e->r->nodes->size);
	mcsvs[1] = QVec_new(e->r->nodes->size);

	e->off = WConv_off(&e->wconv);
	WChar wcnext = WConv_nextchr(&e->wconv);

	if ((e->flags & MINRX_REG_RESUME) != 0 && rm && rm[0].rm_eo > 0)
	{
		while (wcnext != WConv_End && (ptrdiff_t) e->off < rm[0].rm_eo)
		{
			e->wcprev = wcnext;
			e->off = WConv_off(&e->wconv);
			wcnext = WConv_nextchr(&e->wconv);
		}
	}

	NState nsinit;

	nsinit.gen = 0;
	nsinit.boff = 0;
	nsinit.substack = COWVec_new(e->allocator);

	if ((e->flags & MINRX_REG_NOFIRSTBYTES) == 0 && e->r->has_firstbytes
	    && !CSet_test(e->r->firstcset, wcnext))
	{
	      zoom:
		{
			const char *cp = e->wconv.cp;
			const char *ep = e->wconv.ep;

			if (e->r->has_firstunique)
			{
				cp = (const char *) memchr(cp,
							   e->r->firstunique,
							   ep - cp);
				if (cp == NULL)
					goto exit_label;
			}
			else
			{
				while (cp != ep
				       && !e->r->
				       firstbytes[(unsigned char) *cp])
					++cp;
				if (cp == ep)
					goto exit_label;
			}
			if (cp != e->wconv.cp)
			{
				if (e->r->enc == WConv_Encoding_UTF8)
				{
					const char *bp = cp;

					while (bp != e->wconv.cp && cp - bp < 8
					       && (unsigned char) *--bp >=
					       0x80)
						;
					e->wconv.cp =
					    (unsigned char) *bp >=
					    0x80 ? cp - 1 : bp;
				}
				else
				{
					e->wconv.cp = cp - 1;
				}
				wcnext = WConv_nextchr(&e->wconv);
			}
			++e->gen;
			e->wcprev = wcnext;
			e->off = WConv_off(&e->wconv);
			wcnext = WConv_nextchr(&e->wconv);
		}
	}

	nsinit.boff = e->off;
	for (size_t i = 0; i < e->r->nmin; ++i)
		COWVec_put(&nsinit.substack, e->nestoff + i, 0);

	Execute_add(e, mcsvs[0], 0, 0, &nsinit, wcnext);
	if (!QSet_empty(e->epsq))
		Execute_epsclosure(e, mcsvs[0], wcnext);

	for (;;)
	{
		if (wcnext == WConv_End)
			break;
		++e->gen;
		e->wcprev = wcnext;
		e->off = WConv_off(&e->wconv);
		wcnext = WConv_nextchr(&e->wconv);

		while (!QVec_empty(mcsvs[0]))
		{
			size_t n;
			NState ns;
			int newly __attribute__((unused)) = 0;

			// Remove from queue
			n = QSet_remove(mcsvs[0]->qset);
			ns = mcsvs[0]->storage[n];
			Execute_add(e, mcsvs[1], n + 1, e->nodes[n].nstk, &ns,
				    wcnext);
			COWVec_free(&ns.substack);
		}

		if (!e->has_best)
		{
			nsinit.boff = e->off;
			Execute_add(e, mcsvs[1], 0, 0, &nsinit, wcnext);
		}

		if (!QSet_empty(e->epsq))
			Execute_epsclosure(e, mcsvs[1], wcnext);

		if (QVec_empty(mcsvs[1]))
		{
			if (e->has_best)
				break;
			if ((e->flags & MINRX_REG_NOFIRSTBYTES) == 0
			    && e->r->has_firstbytes)
				goto zoom;
		}

		if (wcnext == WConv_End)
			break;
		++e->gen;
		e->wcprev = wcnext;
		e->off = WConv_off(&e->wconv);
		wcnext = WConv_nextchr(&e->wconv);

		while (!QVec_empty(mcsvs[1]))
		{
			size_t n;
			NState ns;

			n = QSet_remove(mcsvs[1]->qset);
			ns = mcsvs[1]->storage[n];
			Execute_add(e, mcsvs[0], n + 1, e->nodes[n].nstk, &ns,
				    wcnext);
			COWVec_free(&ns.substack);
		}

		if (!e->has_best)
		{
			nsinit.boff = e->off;
			Execute_add(e, mcsvs[0], 0, 0, &nsinit, wcnext);
		}

		if (!QSet_empty(e->epsq))
			Execute_epsclosure(e, mcsvs[0], wcnext);

		if (QVec_empty(mcsvs[0]))
		{
			if (e->has_best)
				break;
			if ((e->flags & MINRX_REG_NOFIRSTBYTES) == 0
			    && e->r->has_firstbytes)
				goto zoom;
		}
	}

      exit_label:
	COWVec_free(&nsinit.substack);
	QVec_free(mcsvs[0]);
	QVec_free(mcsvs[1]);

	if (e->has_best)
	{
		if (rm)
		{
			size_t nsub = nm < e->r->nsub ? nm : e->r->nsub;
			size_t i;

			for (i = 0; i < nsub; ++i)
			{
				rm[i].rm_so =
				    (ptrdiff_t) COWVec_get(&e->best,
							   e->suboff + i * 2);
				rm[i].rm_eo =
				    (ptrdiff_t) COWVec_get(&e->best,
							   e->suboff + i * 2 +
							   1);
			}
			for (; i < nm; ++i)
				rm[i].rm_so = rm[i].rm_eo = -1;
		}
		return 0;
	}
	else
	{
		if (rm)
		{
			for (size_t i = 0; i < nm; ++i)
				rm[i].rm_so = rm[i].rm_eo = -1;
		}
		return MINRX_REG_NOMATCH;
	}
}

// Public API functions
int
minrx_regcomp(minrx_regex_t *rx, const char *s, int flags)
{
	return minrx_regncomp(rx, strlen(s), s, flags);
}

int
minrx_regexec(minrx_regex_t *rx, const char *s, size_t nm,
	      minrx_regmatch_t *rm, int flags)
{
	return minrx_regnexec(rx, strlen(s), s, nm, rm, flags);
}

int
minrx_regncomp(minrx_regex_t *rx, size_t ns, const char *s, int flags)
{
	WConv_Encoding enc = WConv_Encoding_MBtoWC;
	const char *loc = setlocale(LC_CTYPE, NULL);

	if ((strcmp(loc, "C") == 0 || (flags & MINRX_REG_NATIVE1B) != 0)
	    && MB_CUR_MAX == 1)
		enc = WConv_Encoding_Byte;
	else if (strcmp(nl_langinfo(CODESET), "UTF-8") == 0)
		enc = WConv_Encoding_UTF8;

	Compile c;

	c.flags = (minrx_regcomp_flags_t) flags;
	c.enc = enc;
	WConv_init(&c.wconv, enc, s, s + ns);
	c.wc = WConv_nextchr(&c.wconv);
	c.csets = CSetArray_new();
	c.has_dot = 0;
	c.has_esc_s = 0;
	c.has_esc_S = 0;
	c.has_esc_w = 0;
	c.has_esc_W = 0;
	c.icmap = NULL;
	c.nmin = 0;
	c.nsub = 0;

	Regexp *r = Compile_compile(&c);

	if (!r)
	{
		r = Regexp_new();
		r->err = MINRX_REG_ESPACE;
	}

	// Clean up icmap and csets
	Icmap_free(c.icmap);
	CSetArray_free(c.csets);

	rx->re_regexp = r;
	rx->re_nsub = r->nsub > 0 ? r->nsub - 1 : 0;
	rx->re_compflags = (minrx_regcomp_flags_t) flags;
	return r->err;
}

int
minrx_regnexec(minrx_regex_t *rx, size_t ns, const char *s, size_t nm,
	       minrx_regmatch_t *rm, int flags)
{
	Regexp *r = (Regexp *) rx->re_regexp;
	Execute e;

	e.r = r;
	e.flags = (minrx_regexec_flags_t) flags;
	e.suboff = r->nmin + r->nstk;
	e.minvcnt = (r->nmin + SIZE_BITS - 1) / SIZE_BITS;
	e.minvoff = e.suboff + 2 * r->nsub;
	e.nestoff = e.minvoff + e.minvcnt;
	e.gen = 0;
	e.off = 0;
	WConv_init(&e.wconv, r->enc, s, s + ns);
	e.wcprev = WConv_End;
	e.allocator = COWVec_Allocator_new(e.nestoff + r->nmin);
	e.has_best = 0;
	e.bestmincount = 0;
	e.epsv = QVec_new(r->nodes->size);
	e.epsq = QSet_new(r->nodes->size);
	e.nodes = r->nodes->data;

	int result = Execute_execute(&e, nm, rm);

	// Clean up
	if (e.has_best)
		COWVec_free(&e.best);
	QVec_free(e.epsv);
	QSet_free(e.epsq);
	COWVec_Allocator_free(e.allocator);

	return result;
}

void
minrx_regfree(minrx_regex_t *rx)
{
	Regexp_free((Regexp *) rx->re_regexp);
	rx->re_regexp = NULL;
}

size_t
minrx_regerror(int errcode, const minrx_regex_t *rx
	       __attribute__((unused)), char *errbuf, size_t errsize)
{
	static const char *const messages[] = {
		N_("success"),
		N_("bad pattern"),
		N_("invalid contents of {}"),
		N_("? * + or {interval} not preceded by valid subpattern"),
		N_("unbalanced {"),
		N_("unbalanced ["),
		N_("invalid collating element"),
		N_("invalid character class name"),
		N_("invalid trailing backslash"),
		N_("unbalanced ("),
		N_("invalid range endpoint"),
		N_("memory allocation failed"),
		N_("invalid \\digit"),
		N_("match not found"),
		N_("unknown error code"),
	};
	if (errcode < 0 || errcode > MINRX_REG_UNKNOWN)
		errcode = MINRX_REG_UNKNOWN;
	size_t size = snprintf(errbuf, errsize, "%s", _(messages[errcode]));

	if (errsize != 0 && size >= errsize)
		errbuf[errsize - 1] = '\0';
	return size + 1;
}
