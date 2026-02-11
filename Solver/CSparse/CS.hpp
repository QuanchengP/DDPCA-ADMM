#ifndef _CS_hpp
#define _CS_hpp

#include "../../General/General.hpp"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstddef>
#include <cstring>

namespace Ddpca{

// #define csi ptrdiff_t //substituted by I64

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct CSSparse    /* matrix in compressed-column or triplet form */
{
    I64 nzmax ;     /* maximum number of entries */
    I64 m ;         /* number of rows */
    I64 n ;         /* number of columns */
    I64 *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    I64 *i ;        /* row indices, size nzmax */
    Real *x ;     /* numerical values, size nzmax */
    I64 nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} CS ;

CS *CSAdd (const CS *A, const CS *B, Real alpha, Real beta) ;
I64 CSCholsol (I64 order, const CS *A, Real *b) ;
CS *CSCompress (const CS *T) ;
I64 CSDupl (CS *A) ;
I64 CSEntry (CS *T, I64 i, I64 j, Real x) ;
I64 CSGaxpy (const CS *A, const Real *x, Real *y) ;
CS *CSLoad (FILE *f) ;
I64 CSLusol (I64 order, const CS *A, Real *b, Real tol) ;
CS *CSMultiply (const CS *A, const CS *B) ;
Real CSNorm (const CS *A) ;
I64 CSPrint (const CS *A, I64 brief) ;
I64 CSQrsol (I64 order, const CS *A, Real *b) ;
CS *CSTranspose (const CS *A, I64 values) ;
/* utilities */
void *CSMalloc (I64 n, I64 size) ;
void *CSCalloc (I64 n, I64 size) ;
void *CSFree (void *p) ;
void *CSRealloc (void *p, I64 n, I64 size, I64 *ok, I64 oldN) ;
CS *CSSpalloc (I64 m, I64 n, I64 nzmax, I64 values, I64 triplet) ;
CS *CSSpfree (CS *A) ;
I64 CSSprealloc (CS *A, I64 nzmax) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct CSSymbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    I64 *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    I64 *q ;        /* fill-reducing column permutation for LU and QR */
    I64 *parent ;   /* elimination tree for Cholesky and QR */
    I64 *cp ;       /* column pointers for Cholesky, row counts for QR */
    I64 *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    I64 m2 ;        /* # of rows for QR, after adding fictitious rows */
    Real lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    Real unz ;    /* # entries in U for LU; in R for QR */
} css ;

typedef struct CSNumeric   /* numeric Cholesky, LU, or QR factorization */
{
    CS *L ;         /* L for LU and Cholesky, V for QR */
    CS *U ;         /* U for LU, R for QR, not used for Cholesky */
    I64 *pinv ;     /* partial pivoting for LU */
    Real *B ;     /* beta [0..n-1] for QR */
} CSN ;

typedef struct CSDmpermResults    /* CSDmperm or CSScc output */
{
    I64 *p ;        /* size m, row permutation */
    I64 *q ;        /* size n, column permutation */
    I64 *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    I64 *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    I64 nb ;        /* # of blocks in fine dmperm decomposition */
    I64 rr [5] ;    /* coarse row decomposition */
    I64 cc [5] ;    /* coarse column decomposition */
} CSD ;

I64 *CSAmd (I64 order, const CS *A) ;
CSN *CSChol (const CS *A, const css *S) ;
CSD *CSDmperm (const CS *A, I64 seed) ;
I64 CSDroptol (CS *A, Real tol) ;
I64 CSDropzeros (CS *A) ;
I64 CSHapply (const CS *V, I64 i, Real beta, Real *x) ;
I64 CSIpvec (const I64 *p, const Real *b, Real *x, I64 n) ;
I64 CSLsolve (const CS *L, Real *x) ;
I64 CSLtsolve (const CS *L, Real *x) ;
CSN *CSLu (const CS *A, const css *S, Real tol) ;
CS *CSPermute (const CS *A, const I64 *pinv, const I64 *q, I64 values) ;
I64 *CSPinv (const I64 *p, I64 n) ;
I64 CSPvec (const I64 *p, const Real *b, Real *x, I64 n) ;
CSN *CSQr (const CS *A, const css *S) ;
css *CSSchol (I64 order, const CS *A) ;
css *CSSqr (I64 order, const CS *A, I64 qr) ;
CS *CSSymperm (const CS *A, const I64 *pinv, I64 values) ;
I64 CSUpdown (CS *L, I64 sigma, const CS *C, const I64 *parent) ;
I64 CSUsolve (const CS *U, Real *x) ;
I64 CSUtsolve (const CS *U, Real *x) ;
/* utilities */
css *CSSfree (css *S) ;
CSN *CSNfree (CSN *N) ;
CSD *CSDfree (CSD *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
I64 *CSCounts (const CS *A, const I64 *parent, const I64 *post, I64 ata) ;
Real CSCumsum (I64 *p, I64 *c, I64 n) ;
I64 CSDfs (I64 j, CS *G, I64 top, I64 *xi, I64 *pstack, const I64 *pinv) ;
I64 CSEreach (const CS *A, I64 k, const I64 *parent, I64 *s, I64 *w) ;
I64 *CSEtree (const CS *A, I64 ata) ;
I64 CSFkeep (CS *A, I64 (*fkeep) (I64, I64, Real, void *), void *other) ;
Real CSHouse (Real *x, Real *beta, I64 n) ;
I64 CSLeaf (I64 i, I64 j, const I64 *first, I64 *maxfirst, I64 *prevleaf,
    I64 *ancestor, I64 *jleaf) ;
I64 *CSMaxtrans (const CS *A, I64 seed) ;
I64 *CSPost (const I64 *parent, I64 n) ;
I64 *CSRandperm (I64 n, I64 seed) ;
I64 CSReach (CS *G, const CS *B, I64 k, I64 *xi, const I64 *pinv) ;
I64 CSScatter (const CS *A, I64 j, Real beta, I64 *w, Real *x, I64 mark,
    CS *C, I64 nz) ;
CSD *CSScc (CS *A) ;
I64 CSSpsolve (CS *G, const CS *B, I64 k, I64 *xi, Real *x,
    const I64 *pinv, I64 lo) ;
I64 CSTdfs (I64 j, I64 k, I64 *head, const I64 *next, I64 *post,
    I64 *stack) ;
// /* utilities */
CSD *CSDalloc (I64 m, I64 n) ;
CSD *CSDdone (CSD *D, CS *C, void *w, I64 ok) ;
CS *CSDone (CS *C, void *w, void *x, I64 ok) ;
I64 *CSIdone (I64 *p, CS *C, void *w, I64 ok) ;
CSN *CSNdone (CSN *N, CS *C, void *w, void *x, I64 ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

}

#endif // _CS_hpp