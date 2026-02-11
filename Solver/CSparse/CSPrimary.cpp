
#include "CS.hpp"

#define CS_VER 3                    /* CSparse Version */
#define CS_SUBVER 1
#define CS_SUBSUB 4
#define CS_DATE "Oct 10, 2014"    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2014"

namespace Ddpca{

/* C = alpha*A + beta*B */
CS *CSAdd (const CS *A, const CS *B, Real alpha, Real beta)
{
    I64 p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
    Real *x, *Bx, *Cx ;
    CS *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (nullptr) ;         /* check inputs */
    if (A->m != B->m || A->n != B->n) return (nullptr) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
    w = (I64 *)CSCalloc (m, sizeof (I64)) ;                       /* get workspace */
    values = (A->x != nullptr) && (Bx != nullptr) ;
    x = values ? (Real *)CSMalloc (m, sizeof (Real)) : nullptr ;    /* get workspace */
    C = CSSpalloc (m, n, anz + bnz, values, 0) ;           /* allocate result*/
    if (!C || !w || (values && !x)) return (CSDone (C, w, x, 0)) ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)
    {
        Cp [j] = nz ;                   /* column j of C starts here */
        nz = CSScatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = CSScatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    CSSprealloc (C, 0) ;               /* remove extra space from C */
    return (CSDone (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
I64 CSCholsol (I64 order, const CS *A, Real *b)
{
    Real *x ;
    css *S ;
    CSN *N ;
    I64 n, ok ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;
    S = CSSchol (order, A) ;               /* ordering and symbolic analysis */
    N = CSChol (A, S) ;                    /* numeric Cholesky factorization */
    x = (Real *)CSMalloc (n, sizeof (Real)) ;    /* get workspace */
    ok = (S && N && x) ;
    if (ok)
    {
        CSIpvec (S->pinv, b, x, n) ;   /* x = P*b */
        CSLsolve (N->L, x) ;           /* x = L\x */
        CSLtsolve (N->L, x) ;          /* x = L'\x */
        CSPvec (S->pinv, x, b, n) ;    /* b = P'*x */
    }
    CSFree (x) ;
    CSSfree (S) ;
    CSNfree (N) ;
    return (ok) ;
}

/* C = compressed-column form of a triplet matrix T */
CS *CSCompress (const CS *T)
{
    I64 m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    Real *Cx, *Tx ;
    CS *C ;
    if (!CS_TRIPLET (T)) return (nullptr) ;                /* check inputs */
    m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
    C = CSSpalloc (m, n, nz, Tx != nullptr, 0) ;          /* allocate result */
    w = (I64*)CSCalloc (n, sizeof (I64)) ;             /* get workspace */
    if (!C || !w) return (CSDone (C, w, nullptr, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    CSCumsum (Cp, w, n) ;                              /* column pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (CSDone (C, w, nullptr, 1)) ;      /* success; free w and return C */
}

/* remove duplicate entries from A */
I64 CSDupl (CS *A)
{
    I64 i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
    Real *Ax ;
    if (!CS_CSC (A)) return (0) ;               /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    w = (I64*)CSMalloc (m, sizeof (I64)) ;     /* get workspace */
    if (!w) return (0) ;                        /* out of memory */
    for (i = 0 ; i < m ; i++) w [i] = -1 ;      /* row i not yet seen */
    for (j = 0 ; j < n ; j++)
    {
        q = nz ;                                /* column j will start at q */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;                        /* A(i,j) is nonzero */
            if (w [i] >= q)
            {
                Ax [w [i]] += Ax [p] ;          /* A(i,j) is a duplicate */
            }
            else
            {
                w [i] = nz ;                    /* record where row i occurs */
                Ai [nz] = i ;                   /* keep A(i,j) */
                Ax [nz++] = Ax [p] ;
            }
        }
        Ap [j] = q ;                            /* record start of column j */
    }
    Ap [n] = nz ;                               /* finalize A */
    CSFree (w) ;                               /* free workspace */
    return (CSSprealloc (A, 0)) ;              /* remove extra space from A */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
I64 CSEntry (CS *T, I64 i, I64 j, Real x)
{
    if (!CS_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
    if (T->nz >= T->nzmax && !CSSprealloc (T,2*(T->nzmax))) return (0) ;
    if (T->x) T->x [T->nz] = x ;
    T->i [T->nz] = i ;
    T->p [T->nz++] = j ;
    T->m = CS_MAX (T->m, i+1) ;
    T->n = CS_MAX (T->n, j+1) ;
    return (1) ;
}

/* y = A*x+y */
I64 CSGaxpy (const CS *A, const Real *x, Real *y)
{
    I64 p, j, n, *Ap, *Ai ;
    Real *Ax ;
    if (!CS_CSC (A) || !x || !y) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            y [Ai [p]] += Ax [p] * x [j] ;
        }
    }
    return (1) ;
}

/* load a triplet matrix from a file */
CS *CSLoad (FILE *f)
{
    Real i, j ;   /* use Real for integers to avoid I64 conflicts */
    Real x ;
    CS *T ;
    if (!f) return (nullptr) ;                             /* check inputs */
    T = CSSpalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
    while (fscanf (f, "%lg %lg %lg\n", &i, &j, &x) == 3)
    {
        if (!CSEntry (T, (I64) i, (I64) j, x)) return (CSSpfree (T)) ;
    }
    return (T) ;
}

/* x=A\b where A is unsymmetric; b overwritten with solution */
I64 CSLusol (I64 order, const CS *A, Real *b, Real tol)
{
    Real *x ;
    css *S ;
    CSN *N ;
    I64 n, ok ;
    if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
    n = A->n ;
    S = CSSqr (order, A, 0) ;              /* ordering and symbolic analysis */
    N = CSLu (A, S, tol) ;                 /* numeric LU factorization */
    x = (Real *)CSMalloc (n, sizeof (Real)) ;    /* get workspace */
    ok = (S && N && x) ;
    if (ok)
    {
        CSIpvec (N->pinv, b, x, n) ;       /* x = b(p) */
        CSLsolve (N->L, x) ;               /* x = L\x */
        CSUsolve (N->U, x) ;               /* x = U\x */
        CSIpvec (S->q, x, b, n) ;          /* b(q) = x */
    }
    CSFree (x) ;
    CSSfree (S) ;
    CSNfree (N) ;
    return (ok) ;
}

/* C = A*B */
CS *CSMultiply (const CS *A, const CS *B)
{
    I64 p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    Real *x, *Bx, *Cx ;
    CS *C ;
    if (!CS_CSC (A) || !CS_CSC (B)) return (nullptr) ;      /* check inputs */
    if (A->n != B->m) return (nullptr) ;
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    w = (I64 *)CSCalloc (m, sizeof (I64)) ;                    /* get workspace */
    values = (A->x != nullptr) && (Bx != nullptr) ;
    x = values ? (Real *)CSMalloc (m, sizeof (Real)) : nullptr ; /* get workspace */
    C = CSSpalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
    if (!C || !w || (values && !x)) return (CSDone (C, w, x, 0)) ;
    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !CSSprealloc (C, 2*(C->nzmax)+m))
        {
            return (CSDone (C, w, x, 0)) ;             /* out of memory */
        } 
        Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = CSScatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    CSSprealloc (C, 0) ;               /* remove extra space from C */
    return (CSDone (C, w, x, 1)) ;     /* success; free workspace, return C */
}

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
Real CSNorm (const CS *A)
{
    I64 p, j, n, *Ap ;
    Real *Ax,  norm = 0, s ;
    if (!CS_CSC (A) || !A->x) return (-1) ;             /* check inputs */
    n = A->n ; Ap = A->p ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]) ;
        norm = CS_MAX (norm, s) ;
    }
    return (norm) ;
}

/* print a sparse matrix; use %g for integers to avoid differences with I64 */
I64 CSPrint (const CS *A, I64 brief)
{
    I64 p, j, m, n, nzmax, nz, *Ap, *Ai ;
    Real *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    printf ("CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (Real) m,
            (Real) n, (Real) nzmax, (Real) (Ap [n]), CSNorm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col %g : locations %g to %g\n", (Real) j, 
                (Real) (Ap [j]), (Real) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : %g\n", (Real) (Ai [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (Real) m,
            (Real) n, (Real) nzmax, (Real) nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            printf ("    %g %g : %g\n", (Real) (Ai [p]), (Real) (Ap [p]),
                Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}

/* x=A\b where A can be rectangular; b overwritten with solution */
I64 CSQrsol (I64 order, const CS *A, Real *b)
{
    Real *x ;
    css *S ;
    CSN *N ;
    CS *AT = nullptr ;
    I64 k, m, n, ok ;
    if (!CS_CSC (A) || !b) return (0) ; /* check inputs */
    n = A->n ;
    m = A->m ;
    if (m >= n)
    {
        S = CSSqr (order, A, 1) ;          /* ordering and symbolic analysis */
        N = CSQr (A, S) ;                  /* numeric QR factorization */
        x = (Real *)CSCalloc (S ? S->m2 : 1, sizeof (Real)) ;    /* get workspace */
        ok = (S && N && x) ;
        if (ok)
        {
            CSIpvec (S->pinv, b, x, m) ;   /* x(0:m-1) = b(p(0:m-1) */
            for (k = 0 ; k < n ; k++)       /* apply Householder refl. to x */
            {
                CSHapply (N->L, k, N->B [k], x) ;
            }
            CSUsolve (N->U, x) ;           /* x = R\x */
            CSIpvec (S->q, x, b, n) ;      /* b(q(0:n-1)) = x(0:n-1) */
        }
    }
    else
    {
        AT = CSTranspose (A, 1) ;          /* Ax=b is underdetermined */
        S = CSSqr (order, AT, 1) ;         /* ordering and symbolic analysis */
        N = CSQr (AT, S) ;                 /* numeric QR factorization of A' */
        x = (Real *)CSCalloc (S ? S->m2 : 1, sizeof (Real)) ;    /* get workspace */
        ok = (AT && S && N && x) ;
        if (ok)
        {
            CSPvec (S->q, b, x, m) ;       /* x(q(0:m-1)) = b(0:m-1) */
            CSUtsolve (N->U, x) ;          /* x = R'\x */
            for (k = m-1 ; k >= 0 ; k--)    /* apply Householder refl. to x */
            {
                CSHapply (N->L, k, N->B [k], x) ;
            }
            CSPvec (S->pinv, x, b, n) ;    /* b(0:n-1) = x(p(0:n-1)) */
        }
    }
    CSFree (x) ;
    CSSfree (S) ;
    CSNfree (N) ;
    CSSpfree (AT) ;
    return (ok) ;
}

/* C = A' */
CS *CSTranspose (const CS *A, I64 values)
{
    I64 p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    Real *Cx, *Ax ;
    CS *C ;
    if (!CS_CSC (A)) return (nullptr) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = CSSpalloc (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    w = (I64*)CSCalloc (m, sizeof (I64)) ;                /* get workspace */
    if (!C || !w) return (CSDone (C, w, nullptr, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    CSCumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (CSDone (C, w, nullptr, 1)) ;  /* success; free w and return C */
}

/* wrapper for malloc */
void *CSMalloc (I64 n, I64 size)
{
    // return (malloc (CS_MAX (n,1) * size)) ;
    //********************20260102********************
    const I64 totalSize = CS_MAX (n,1) * size;
    const I64 alignedSize = ((totalSize + nfsAlign - 1) / nfsAlign) * nfsAlign;
    return std::aligned_alloc(nfsAlign, alignedSize);
    //********************20260102********************
}

/* wrapper for calloc */
void *CSCalloc (I64 n, I64 size)
{
    // return (calloc (CS_MAX (n,1), size)) ;
    //********************20260102********************
    const I64 totalSize = CS_MAX (n,1) * size;
    const I64 alignedSize = ((totalSize + nfsAlign - 1) / nfsAlign) * nfsAlign;
    void* ptr = std::aligned_alloc(nfsAlign, alignedSize);
    std::memset(ptr, 0, alignedSize);
    return ptr;
    //********************20260102********************
}

/* wrapper for free */
void *CSFree (void *p)
{
    if (p) std::free (p) ;       /* free p if it is not already nullptr */
    return (nullptr) ;         /* return nullptr to simplify the use of CSFree */
}

/* wrapper for realloc */
void *CSRealloc (void *p, I64 n, I64 size, I64 *ok, I64 oldN)
{
    // void *pnew ;
    // pnew = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
    // *ok = (pnew != nullptr) ;                  /* realloc fails if pnew is nullptr */
    // return pnew;//((*ok) ? pnew : p) ;       /* return original p if failure */
    //********************20260102********************
    const I64 totalSize = CS_MAX (n,1) * size;
    const I64 alignedSize = ((totalSize + nfsAlign - 1) / nfsAlign) * nfsAlign;
    void* pnew = std::aligned_alloc(nfsAlign, alignedSize);
    std::memcpy(pnew, p, oldN * size);
    std::free(p);
    *ok = (pnew != nullptr) ;
    return pnew;
    //********************20260102********************
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
CS *CSSpalloc (I64 m, I64 n, I64 nzmax, I64 values, I64 triplet)
{
    CS *A = (CS *)CSCalloc (1, sizeof (CS)) ;    /* allocate the CS struct */
    if (!A) return (nullptr) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = (I64 *)CSMalloc (triplet ? nzmax : n+1, sizeof (I64)) ;
    A->i = (I64 *)CSMalloc (nzmax, sizeof (I64)) ;
    A->x = values ? (Real *)CSMalloc (nzmax, sizeof (Real)) : nullptr ;
    return ((!A->p || !A->i || (values && !A->x)) ? CSSpfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
I64 CSSprealloc (CS *A, I64 nzmax)
{
    I64 ok, oki, okj = 1, okx = 1 ;
    if (!A) return (0) ;
    //********************20260102********************
    I64 nzmax_0 = (CS_CSC (A)) ? (A->p [A->n]) : A->nz ;//20260102
    if (nzmax <= 0) nzmax = nzmax_0;
    A->i = (I64 *)CSRealloc (A->i, nzmax, sizeof (I64), &oki, nzmax_0) ;
    if (CS_TRIPLET (A)) A->p = (I64 *)CSRealloc (A->p, nzmax, sizeof (I64), &okj, nzmax_0) ;
    if (A->x) A->x = (Real *)CSRealloc (A->x, nzmax, sizeof (Real), &okx, nzmax_0) ;
    //********************20260102********************
    ok = (oki && okj && okx) ;
    if (ok) A->nzmax = nzmax ;
    return (ok) ;
}

/* free a sparse matrix */
CS *CSSpfree (CS *A)
{
    if (!A) return (nullptr) ;     /* do nothing if A already nullptr */
    CSFree (A->p) ;
    CSFree (A->i) ;
    CSFree (A->x) ;
    return ((CS *) CSFree (A)) ;   /* free the CS struct and return nullptr */
}

}
