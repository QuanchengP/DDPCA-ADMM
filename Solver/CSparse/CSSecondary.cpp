
#include "CS.hpp"

namespace Ddpca{

/* clear w */
static I64 cs_wclear (I64 mark, I64 lemax, I64 *w, I64 n)
{
    I64 k ;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
    }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}

/* keep off-diagonal entries; drop diagonal entries */
static I64 cs_diag (I64 i, I64 j, Real, void *) { return (i != j) ; }

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
I64 *CSAmd (I64 order, const CS *A)  /* order 0:natural, 1:Chol, 2:LU, 3:QR */
{
    CS *C, *A2, *AT ;
    I64 *Cp, *Ci, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
        *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
        k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
        ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
    I64 h ;
    /* --- Construct matrix C ----------------------------------------------- */
    if (!CS_CSC (A) || order <= 0 || order > 3) return (nullptr) ; /* check */
    AT = CSTranspose (A, 0) ;              /* compute A' */
    if (!AT) return (nullptr) ;
    m = A->m ; n = A->n ;
    dense = CS_MAX (16, 10 * sqrt ((Real) n)) ;   /* find dense threshold */
    dense = CS_MIN (n-2, dense) ;
    if (order == 1 && n == m)
    {
        C = CSAdd (A, AT, 0, 0) ;          /* C = A+A' */
    }
    else if (order == 2)
    {
        ATp = AT->p ;                       /* drop dense columns from AT */
        ATi = AT->i ;
        for (p2 = 0, j = 0 ; j < m ; j++)
        {
            p = ATp [j] ;                   /* column j of AT starts here */
            ATp [j] = p2 ;                  /* new column j starts here */
            if (ATp [j+1] - p > dense) continue ;   /* skip dense col j */
            for ( ; p < ATp [j+1] ; p++) ATi [p2++] = ATi [p] ;
        }
        ATp [m] = p2 ;                      /* finalize AT */
        A2 = CSTranspose (AT, 0) ;         /* A2 = AT' */
        C = A2 ? CSMultiply (AT, A2) : nullptr ;  /* C=A'*A with no dense rows */
        CSSpfree (A2) ;
    }
    else
    {
        C = CSMultiply (AT, A) ;           /* C=A'*A */
    }
    CSSpfree (AT) ;
    if (!C) return (nullptr) ;
    CSFkeep (C, &cs_diag, nullptr) ;          /* drop diagonal entries */
    Cp = C->p ;
    cnz = Cp [n] ;
    P = (I64 *)CSMalloc (n+1, sizeof (I64)) ;     /* allocate result */
    W = (I64 *)CSMalloc (8*(n+1), sizeof (I64)) ; /* get workspace */
    t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
    if (!P || !W || !CSSprealloc (C, t)) return (CSIdone (P, C, W, 0)) ;
    len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
    head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
    w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
    last = P ;                              /* use P as workspace for last */
    /* --- Initialize quotient graph ---------------------------------------- */
    for (k = 0 ; k < n ; k++) len [k] = Cp [k+1] - Cp [k] ;
    len [n] = 0 ;
    nzmax = C->nzmax ;
    Ci = C->i ;
    for (i = 0 ; i <= n ; i++)
    {
        head [i] = -1 ;                     /* degree list i is empty */
        last [i] = -1 ;
        next [i] = -1 ;
        hhead [i] = -1 ;                    /* hash list i is empty */
        nv [i] = 1 ;                        /* node i is just one node */
        w [i] = 1 ;                         /* node i is alive */
        elen [i] = 0 ;                      /* Ek of node i is empty */
        degree [i] = len [i] ;              /* degree of node i */
    }
    mark = cs_wclear (0, 0, w, n) ;         /* clear w */
    elen [n] = -2 ;                         /* n is a dead element */
    Cp [n] = -1 ;                           /* n is a root of assembly tree */
    w [n] = 0 ;                             /* n is a dead element */
    /* --- Initialize degree lists ------------------------------------------ */
    for (i = 0 ; i < n ; i++)
    {
        d = degree [i] ;
        if (d == 0)                         /* node i is empty */
        {
            elen [i] = -2 ;                 /* element i is dead */
            nel++ ;
            Cp [i] = -1 ;                   /* i is a root of assembly tree */
            w [i] = 0 ;
        }
        else if (d > dense)                 /* node i is dense */
        {
            nv [i] = 0 ;                    /* absorb i into element n */
            elen [i] = -1 ;                 /* node i is dead */
            nel++ ;
            Cp [i] = CS_FLIP (n) ;
            nv [n]++ ;
        }
        else
        {
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;           /* put node i in degree list d */
            head [d] = i ;
        }
    }
    while (nel < n)                         /* while (selecting pivots) do */
    {
        /* --- Select node of minimum approximate degree -------------------- */
        for (k = -1 ; mindeg < n && (k = head [mindeg]) == -1 ; mindeg++) ;
        if (next [k] != -1) last [next [k]] = -1 ;
        head [mindeg] = next [k] ;          /* remove k from degree list */
        elenk = elen [k] ;                  /* elenk = |Ek| */
        nvk = nv [k] ;                      /* # of nodes k represents */
        nel += nvk ;                        /* nv[k] nodes of A eliminated */
        /* --- Garbage collection ------------------------------------------- */
        if (elenk > 0 && cnz + mindeg >= nzmax)
        {
            for (j = 0 ; j < n ; j++)
            {
                if ((p = Cp [j]) >= 0)      /* j is a live node or element */
                {
                    Cp [j] = Ci [p] ;       /* save first entry of object */
                    Ci [p] = CS_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
                }
            }
            for (q = 0, p = 0 ; p < cnz ; ) /* scan all of memory */
            {
                if ((j = CS_FLIP (Ci [p++])) >= 0)  /* found object j */
                {
                    Ci [q] = Cp [j] ;       /* restore first entry of object */
                    Cp [j] = q++ ;          /* new pointer to object j */
                    for (k3 = 0 ; k3 < len [j]-1 ; k3++) Ci [q++] = Ci [p++] ;
                }
            }
            cnz = q ;                       /* Ci [cnz...nzmax-1] now free */
        }
        /* --- Construct new element ---------------------------------------- */
        dk = 0 ;
        nv [k] = -nvk ;                     /* flag k as in Lk */
        p = Cp [k] ;
        pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
        pk2 = pk1 ;
        for (k1 = 1 ; k1 <= elenk + 1 ; k1++)
        {
            if (k1 > elenk)
            {
                e = k ;                     /* search the nodes in k */
                pj = p ;                    /* list of nodes starts at Ci[pj]*/
                ln = len [k] - elenk ;      /* length of list of nodes in k */
            }
            else
            {
                e = Ci [p++] ;              /* search the nodes in e */
                pj = Cp [e] ;
                ln = len [e] ;              /* length of list of nodes in e */
            }
            for (k2 = 1 ; k2 <= ln ; k2++)
            {
                i = Ci [pj++] ;
                if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
                dk += nvi ;                 /* degree[Lk] += size of node i */
                nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
                Ci [pk2++] = i ;            /* place i in Lk */
                if (next [i] != -1) last [next [i]] = last [i] ;
                if (last [i] != -1)         /* remove i from degree list */
                {
                    next [last [i]] = next [i] ;
                }
                else
                {
                    head [degree [i]] = next [i] ;
                }
            }
            if (e != k)
            {
                Cp [e] = CS_FLIP (k) ;      /* absorb e into k */
                w [e] = 0 ;                 /* e is now a dead element */
            }
        }
        if (elenk != 0) cnz = pk2 ;         /* Ci [cnz...nzmax] is free */
        degree [k] = dk ;                   /* external degree of k - |Lk\i| */
        Cp [k] = pk1 ;                      /* element k is in Ci[pk1..pk2-1] */
        len [k] = pk2 - pk1 ;
        elen [k] = -2 ;                     /* k is now an element */
        /* --- Find set differences ----------------------------------------- */
        mark = cs_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan 1: find |Le\Lk| */
        {
            i = Ci [pk] ;
            if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
            nvi = -nv [i] ;                      /* nv [i] was negated */
            wnvi = mark - nvi ;
            for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++)  /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] >= mark)
                {
                    w [e] -= nvi ;          /* decrement |Le\Lk| */
                }
                else if (w [e] != 0)        /* ensure e is a live element */
                {
                    w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
                }
            }
        }
        /* --- Degree update ------------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan2: degree update */
        {
            i = Ci [pk] ;                   /* consider node i in Lk */
            p1 = Cp [i] ;
            p2 = p1 + elen [i] - 1 ;
            pn = p1 ;
            for (h = 0, d = 0, p = p1 ; p <= p2 ; p++)    /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] != 0)             /* e is an unabsorbed element */
                {
                    dext = w [e] - mark ;   /* dext = |Le\Lk| */
                    if (dext > 0)
                    {
                        d += dext ;         /* sum up the set differences */
                        Ci [pn++] = e ;     /* keep e in Ei */
                        h += e ;            /* compute the hash of node i */
                    }
                    else
                    {
                        Cp [e] = CS_FLIP (k) ;  /* aggressive absorb. e->k */
                        w [e] = 0 ;             /* e is a dead element */
                    }
                }
            }
            elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
            p3 = pn ;
            p4 = p1 + len [i] ;
            for (p = p2 + 1 ; p < p4 ; p++) /* prune edges in Ai */
            {
                j = Ci [p] ;
                if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
                d += nvj ;                  /* degree(i) += |j| */
                Ci [pn++] = j ;             /* place j in node list of i */
                h += j ;                    /* compute hash for node i */
            }
            if (d == 0)                     /* check for mass elimination */
            {
                Cp [i] = CS_FLIP (k) ;      /* absorb i into k */
                nvi = -nv [i] ;
                dk -= nvi ;                 /* |Lk| -= |i| */
                nvk += nvi ;                /* |k| += nv[i] */
                nel += nvi ;
                nv [i] = 0 ;
                elen [i] = -1 ;             /* node i is dead */
            }
            else
            {
                degree [i] = CS_MIN (degree [i], d) ;   /* update degree(i) */
                Ci [pn] = Ci [p3] ;         /* move first node to end */
                Ci [p3] = Ci [p1] ;         /* move 1st el. to end of Ei */
                Ci [p1] = k ;               /* add k as 1st element in of Ei */
                len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
                h = ((h<0) ? (-h):h) % n ;  /* finalize hash of i */
                next [i] = hhead [h] ;      /* place i in hash bucket */
                hhead [h] = i ;
                last [i] = h ;              /* save hash of i in last[i] */
            }
        }                                   /* scan2 is done */
        degree [k] = dk ;                   /* finalize |Lk| */
        lemax = CS_MAX (lemax, dk) ;
        mark = cs_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
        /* --- Supernode detection ------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)
        {
            i = Ci [pk] ;
            if (nv [i] >= 0) continue ;         /* skip if i is dead */
            h = last [i] ;                      /* scan hash bucket of node i */
            i = hhead [h] ;
            hhead [h] = -1 ;                    /* hash bucket will be empty */
            for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++)
            {
                ln = len [i] ;
                eln = elen [i] ;
                for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++) w [Ci [p]] = mark;
                jlast = i ;
                for (j = next [i] ; j != -1 ; ) /* compare i with all j */
                {
                    ok = (len [j] == ln) && (elen [j] == eln) ;
                    for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++)
                    {
                        if (w [Ci [p]] != mark) ok = 0 ;    /* compare i and j*/
                    }
                    if (ok)                     /* i and j are identical */
                    {
                        Cp [j] = CS_FLIP (i) ;  /* absorb j into i */
                        nv [i] += nv [j] ;
                        nv [j] = 0 ;
                        elen [j] = -1 ;         /* node j is dead */
                        j = next [j] ;          /* delete j from hash bucket */
                        next [jlast] = j ;
                    }
                    else
                    {
                        jlast = j ;             /* j and i are different */
                        j = next [j] ;
                    }
                }
            }
        }
        /* --- Finalize new element------------------------------------------ */
        for (p = pk1, pk = pk1 ; pk < pk2 ; pk++)   /* finalize Lk */
        {
            i = Ci [pk] ;
            if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
            nv [i] = nvi ;                      /* restore nv[i] */
            d = degree [i] + dk - nvi ;         /* compute external degree(i) */
            d = CS_MIN (d, n - nel - nvi) ;
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;               /* put i back in degree list */
            last [i] = -1 ;
            head [d] = i ;
            mindeg = CS_MIN (mindeg, d) ;       /* find new minimum degree */
            degree [i] = d ;
            Ci [p++] = i ;                      /* place i in Lk */
        }
        nv [k] = nvk ;                      /* # nodes absorbed into k */
        if ((len [k] = p-pk1) == 0)         /* length of adj list of element k*/
        {
            Cp [k] = -1 ;                   /* k is a root of the tree */
            w [k] = 0 ;                     /* k is now a dead element */
        }
        if (elenk != 0) cnz = p ;           /* free unused space in Lk */
    }
    /* --- Postordering ----------------------------------------------------- */
    for (i = 0 ; i < n ; i++) Cp [i] = CS_FLIP (Cp [i]) ;/* fix assembly tree */
    for (j = 0 ; j <= n ; j++) head [j] = -1 ;
    for (j = n ; j >= 0 ; j--)              /* place unordered nodes in lists */
    {
        if (nv [j] > 0) continue ;          /* skip if j is an element */
        next [j] = head [Cp [j]] ;          /* place j in list of its parent */
        head [Cp [j]] = j ;
    }
    for (e = n ; e >= 0 ; e--)              /* place elements in lists */
    {
        if (nv [e] <= 0) continue ;         /* skip unless e is an element */
        if (Cp [e] != -1)
        {
            next [e] = head [Cp [e]] ;      /* place e in list of its parent */
            head [Cp [e]] = e ;
        }
    }
    for (k = 0, i = 0 ; i <= n ; i++)       /* postorder the assembly tree */
    {
        if (Cp [i] == -1) k = CSTdfs (i, k, head, next, P, w) ;
    }
    return (CSIdone (P, C, W, 1)) ;
}

/* L = chol (A, [pinv parent cp]), pinv is optional */
CSN *CSChol (const CS *A, const css *S)
{
    Real d, lki, *Lx, *x, *Cx ;
    I64 top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
    CS *L, *C, *E ;
    CSN *N ;
    if (!CS_CSC (A) || !S || !S->cp || !S->parent) return (nullptr) ;
    n = A->n ;
    N = (CSN *)CSCalloc (1, sizeof (CSN)) ;       /* allocate result */
    c = (I64 *)CSMalloc (2*n, sizeof (I64)) ;     /* get I64 workspace */
    x = (Real *)CSMalloc (n, sizeof (Real)) ;    /* get Real workspace */
    cp = S->cp ; pinv = S->pinv ; parent = S->parent ;
    C = pinv ? CSSymperm (A, pinv, 1) : ((CS *) A) ;
    E = pinv ? C : nullptr ;           /* E is alias for A, or a copy E=A(p,p) */
    if (!N || !c || !x || !C) return (CSNdone (N, E, c, x, 0)) ;
    s = c + n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    N->L = L = CSSpalloc (n, n, cp [n], 1, 0) ;    /* allocate result */
    if (!L) return (CSNdone (N, E, c, x, 0)) ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (k = 0 ; k < n ; k++) Lp [k] = c [k] = cp [k] ;
    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        top = CSEreach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
        x [k] = 0 ;                                 /* x (0:k) is now zero */
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
        }
        d = x [k] ;                     /* d = C(k,k) */
        x [k] = 0 ;                     /* clear x for k+1st iteration */
        /* --- Triangular solve --------------------------------------------- */
        for ( ; top < n ; top++)    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        {
            i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            lki = x [i] / Lx [Lp [i]] ; /* L(k,i) = x (i) / L(i,i) */
            x [i] = 0 ;                 /* clear x for k+1st iteration */
            for (p = Lp [i] + 1 ; p < c [i] ; p++)
            {
                x [Li [p]] -= Lx [p] * lki ;
            }
            d -= lki * lki ;            /* d = d - L(k,i)*L(k,i) */
            p = c [i]++ ;
            Li [p] = k ;                /* store L(k,i) in column i */
            Lx [p] = lki ;
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        if (d <= 0) return (CSNdone (N, E, c, x, 0)) ; /* not pos def */
        p = c [k]++ ;
        Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
        Lx [p] = sqrt (d) ;
    }
    Lp [n] = cp [n] ;               /* finalize L */
    return (CSNdone (N, E, c, x, 1)) ; /* success: free E,s,x; return N */
}

/* breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3) */
static I64 cs_bfs (const CS *A, I64 n, I64 *wi, I64 *wj, I64 *queue,
    const I64 *imatch, const I64 *jmatch, I64 mark)
{
    I64 *Ap, *Ai, head = 0, tail = 0, j, i, p, j2 ;
    CS *C ;
    for (j = 0 ; j < n ; j++)           /* place all unmatched nodes in queue */
    {
        if (imatch [j] >= 0) continue ; /* skip j if matched */
        wj [j] = 0 ;                    /* j in set C0 (R0 if transpose) */
        queue [tail++] = j ;            /* place unmatched col j in queue */
    }
    if (tail == 0) return (1) ;         /* quick return if no unmatched nodes */
    C = (mark == 1) ? ((CS *) A) : CSTranspose (A, 0) ;
    if (!C) return (0) ;                /* bfs of C=A' to find R3,C3 from R0 */
    Ap = C->p ; Ai = C->i ;
    while (head < tail)                 /* while queue is not empty */
    {
        j = queue [head++] ;            /* get the head of the queue */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (wi [i] >= 0) continue ; /* skip if i is marked */
            wi [i] = mark ;             /* i in set R1 (C3 if transpose) */
            j2 = jmatch [i] ;           /* traverse alternating path to j2 */
            if (wj [j2] >= 0) continue ;/* skip j2 if it is marked */
            wj [j2] = mark ;            /* j2 in set C1 (R3 if transpose) */
            queue [tail++] = j2 ;       /* add j2 to queue */
        }
    }
    if (mark != 1) CSSpfree (C) ;      /* free A' if it was created */
    return (1) ;
}

/* collect matched rows and columns into p and q */
static void cs_matched (I64 n, const I64 *wj, const I64 *imatch, I64 *p, I64 *q,
    I64 *cc, I64 *rr, I64 set, I64 mark)
{
    I64 kc = cc [set], j ;
    I64 kr = rr [set-1] ;
    for (j = 0 ; j < n ; j++)
    {
        if (wj [j] != mark) continue ;      /* skip if j is not in C set */
        p [kr++] = imatch [j] ;
        q [kc++] = j ;
    }
    cc [set+1] = kc ;
    rr [set] = kr ;
}

/* collect unmatched rows into the permutation vector p */
static void cs_unmatched (I64 m, const I64 *wi, I64 *p, I64 *rr, I64 set)
{
    I64 i, kr = rr [set] ;
    for (i = 0 ; i < m ; i++) if (wi [i] == 0) p [kr++] = i ;
    rr [set+1] = kr ;
}

/* return 1 if row i is in R2 */
static I64 cs_rprune (I64 i, I64, Real, void *other)
{
    I64 *rr = (I64 *) other ;
    return (i >= rr [1] && i < rr [2]) ;
}

/* Given A, compute coarse and then fine dmperm */
CSD *CSDmperm (const CS *A, I64 seed)
{
    I64 m, n, i, j, k, cnz, nc, *jmatch, *imatch, *wi, *wj, *pinv, *Cp, *Ci,
        *ps, *rs, nb1, nb2, *p, *q, *cc, *rr, *r, *s, ok ;
    CS *C ;
    CSD *D, *scc ;
    /* --- Maximum matching ------------------------------------------------- */
    if (!CS_CSC (A)) return (nullptr) ;            /* check inputs */
    m = A->m ; n = A->n ;
    D = CSDalloc (m, n) ;                      /* allocate result */
    if (!D) return (nullptr) ;
    p = D->p ; q = D->q ; r = D->r ; s = D->s ; cc = D->cc ; rr = D->rr ;
    jmatch = CSMaxtrans (A, seed) ;            /* max transversal */
    imatch = jmatch + m ;                       /* imatch = inverse of jmatch */
    if (!jmatch) return (CSDdone (D, nullptr, jmatch, 0)) ;
    /* --- Coarse decomposition --------------------------------------------- */
    wi = r ; wj = s ;                           /* use r and s as workspace */
    for (j = 0 ; j < n ; j++) wj [j] = -1 ;     /* unmark all cols for bfs */
    for (i = 0 ; i < m ; i++) wi [i] = -1 ;     /* unmark all rows for bfs */
    cs_bfs (A, n, wi, wj, q, imatch, jmatch, 1) ;       /* find C1, R1 from C0*/
    ok = cs_bfs (A, m, wj, wi, p, jmatch, imatch, 3) ;  /* find R3, C3 from R0*/
    if (!ok) return (CSDdone (D, nullptr, jmatch, 0)) ;
    cs_unmatched (n, wj, q, cc, 0) ;                    /* unmatched set C0 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 1, 1) ;    /* set R1 and C1 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 2, -1) ;   /* set R2 and C2 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 3, 3) ;    /* set R3 and C3 */
    cs_unmatched (m, wi, p, rr, 3) ;                    /* unmatched set R0 */
    CSFree (jmatch) ;
    /* --- Fine decomposition ----------------------------------------------- */
    pinv = CSPinv (p, m) ;         /* pinv=p' */
    if (!pinv) return (CSDdone (D, nullptr, nullptr, 0)) ;
    C = CSPermute (A, pinv, q, 0) ;/* C=A(p,q) (it will hold A(R2,C2)) */
    CSFree (pinv) ;
    if (!C) return (CSDdone (D, nullptr, nullptr, 0)) ;
    Cp = C->p ;
    nc = cc [3] - cc [2] ;          /* delete cols C0, C1, and C3 from C */
    if (cc [2] > 0) for (j = cc [2] ; j <= cc [3] ; j++) Cp [j-cc[2]] = Cp [j] ;
    C->n = nc ;
    if (rr [2] - rr [1] < m)        /* delete rows R0, R1, and R3 from C */
    {
        CSFkeep (C, cs_rprune, rr) ;
        cnz = Cp [nc] ;
        Ci = C->i ;
        if (rr [1] > 0) for (k = 0 ; k < cnz ; k++) Ci [k] -= rr [1] ;
    }
    C->m = nc ;
    scc = CSScc (C) ;              /* find strongly connected components of C*/
    if (!scc) return (CSDdone (D, C, nullptr, 0)) ;
    /* --- Combine coarse and fine decompositions --------------------------- */
    ps = scc->p ;                   /* C(ps,ps) is the permuted matrix */
    rs = scc->r ;                   /* kth block is rs[k]..rs[k+1]-1 */
    nb1 = scc->nb  ;                /* # of blocks of A(R2,C2) */
    for (k = 0 ; k < nc ; k++) wj [k] = q [ps [k] + cc [2]] ;
    for (k = 0 ; k < nc ; k++) q [k + cc [2]] = wj [k] ;
    for (k = 0 ; k < nc ; k++) wi [k] = p [ps [k] + rr [1]] ;
    for (k = 0 ; k < nc ; k++) p [k + rr [1]] = wi [k] ;
    nb2 = 0 ;                       /* create the fine block partitions */
    r [0] = s [0] = 0 ;
    if (cc [2] > 0) nb2++ ;         /* leading coarse block A (R1, [C0 C1]) */
    for (k = 0 ; k < nb1 ; k++)     /* coarse block A (R2,C2) */
    {
        r [nb2] = rs [k] + rr [1] ; /* A (R2,C2) splits into nb1 fine blocks */
        s [nb2] = rs [k] + cc [2] ;
        nb2++ ;
    }
    if (rr [2] < m)
    {
        r [nb2] = rr [2] ;          /* trailing coarse block A ([R3 R0], C3) */
        s [nb2] = cc [3] ;
        nb2++ ;
    }
    r [nb2] = m ;
    s [nb2] = n ;
    D->nb = nb2 ;
    CSDfree (scc) ;
    return (CSDdone (D, C, nullptr, 1)) ;
}

static I64 cs_tol (I64, I64, Real aij, void *tol)
{
    return (fabs (aij) > *((Real *) tol)) ;
}
I64 CSDroptol (CS *A, Real tol)
{
    return (CSFkeep (A, &cs_tol, &tol)) ;    /* keep all large entries */
}

static I64 cs_nonzero (I64, I64, Real aij, void *)
{
    return (aij != 0) ;
}
I64 CSDropzeros (CS *A)
{
    return (CSFkeep (A, &cs_nonzero, nullptr)) ;  /* keep all nonzero entries */
} 

/* apply the ith Householder vector to x */
I64 CSHapply (const CS *V, I64 i, Real beta, Real *x)
{
    I64 p, *Vp, *Vi ;
    Real *Vx, tau = 0 ;
    if (!CS_CSC (V) || !x) return (0) ;     /* check inputs */
    Vp = V->p ; Vi = V->i ; Vx = V->x ;
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* tau = v'*x */
    {
        tau += Vx [p] * x [Vi [p]] ;
    }
    tau *= beta ;                           /* tau = beta*(v'*x) */
    for (p = Vp [i] ; p < Vp [i+1] ; p++)   /* x = x - v*tau */
    {
        x [Vi [p]] -= Vx [p] * tau ;
    }
    return (1) ;
}

/* x(p) = b, for dense vectors x and b; p=nullptr denotes identity */
I64 CSIpvec (const I64 *p, const Real *b, Real *x, I64 n)
{
    I64 k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [p ? p [k] : k] = b [k] ;
    return (1) ;
}

/* solve Lx=b where x and b are dense.  x=b on input, solution on output. */
I64 CSLsolve (const CS *L, Real *x)
{
    I64 p, j, n, *Lp, *Li ;
    Real *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] /= Lx [Lp [j]] ;
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}

/* solve L'x=b where x and b are dense.  x=b on input, solution on output. */
I64 CSLtsolve (const CS *L, Real *x)
{
    I64 p, j, n, *Lp, *Li ;
    Real *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= Lx [p] * x [Li [p]] ;
        }
        x [j] /= Lx [Lp [j]] ;
    }
    return (1) ;
}

/* [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess */
CSN *CSLu (const CS *A, const css *S, Real tol)
{
    CS *L, *U ;
    CSN *N ;
    Real pivot, *Lx, *Ux, *x,  a, t ;
    I64 *Lp, *Li, *Up, *Ui, *pinv, *xi, *q, n, ipiv, k, top, p, i, col, lnz,unz;
    if (!CS_CSC (A) || !S) return (nullptr) ;          /* check inputs */
    n = A->n ;
    q = S->q ; lnz = S->lnz ; unz = S->unz ;
    x = (Real *)CSMalloc (n, sizeof (Real)) ;            /* get Real workspace */
    xi = (I64 *)CSMalloc (2*n, sizeof (I64)) ;            /* get I64 workspace */
    N = (CSN *)CSCalloc (1, sizeof (CSN)) ;               /* allocate result */
    if (!x || !xi || !N) return (CSNdone (N, nullptr, xi, x, 0)) ;
    N->L = L = CSSpalloc (n, n, lnz, 1, 0) ;       /* allocate result L */
    N->U = U = CSSpalloc (n, n, unz, 1, 0) ;       /* allocate result U */
    N->pinv = pinv = (I64 *)CSMalloc (n, sizeof (I64)) ;  /* allocate result pinv */
    if (!L || !U || !pinv) return (CSNdone (N, nullptr, xi, x, 0)) ;
    Lp = L->p ; Up = U->p ;
    for (i = 0 ; i < n ; i++) x [i] = 0 ;           /* clear workspace */
    for (i = 0 ; i < n ; i++) pinv [i] = -1 ;       /* no rows pivotal yet */
    for (k = 0 ; k <= n ; k++) Lp [k] = 0 ;         /* no cols of L yet */
    lnz = unz = 0 ;
    for (k = 0 ; k < n ; k++)       /* compute L(:,k) and U(:,k) */
    {
        /* --- Triangular solve --------------------------------------------- */
        Lp [k] = lnz ;              /* L(:,k) starts here */
        Up [k] = unz ;              /* U(:,k) starts here */
        if ((lnz + n > L->nzmax && !CSSprealloc (L, 2*L->nzmax + n)) ||
            (unz + n > U->nzmax && !CSSprealloc (U, 2*U->nzmax + n)))
        {
            return (CSNdone (N, nullptr, xi, x, 0)) ;
        }
        Li = L->i ; Lx = L->x ; Ui = U->i ; Ux = U->x ;
        col = q ? (q [k]) : k ;
        top = CSSpsolve (L, A, col, xi, x, pinv, 1) ;  /* x = L\A(:,col) */
        /* --- Find pivot --------------------------------------------------- */
        ipiv = -1 ;
        a = -1 ;
        for (p = top ; p < n ; p++)
        {
            i = xi [p] ;            /* x(i) is nonzero */
            if (pinv [i] < 0)       /* row i is not yet pivotal */
            {
                if ((t = fabs (x [i])) > a)
                {
                    a = t ;         /* largest pivot candidate so far */
                    ipiv = i ;
                }
            }
            else                    /* x(i) is the entry U(pinv[i],k) */
            {
                Ui [unz] = pinv [i] ;
                Ux [unz++] = x [i] ;
            }
        }
        if (ipiv == -1 || a <= 0) return (CSNdone (N, nullptr, xi, x, 0)) ;
        if (pinv [col] < 0 && fabs (x [col]) >= a*tol) ipiv = col ;
        /* --- Divide by pivot ---------------------------------------------- */
        pivot = x [ipiv] ;          /* the chosen pivot */
        Ui [unz] = k ;              /* last entry in U(:,k) is U(k,k) */
        Ux [unz++] = pivot ;
        pinv [ipiv] = k ;           /* ipiv is the kth pivot row */
        Li [lnz] = ipiv ;           /* first entry in L(:,k) is L(k,k) = 1 */
        Lx [lnz++] = 1 ;
        for (p = top ; p < n ; p++) /* L(k+1:n,k) = x / pivot */
        {
            i = xi [p] ;
            if (pinv [i] < 0)       /* x(i) is an entry in L(:,k) */
            {
                Li [lnz] = i ;      /* save unpermuted row in L */
                Lx [lnz++] = x [i] / pivot ;    /* scale pivot column */
            }
            x [i] = 0 ;             /* x [0..n-1] = 0 for next k */
        }
    }
    /* --- Finalize L and U ------------------------------------------------- */
    Lp [n] = lnz ;
    Up [n] = unz ;
    Li = L->i ;                     /* fix row indices of L for final pinv */
    for (p = 0 ; p < lnz ; p++) Li [p] = pinv [Li [p]] ;
    CSSprealloc (L, 0) ;           /* remove extra space from L and U */
    CSSprealloc (U, 0) ;
    return (CSNdone (N, nullptr, xi, x, 1)) ;     /* success */
}

/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1. */
CS *CSPermute (const CS *A, const I64 *pinv, const I64 *q, I64 values)
{
    I64 t, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
    Real *Cx, *Ax ;
    CS *C ;
    if (!CS_CSC (A)) return (nullptr) ;    /* check inputs */
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = CSSpalloc (m, n, Ap [n], values && Ax != nullptr, 0) ;  /* alloc result */
    if (!C) return (CSDone (C, nullptr, nullptr, 0)) ;   /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k = 0 ; k < n ; k++)
    {
        Cp [k] = nz ;                   /* column k of C is column q[k] of A */
        j = q ? (q [k]) : k ;
        for (t = Ap [j] ; t < Ap [j+1] ; t++)
        {
            if (Cx) Cx [nz] = Ax [t] ;  /* row i of A is row pinv[i] of C */
            Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
        }
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    return (CSDone (C, nullptr, nullptr, 1)) ;
}

/* pinv = p', or p = pinv' */
I64 *CSPinv (I64 const *p, I64 n)
{
    I64 k, *pinv ;
    if (!p) return (nullptr) ;                     /* p = nullptr denotes identity */
    pinv = (I64 *)CSMalloc (n, sizeof (I64)) ;        /* allocate result */
    if (!pinv) return (nullptr) ;                  /* out of memory */
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;/* invert the permutation */
    return (pinv) ;                             /* return result */
}

/* x = b(p), for dense vectors x and b; p=nullptr denotes identity */
I64 CSPvec (const I64 *p, const Real *b, Real *x, I64 n)
{
    I64 k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [k] = b [p ? p [k] : k] ;
    return (1) ;
}

/* sparse QR factorization [V,beta,pinv,R] = qr (A) */
CSN *CSQr (const CS *A, const css *S)
{
    Real *Rx, *Vx, *Ax, *x,  *Beta ;
    I64 i, k, p, /*m,*/ n, vnz, p1, top, m2, len, col, rnz, *s, *leftmost, *Ap, *Ai,
        *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q ;
    CS *R, *V ;
    CSN *N ;
    if (!CS_CSC (A) || !S) return (nullptr) ;
    /*m = A->m ;*/ n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    q = S->q ; parent = S->parent ; pinv = S->pinv ; m2 = S->m2 ;
    vnz = S->lnz ; rnz = S->unz ; leftmost = S->leftmost ;
    w = (I64 *)CSMalloc (m2+n, sizeof (I64)) ;            /* get I64 workspace */
    x = (Real *)CSMalloc (m2, sizeof (Real)) ;           /* get Real workspace */
    N = (CSN *)CSCalloc (1, sizeof (CSN)) ;               /* allocate result */
    if (!w || !x || !N) return (CSNdone (N, nullptr, w, x, 0)) ;
    s = w + m2 ;                                    /* s is size n */
    for (k = 0 ; k < m2 ; k++) x [k] = 0 ;          /* clear workspace x */
    N->L = V = CSSpalloc (m2, n, vnz, 1, 0) ;      /* allocate result V */
    N->U = R = CSSpalloc (m2, n, rnz, 1, 0) ;      /* allocate result R */
    N->B = Beta = (Real *)CSMalloc (n, sizeof (Real)) ;  /* allocate result Beta */
    if (!R || !V || !Beta) return (CSNdone (N, nullptr, w, x, 0)) ;
    Rp = R->p ; Ri = R->i ; Rx = R->x ;
    Vp = V->p ; Vi = V->i ; Vx = V->x ;
    for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */
    rnz = 0 ; vnz = 0 ;
    for (k = 0 ; k < n ; k++)               /* compute V and R */
    {
        Rp [k] = rnz ;                      /* R(:,k) starts here */
        Vp [k] = p1 = vnz ;                 /* V(:,k) starts here */
        w [k] = k ;                         /* add V(k,k) to pattern of V */
        Vi [vnz++] = k ;
        top = n ;
        col = q ? q [k] : k ;
        for (p = Ap [col] ; p < Ap [col+1] ; p++)   /* find R(:,k) pattern */
        {
            i = leftmost [Ai [p]] ;         /* i = min(find(A(i,q))) */
            for (len = 0 ; w [i] != k ; i = parent [i]) /* traverse up to k */
            {
                s [len++] = i ;
                w [i] = k ;
            }
            while (len > 0) s [--top] = s [--len] ; /* push path on stack */
            i = pinv [Ai [p]] ;             /* i = permuted row of A(:,col) */
            x [i] = Ax [p] ;                /* x (i) = A(:,col) */
            if (i > k && w [i] < k)         /* pattern of V(:,k) = x (k+1:m) */
            {
                Vi [vnz++] = i ;            /* add i to pattern of V(:,k) */
                w [i] = k ;
            }
        }
        for (p = top ; p < n ; p++) /* for each i in pattern of R(:,k) */
        {
            i = s [p] ;                     /* R(i,k) is nonzero */
            CSHapply (V, i, Beta [i], x) ; /* apply (V(i),Beta(i)) to x */
            Ri [rnz] = i ;                  /* R(i,k) = x(i) */
            Rx [rnz++] = x [i] ;
            x [i] = 0 ;
            if (parent [i] == k) vnz = CSScatter (V, i, 0, w, nullptr, k, V, vnz);
        }
        for (p = p1 ; p < vnz ; p++)        /* gather V(:,k) = x */
        {
            Vx [p] = x [Vi [p]] ;
            x [Vi [p]] = 0 ;
        }
        Ri [rnz] = k ;                     /* R(k,k) = norm (x) */
        Rx [rnz++] = CSHouse (Vx+p1, Beta+k, vnz-p1) ; /* [v,beta]=house(x) */
    }
    Rp [n] = rnz ;                          /* finalize R */
    Vp [n] = vnz ;                          /* finalize V */
    return (CSNdone (N, nullptr, w, x, 1)) ;  /* success */
}

/* ordering and symbolic analysis for a Cholesky factorization */
css *CSSchol (I64 order, const CS *A)
{
    I64 n, *c, *post, *P ;
    CS *C ;
    css *S ;
    if (!CS_CSC (A)) return (nullptr) ;        /* check inputs */
    n = A->n ;
    S = (css *)CSCalloc (1, sizeof (css)) ;       /* allocate result S */
    if (!S) return (nullptr) ;                 /* out of memory */
    P = CSAmd (order, A) ;                 /* P = amd(A+A'), or natural */
    S->pinv = CSPinv (P, n) ;              /* find inverse permutation */
    CSFree (P) ;
    if (order && !S->pinv) return (CSSfree (S)) ;
    C = CSSymperm (A, S->pinv, 0) ;        /* C = spones(triu(A(P,P))) */
    S->parent = CSEtree (C, 0) ;           /* find etree of C */
    post = CSPost (S->parent, n) ;         /* postorder the etree */
    c = CSCounts (C, S->parent, post, 0) ; /* find column counts of chol(C) */
    CSFree (post) ;
    CSSpfree (C) ;
    S->cp = (I64 *)CSMalloc (n+1, sizeof (I64)) ; /* allocate result S->cp */
    S->unz = S->lnz = CSCumsum (S->cp, c, n) ; /* find column pointers for L */
    CSFree (c) ;
    return ((S->lnz >= 0) ? S : CSSfree (S)) ;
}

/* compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent */
static I64 cs_vcount (const CS *A, css *S)
{
    I64 i, k, p, pa, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i, *next, *head,
        *tail, *nque, *pinv, *leftmost, *w, *parent = S->parent ;
    S->pinv = pinv = (I64 *)CSMalloc (m+n, sizeof (I64)) ;        /* allocate pinv, */
    S->leftmost = leftmost = (I64 *)CSMalloc (m, sizeof (I64)) ;  /* and leftmost */
    w = (I64 *)CSMalloc (m+3*n, sizeof (I64)) ;   /* get workspace */
    if (!pinv || !w || !leftmost)
    {
        CSFree (w) ;                       /* pinv and leftmost freed later */
        return (0) ;                        /* out of memory */
    }
    next = w ; head = w + m ; tail = w + m + n ; nque = w + m + 2*n ;
    for (k = 0 ; k < n ; k++) head [k] = -1 ;   /* queue k is empty */
    for (k = 0 ; k < n ; k++) tail [k] = -1 ;
    for (k = 0 ; k < n ; k++) nque [k] = 0 ;
    for (i = 0 ; i < m ; i++) leftmost [i] = -1 ;
    for (k = n-1 ; k >= 0 ; k--)
    {
        for (p = Ap [k] ; p < Ap [k+1] ; p++)
        {
            leftmost [Ai [p]] = k ;         /* leftmost[i] = min(find(A(i,:)))*/
        }
    }
    for (i = m-1 ; i >= 0 ; i--)            /* scan rows in reverse order */
    {
        pinv [i] = -1 ;                     /* row i is not yet ordered */
        k = leftmost [i] ;
        if (k == -1) continue ;             /* row i is empty */
        if (nque [k]++ == 0) tail [k] = i ; /* first row in queue k */
        next [i] = head [k] ;               /* put i at head of queue k */
        head [k] = i ;
    }
    S->lnz = 0 ;
    S->m2 = m ;
    for (k = 0 ; k < n ; k++)               /* find row permutation and nnz(V)*/
    {
        i = head [k] ;                      /* remove row i from queue k */
        S->lnz++ ;                          /* count V(k,k) as nonzero */
        if (i < 0) i = S->m2++ ;            /* add a fictitious row */
        pinv [i] = k ;                      /* associate row i with V(:,k) */
        if (--nque [k] <= 0) continue ;     /* skip if V(k+1:m,k) is empty */
        S->lnz += nque [k] ;                /* nque [k] is nnz (V(k+1:m,k)) */
        if ((pa = parent [k]) != -1)        /* move all rows to parent of k */
        {
            if (nque [pa] == 0) tail [pa] = tail [k] ;
            next [tail [k]] = head [pa] ;
            head [pa] = next [i] ;
            nque [pa] += nque [k] ;
        }
    }
    for (i = 0 ; i < m ; i++) if (pinv [i] < 0) pinv [i] = k++ ;
    CSFree (w) ;
    return (1) ;
}

/* symbolic ordering and analysis for QR or LU */
css *CSSqr (I64 order, const CS *A, I64 qr)
{
    I64 n, k, ok = 1, *post ;
    css *S ;
    if (!CS_CSC (A)) return (nullptr) ;        /* check inputs */
    n = A->n ;
    S = (css *)CSCalloc (1, sizeof (css)) ;       /* allocate result S */
    if (!S) return (nullptr) ;                 /* out of memory */
    S->q = CSAmd (order, A) ;              /* fill-reducing ordering */
    if (order && !S->q) return (CSSfree (S)) ;
    if (qr)                                 /* QR symbolic analysis */
    {
        CS *C = order ? CSPermute (A, nullptr, S->q, 0) : ((CS *) A) ;
        S->parent = CSEtree (C, 1) ;       /* etree of C'*C, where C=A(:,q) */
        post = CSPost (S->parent, n) ;
        S->cp = CSCounts (C, S->parent, post, 1) ;  /* col counts chol(C'*C) */
        CSFree (post) ;
        ok = C && S->parent && S->cp && cs_vcount (C, S) ;
        if (ok) for (S->unz = 0, k = 0 ; k < n ; k++) S->unz += S->cp [k] ;
        if (order) CSSpfree (C) ;
    }
    else
    {
        S->unz = 4*(A->p [n]) + n ;         /* for LU factorization only, */
        S->lnz = S->unz ;                   /* guess nnz(L) and nnz(U) */
    }
    return (ok ? S : CSSfree (S)) ;        /* return result S */
}

/* C = A(p,p) where A and C are symmetric the upper part stored; pinv not p */
CS *CSSymperm (const CS *A, const I64 *pinv, I64 values)
{
    I64 i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w ;
    Real *Cx, *Ax ;
    CS *C ;
    if (!CS_CSC (A)) return (nullptr) ;                    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = CSSpalloc (n, n, Ap [n], values && (Ax != nullptr), 0) ; /* alloc result*/
    w = (I64 *)CSCalloc (n, sizeof (I64)) ;                   /* get workspace */
    if (!C || !w) return (CSDone (C, w, nullptr, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (j = 0 ; j < n ; j++)           /* count entries in each column of C */
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (i > j) continue ;       /* skip lower triangular part of A */
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            w [CS_MAX (i2, j2)]++ ;     /* column count of C */
        }
    }
    CSCumsum (Cp, w, n) ;              /* compute column pointers of C */
    for (j = 0 ; j < n ; j++)
    {
        j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (i > j) continue ;       /* skip lower triangular part of A*/
            i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
            Ci [q = w [CS_MAX (i2, j2)]++] = CS_MIN (i2, j2) ;
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (CSDone (C, w, nullptr, 1)) ;  /* success; free workspace, return C */
}

/* sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1) */
I64 CSUpdown (CS *L, I64 sigma, const CS *C, const I64 *parent)
{
    I64 n, p, f, j, *Lp, *Li, *Cp, *Ci ;
    Real *Lx, *Cx, alpha, beta = 1, delta, gamma, w1, w2, *w, beta2 = 1 ;
    if (!CS_CSC (L) || !CS_CSC (C) || !parent) return (0) ;  /* check inputs */
    Lp = L->p ; Li = L->i ; Lx = L->x ; n = L->n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    if ((p = Cp [0]) >= Cp [1]) return (1) ;        /* return if C empty */
    w = (Real *)CSMalloc (n, sizeof (Real)) ;            /* get workspace */
    if (!w) return (0) ;                            /* out of memory */
    f = Ci [p] ;
    for ( ; p < Cp [1] ; p++) f = CS_MIN (f, Ci [p]) ;  /* f = min (find (C)) */
    for (j = f ; j != -1 ; j = parent [j]) w [j] = 0 ;  /* clear workspace w */
    for (p = Cp [0] ; p < Cp [1] ; p++) w [Ci [p]] = Cx [p] ; /* w = C */
    for (j = f ; j != -1 ; j = parent [j])          /* walk path f up to root */
    {
        p = Lp [j] ;
        alpha = w [j] / Lx [p] ;                    /* alpha = w(j) / L(j,j) */
        beta2 = beta*beta + sigma*alpha*alpha ;
        if (beta2 <= 0) break ;                     /* not positive definite */
        beta2 = sqrt (beta2) ;
        delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta) ;
        gamma = sigma * alpha / (beta2 * beta) ;
        Lx [p] = delta * Lx [p] + ((sigma > 0) ? (gamma * w [j]) : 0) ;
        beta = beta2 ;
        for (p++ ; p < Lp [j+1] ; p++)
        {
            w1 = w [Li [p]] ;
            w [Li [p]] = w2 = w1 - alpha * Lx [p] ;
            Lx [p] = delta * Lx [p] + gamma * ((sigma > 0) ? w1 : w2) ;
        }
    }
    CSFree (w) ;
    return (beta2 > 0) ;
}

/* solve Ux=b where x and b are dense.  x=b on input, solution on output. */
I64 CSUsolve (const CS *U, Real *x)
{
    I64 p, j, n, *Up, *Ui ;
    Real *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        x [j] /= Ux [Up [j+1]-1] ;
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [Ui [p]] -= Ux [p] * x [j] ;
        }
    }
    return (1) ;
}

/* solve U'x=b where x and b are dense.  x=b on input, solution on output. */
I64 CSUtsolve (const CS *U, Real *x)
{
    I64 p, j, n, *Up, *Ui ;
    Real *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [j] -= Ux [p] * x [Ui [p]] ;
        }
        x [j] /= Ux [Up [j+1]-1] ;
    }
    return (1) ;
}

/* free a numeric factorization */
CSN *CSNfree (CSN *N)
{
    if (!N) return (nullptr) ;     /* do nothing if N already nullptr */
    CSSpfree (N->L) ;
    CSSpfree (N->U) ;
    CSFree (N->pinv) ;
    CSFree (N->B) ;
    return ((CSN *) CSFree (N)) ;  /* free the CSN struct and return nullptr */
}

/* free a symbolic factorization */
css *CSSfree (css *S)
{
    if (!S) return (nullptr) ;     /* do nothing if S already nullptr */
    CSFree (S->pinv) ;
    CSFree (S->q) ;
    CSFree (S->parent) ;
    CSFree (S->cp) ;
    CSFree (S->leftmost) ;
    return ((css *) CSFree (S)) ;  /* free the css struct and return nullptr */
}

/* free a CSDmperm or CSScc result */
CSD *CSDfree (CSD *D)
{
    if (!D) return (nullptr) ;     /* do nothing if D already nullptr */
    CSFree (D->p) ;
    CSFree (D->q) ;
    CSFree (D->r) ;
    CSFree (D->s) ;
    return ((CSD *) CSFree (D)) ;  /* free the CSD struct and return nullptr */
}

}
