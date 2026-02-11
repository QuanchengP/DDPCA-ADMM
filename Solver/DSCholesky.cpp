#include "DSCholesky.hpp"

#include <cassert>
#include <execution>

namespace Ddpca{

void DSCholesky::Establish(const SparseMatrix& solverMatrix/*, const I64 numbThreads = 1*/){
    Log("    DSCholesky::Establish");
    //
    miniMatr = solverMatrix.ToCSparse();
    cholSymb = CSSchol (1, miniMatr) ;           /* ordering and symbolic analysis */
    assert(cholSymb != nullptr && "cholSymb is nullptr");
    cholNume = CSChol (miniMatr, cholSymb) ;     /* numeric Cholesky factorization */
    assert(cholNume != nullptr && "cholNume is nullptr");
    I64 n = solverMatrix.M;
    b = (double *)CSMalloc (n, sizeof (double));
    x = (double *)CSMalloc (n, sizeof (double)) ;   /* get workspace */
}

void DSCholesky::Solve(
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution/*, const I64 numbThreads = 1*/){
    //
    I64 n = rhs.size();
    std::copy(std::execution::unseq, rhs.begin(), rhs.end(), b);
    //
    CSIpvec (cholSymb->pinv, b, x, n) ;             /* x = P*b */
    CSLsolve (cholNume->L, x) ;                     /* x = L\x */
    CSLtsolve (cholNume->L, x) ;                    /* x = L'\x */
    CSPvec (cholSymb->pinv, x, b, n) ;              /* b = P'*x */
    //
    std::copy(std::execution::unseq, b, b + n, solution.begin());
}

} // namespace Ddpca