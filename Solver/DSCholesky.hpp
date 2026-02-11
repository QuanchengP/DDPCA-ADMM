#ifndef _DSCholesky_hpp
#define _DSCholesky_hpp
#include "CSparse/CS.hpp"

#include "../General/SparseMatrix.hpp"

namespace Ddpca{

class DSCholesky{

public:
    CS *miniMatr = nullptr;
    css *cholSymb = nullptr;
    CSN *cholNume = nullptr;
    Real *x = nullptr;
    Real *b = nullptr;

    ~DSCholesky(){
        if(!miniMatr) CSSpfree (miniMatr) ;
        if(!cholSymb) CSSfree (cholSymb) ;
        if(!cholNume) CSNfree (cholNume) ;
        if(!x) CSFree (x) ;
        if(!b) CSFree (b) ;
    }

public:

    void Establish(const SparseMatrix& solverMatrix/*, const I64 numbThreads = 1*/);

    void Solve(
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution/*, const I64 numbThreads = 1*/);

}; // class DSCholesky

} // namespace Ddpca

#endif // _DSCholesky_hpp