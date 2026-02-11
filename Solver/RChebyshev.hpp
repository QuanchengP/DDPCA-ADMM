#ifndef _RChebyshev_hpp
#define _RChebyshev_hpp

#include "../General/SparseMatrix.hpp"

namespace Ddpca{

//amgcl/relaxation/chebyshev.hpp
//relaxation: chebyshev
class RChebyshev{
public:
    // Centre of ellipse containing the eigenvalues of A:
    Real d;
    // Semi-major axis of ellipse containing the eigenvalues of A:
    Real c;
    I64 degree;
    Real higher;
    Real lower;
    I64 power_iters;

    RChebyshev() : degree(5), higher(1.0), lower(0.1666666666666666666), power_iters(15) {};
    ~RChebyshev() = default;

public:

    AlignedVectorRx residual;
    AlignedVectorRx p;

public:

    void Establish(
        const SparseMatrix& A, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

    void Apply(
        const SparseMatrix& A, 
        const AlignedVectorRx& rhs, 
        AlignedVectorRx& solution, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
}; // class RChebyshev

} // namespace Ddpca

#endif // _RChebyshev_hpp