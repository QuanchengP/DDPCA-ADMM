#ifndef _ConjugateGradient_hpp
#define _ConjugateGradient_hpp

#include "PCDiagonal.hpp"
#include "PCGeometricMultigrid.hpp"

namespace Ddpca {

static constexpr I64 PreconditionerNone = -1;
static constexpr I64 PreconditionerMultigrid = 0;
static constexpr I64 PreconditionerDiagonal = 1;

class ConjugateGradient{
    
public:

    I64 preconditioner = PreconditionerNone;
    PCGeometricMultigrid pcGM;
    PCDiagonal pcD;

    // x needs initial values
    void Solve(
        const SparseMatrix& tempStiff, const AlignedVectorRx& rhs, 
        AlignedVectorRx& x, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);
    
public:

    AlignedVectorRx residualError;//(numbX)
    AlignedVectorRx searchDirection;//numbX, 0.0
    AlignedVectorRx AxProduct;//numbX
    AlignedVectorRx preconditionedResidual;//numbX

    void Establish(
        const SparseMatrix& tempStiff, 
        const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const I64 nestLevel);

}; // class ConjugateGradient

} // namespace Ddpca

#endif // _ConjugateGradient_hpp