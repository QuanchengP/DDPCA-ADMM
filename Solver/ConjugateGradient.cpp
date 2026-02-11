#include "ConjugateGradient.hpp"
#include "../Mesh/Coordinate.hpp"

#include <cassert>
#include <cstring>

namespace Ddpca {

void ConjugateGradient::Solve(
    const SparseMatrix& tempStiff, 
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& x, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("    ConjugateGradient::Solve");
    I64 maxLevel = (pcGM.realProlong).size();
    const I64 numbX = rhs.size();
    Real tolerance = 1.0E-14 * NRM2(rhs);
    //
    I64 tempIteration = 0;
    MV(-1.0, tempStiff, x, 1.0, rhs, residualError, threadTask, nestLevel);
    switch(preconditioner){
        case PreconditionerNone:
            std::copy(std::execution::unseq, 
                residualError.begin(), residualError.end(), searchDirection.begin());
            break;
        case PreconditionerMultigrid:
            // will residualError be a better initial value for searchDirection
            // MV(inverseDiagonal, residualError, searchDirection, threadTask, nestLevel); //worse
            std::fill(std::execution::unseq, searchDirection.begin(), searchDirection.end(), 0.0);
            pcGM.Apply(maxLevel, residualError, searchDirection, threadTask, nestLevel);
            break;
        case PreconditionerDiagonal:
            pcD.Apply(residualError, searchDirection, threadTask, nestLevel);
            break;
    }
    Real deltaNew = DOT(residualError, searchDirection);
    //
    Real residualNorm, alpha, deltaOld, beta;
    while(tempIteration < numbX){
        //
        residualNorm = NRM2(residualError);
        if(residualNorm <= tolerance || tempIteration % 100 == 99){
            std::cout << "#Iteration: " << tempIteration
                << ", residual: " << residualNorm << "/" << tolerance << "\n";
            if(residualNorm <= tolerance){
                Log("    ConjugateGradient::Solve converged");
                break;
            }
        }
        //
        MV(tempStiff, searchDirection, AxProduct, threadTask, nestLevel);
        alpha = deltaNew / DOT(searchDirection, AxProduct);
        AXPY(alpha, searchDirection, x);
        AXPY(-alpha, AxProduct, residualError);
        switch(preconditioner){
            case PreconditionerNone:
                std::copy(std::execution::unseq, 
                    residualError.begin(), residualError.end(), preconditionedResidual.begin());
                break;
            case PreconditionerMultigrid:
                // memset(preconditionedResidual.data(), 0, numbX * sizeof(Real));// a better initial value?
                // MV(inverseDiagonal, residualError, preconditionedResidual, threadTask, nestLevel); // worse
                std::fill(std::execution::unseq, preconditionedResidual.begin(), preconditionedResidual.end(), 0.0);
                pcGM.Apply(maxLevel, residualError, preconditionedResidual, threadTask, nestLevel);
                break;
            case PreconditionerDiagonal:
                pcD.Apply(residualError, preconditionedResidual, threadTask, nestLevel);
                break;
        }
        deltaOld = deltaNew;
        deltaNew = DOT(residualError, preconditionedResidual);
        beta = deltaNew / deltaOld;
        AXPY(beta, searchDirection, preconditionedResidual, searchDirection);
        ++ tempIteration;
    }
}

void ConjugateGradient::Establish(
    const SparseMatrix& tempStiff, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    switch(preconditioner){
        case PreconditionerNone:
            break;
        case PreconditionerMultigrid:
            pcGM.Establish(threadTask, nestLevel);
            break;
        case PreconditionerDiagonal:
            pcD.Establish(tempStiff, threadTask, nestLevel);
            break;
    }
    const I64 numbX = tempStiff.M;
    //
    residualError.resize(numbX);
    searchDirection.resize(numbX, 0.0);
    AxProduct.resize(numbX);
    preconditionedResidual.resize(numbX);
}

} // namespace Ddpca