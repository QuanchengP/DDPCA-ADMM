#include "RChebyshev.hpp"
#include "../Mesh/Coordinate.hpp"

namespace Ddpca{

void RChebyshev::Establish(
    const SparseMatrix& A, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Real lambda_max = A.SpectralRadius(power_iters, threadTask, nestLevel);
    Real lambda_min = lower * lambda_max;
    lambda_max *= higher;
    d = (lambda_max + lambda_min) * 0.5;
    c = (lambda_max - lambda_min) * 0.5;
    //
    residual.resize(A.M, 0.0);
    p.resize(A.M, 0.0);
}

void RChebyshev::Apply(
    const SparseMatrix& A, 
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(A.M, numbPart, startIndex, endIndex);
    //
    alignas(nfsAlign) Real alpha = 0.0, beta = 0.0;
    for(I64 tk = 0; tk < degree; ++ tk){
        MV(-1.0, A, solution, 1.0, rhs, residual, threadTask, nestLevel);
        switch(tk){
            case 0:
                alpha = 1.0 / d;
                beta = 0.0;
                break;
            case 1:
                alpha = 2.0 * d / (2.0 * d * d - c * c);
                beta = alpha * d - 1.0;
                break;
            default:
                alpha = 1.0 / (d - 0.25 * alpha * c * c);
                beta = alpha * d - 1.0;
                break;
        }
        // SCAL(beta, p);
        // AXPY(alpha, residual, p);
        // XPEY(solution, p);
        std::function<void(I64)> taskFunction = 
            [&](I64 tp){
                I64 start_tp = startIndex[tp];
                I64 end_tp = endIndex[tp];
                for(I64 ti = start_tp; ti < end_tp; ++ ti){
                    p[ti] = beta * p[ti] + alpha * residual[ti];
                    solution[ti] += p[ti];
                }
            };
        switch(nestLevel){
            case 0: case 1: 
                threadManager.RunTask(nestLevel, threadTask, taskFunction);
                break;
            default:
                taskFunction(0);
                break;
        }
    }
}

} // namespace Ddpca