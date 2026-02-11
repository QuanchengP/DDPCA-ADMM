#include "PCGeometricMultigrid.hpp"
#include "../Mesh/Mesh.hpp"

#include <cassert>
#include <numeric>    // for std::reduce, std::transform_reduce

namespace Ddpca {

void PCGeometricMultigrid::Establish(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("PCGeometricMultigrid::Establish");
    const I64 maxLevel = realProlong.size();
    //
    realProlongT.resize(maxLevel);
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction_0 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tl = start_tp; tl < end_tp; ++ tl){
                TRANSPOSE(realProlong[tl], realProlongT[tl]);
            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_0);
            break;
        default:
            taskFunction_0(0);
            break;
    }
    // hierarchyStiffness.resize(maxLevel + 1); //must be deleted
    for(I64 tl = maxLevel - 1; tl >= 0; -- tl){
        SparseMatrix tempIK;
        GEMM(realProlongT[tl], hierarchyStiffness[tl + 1], tempIK, threadTask, nestLevel);
        GEMM(tempIK, realProlong[tl], hierarchyStiffness[tl], threadTask, nestLevel);
    }
    //
    coarsestDs.Establish(hierarchyStiffness[0]);
    //
    y.resize(maxLevel + 1); 
    residualError.resize(maxLevel + 1);
    coarseResiErro.resize(maxLevel + 1);
    coarseSolution.resize(maxLevel + 1);
    for(I64 tl = 1; tl <= maxLevel; ++ tl){
        y[tl].resize(realProlong[tl - 1].M);
        residualError[tl].resize(realProlong[tl - 1].M);
        coarseResiErro[tl].resize(realProlong[tl - 1].N);
        coarseSolution[tl].resize(realProlong[tl - 1].N);
    }

    EvenlyDistribute(maxLevel + 1, numbPart, startIndex, endIndex);
    switch(relaxation){
        case RelaxationSSGS:{
            gaussSeidel.resize(maxLevel + 1);
            std::function<void(I64)> taskFunction_1 = 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 tl = start_tp; tl < end_tp; ++ tl){
                        gaussSeidel[tl].Establish(hierarchyStiffness[tl]);
                    }
                };
            switch(nestLevel){
                case 0: case 1: 
                    threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
                    break;
                default:
                    taskFunction_1(0);
                    break;
            }
            break;
        }
        case RelaxationPSGS:{
            gaussSeidel.resize(maxLevel + 1);
            std::function<void(I64)> taskFunction_1 = 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 tl = start_tp; tl < end_tp; ++ tl){
                        gaussSeidel[tl].Establish(hierarchyStiffness[tl]);
                        gaussSeidel[tl].HierarchySplit(hierarchyStiffness[tl]);
                    }
                };
            switch(nestLevel){
                case 0: case 1: 
                    threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
                    break;
                default:
                    taskFunction_1(0);
                    break;
            }
            break;
        }
        case RelaxationJacobian:{
            jacobian.resize(maxLevel + 1);
            std::function<void(I64)> taskFunction_1 = 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 tl = start_tp; tl < end_tp; ++ tl){
                        jacobian[tl].Establish(hierarchyStiffness[tl]);
                    }
                };
            switch(nestLevel){
                case 0: case 1: 
                    threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
                    break;
                default:
                    taskFunction_1(0);
                    break;
            }
            break;
        }
        case RelaxationChebyshev:{
            chebRela.resize(maxLevel + 1);
            for(I64 tl = 1; tl <= maxLevel; ++ tl){
                chebRela[tl].Establish(hierarchyStiffness[tl], threadTask, nestLevel);
            }
            break;
        }
        case RelaxationSSOR:{
            ssor.resize(maxLevel + 1);
            std::function<void(I64)> taskFunction_1 = 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 tl = start_tp; tl < end_tp; ++ tl){
                        ssor[tl].Establish(hierarchyStiffness[tl]);
                    }
                };
            switch(nestLevel){
                case 0: case 1: 
                    threadManager.RunTask(nestLevel, threadTask, taskFunction_1);
                    break;
                default:
                    taskFunction_1(0);
                    break;
            }
            break;
        }
        default:
            assert(0 && "relaxation not implemented");
            break;
    }
    // conjGrad.InverseDiagonal(numbThreads);
}

void PCGeometricMultigrid::Apply(
    I64 tempLevel, 
    const AlignedVectorRx& rhs, 
    AlignedVectorRx& solution, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    if(tempLevel == 0){
        coarsestDs.Solve(rhs, solution/*, numbThreads*/);
        return;
    }
    
    //pre-smoothing
    switch(relaxation){
        case RelaxationSSGS:
            gaussSeidel[tempLevel].Apply(rhs, solution, y[tempLevel], threadTask, nestLevel);
            break;
        case RelaxationPSGS:
            gaussSeidel[tempLevel].Apply(rhs, solution, threadTask, nestLevel);
            break;
        case RelaxationJacobian:
            jacobian[tempLevel].Apply(hierarchyStiffness[tempLevel], rhs, solution, 4, threadTask, nestLevel);
            break;
        case RelaxationChebyshev:
            chebRela[tempLevel].Apply(hierarchyStiffness[tempLevel], rhs, solution, threadTask, nestLevel);
            break;
        case RelaxationSSOR:
            ssor[tempLevel].Apply(rhs, solution, threadTask, nestLevel);
            break;
        default:
            break;
    }

    //coarse grid correction
    std::fill(std::execution::unseq, coarseSolution[tempLevel].begin(), coarseSolution[tempLevel].end(), 0.0);
    std::copy(std::execution::unseq, rhs.begin(), rhs.end(), residualError[tempLevel].begin());
    switch(relaxation){
        case RelaxationSSGS:
            AXPY(-1.0, y[tempLevel], residualError[tempLevel]);
            break;
        case RelaxationPSGS: case RelaxationJacobian: case RelaxationChebyshev: case RelaxationSSOR:
            MV(-1.0, hierarchyStiffness[tempLevel], solution, 
                1.0, residualError[tempLevel], residualError[tempLevel], threadTask, nestLevel);
            break;
        default:
            break;
    }
    MV(realProlongT[tempLevel - 1], residualError[tempLevel], coarseResiErro[tempLevel], threadTask, nestLevel);
    Apply(tempLevel - 1, coarseResiErro[tempLevel], coarseSolution[tempLevel], threadTask, nestLevel);
    PEMV(realProlong[tempLevel - 1], coarseSolution[tempLevel], solution, threadTask, nestLevel);

    //post-smoothing
    switch(relaxation){
        case RelaxationSSGS:
            gaussSeidel[tempLevel].Apply(rhs, solution, y[tempLevel], threadTask, nestLevel);
            break;
        case RelaxationPSGS:
            gaussSeidel[tempLevel].Apply(rhs, solution, threadTask, nestLevel);
            break;
        case RelaxationJacobian:
            jacobian[tempLevel].Apply(hierarchyStiffness[tempLevel], rhs, solution, 4, threadTask, nestLevel);
            break;
        case RelaxationChebyshev:
            chebRela[tempLevel].Apply(hierarchyStiffness[tempLevel], rhs, solution, threadTask, nestLevel);
            break;
        case RelaxationSSOR:
            ssor[tempLevel].Apply(rhs, solution, threadTask, nestLevel);
            break;
        default:
            break;
    }
}

} //namespace Ddpca