#include "FemStiffnessMatrix.hpp"

#include <ranges>

namespace Ddpca{

void FemStiffnessMatrix(
    Mesh& inputMesh, 
    const Real materialElasticity, 
    const Real materialPoisson, 
    SparseMatrix& outputMatrix, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("FemStiffnessMatrix");
    //
    DenseMatrix elasticity(6, 6);
    CalculateElasticity(materialElasticity, materialPoisson, elasticity);
    //
    I64 M = 3 * inputMesh.node2Coordinate.size();
    I64 N = M;
    I64 nnz = inputMesh.EffectiveElements(threadTask, nestLevel) * 576;
    Log("    FemStiffnessMatrix: M   = " + std::to_string(M));
    Log("    FemStiffnessMatrix: nnz = " + std::to_string(nnz));
    Triplet* coo = AlignedAllocate<Triplet>(nnz);
    //
    I64 elementSize = inputMesh.elements.size();
    I64 numbPart = (nestLevel < threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(elementSize, numbPart, startIndex, endIndex);
    //
    I64 baseJKMN[8][3][8][3];
    for(I64 tj = 0; tj < 8; tj ++){
        I64 base_j = tj * 72;
        for(I64 tk = 0; tk < 3; tk ++){
            I64 base_k = base_j + tk * 24;
            for(I64 tm = 0; tm < 8; tm ++){
                I64 base_m = base_k + tm * 3;
                for(I64 tn = 0; tn  < 3; tn ++){
                    baseJKMN[tj][tk][tm][tn] = base_m + tn;
                }
            }
        }
    }
    TrilinearQuadrature<3> trilQuad = GetTrilinearQuadrature<3>();
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            DenseMatrix exyz(8, 3), shapeDerivatives(3, 8), shapeDerivativesX(3, 8);
            DenseMatrix jacobianJ(3, 3), invJacobianJ(3, 3), K_E(24, 24), B_L(6, 24), B_LTE(24, 6);
            Real jacobianDeterminant;
            auto iota_view = std::views::iota(0, 8);
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                if(ti % 10000 == 0) std::cout << "FemStiffnessMatrix ti = " << ti << "\n";
                const auto& element = inputMesh.elements[ti];
                if(element.children.size() > 0){
                    continue;
                }
                for(I64 tj = 0; tj < 8; ++ tj){
                    auto iteratorImnc = (inputMesh.node2Coordinate).find(
                        element.cornerNodes[tj]);
                    for(I64 tk = 0; tk < 3; ++ tk){
                        exyz(tj,tk) = (iteratorImnc->second)[tk];
                    }
                }
                //
                K_E.Fill(0.0);
                for(I64 tj = 0; tj < trilQuad.numbGaussPoints; ++ tj){
                    shapeDerivatives.data = (trilQuad.shapeDerivatives)[tj].data;
                    GEMM(shapeDerivatives, exyz, jacobianJ);
                    //a11(a22a33-a23a32)-a12(a21a33-a23a31)+a13(a21a32-a22a31)
                    jacobianDeterminant = 
                        jacobianJ(0,0) * (jacobianJ(1,1) * jacobianJ(2,2) - jacobianJ(1,2) * jacobianJ(2,1)) -
                        jacobianJ(0,1) * (jacobianJ(1,0) * jacobianJ(2,2) - jacobianJ(1,2) * jacobianJ(2,0)) +
                        jacobianJ(0,2) * (jacobianJ(1,0) * jacobianJ(2,1) - jacobianJ(1,1) * jacobianJ(2,0));
                    invJacobianJ(0,0) =  jacobianJ(1,1) * jacobianJ(2,2) - jacobianJ(1,2) * jacobianJ(2,1);
                    invJacobianJ(0,1) = -jacobianJ(0,1) * jacobianJ(2,2) + jacobianJ(0,2) * jacobianJ(2,1);
                    invJacobianJ(0,2) =  jacobianJ(0,1) * jacobianJ(1,2) - jacobianJ(0,2) * jacobianJ(1,1);
                    invJacobianJ(1,0) = -jacobianJ(1,0) * jacobianJ(2,2) + jacobianJ(1,2) * jacobianJ(2,0);
                    invJacobianJ(1,1) =  jacobianJ(0,0) * jacobianJ(2,2) - jacobianJ(0,2) * jacobianJ(2,0);
                    invJacobianJ(1,2) = -jacobianJ(0,0) * jacobianJ(1,2) + jacobianJ(0,2) * jacobianJ(1,0);
                    invJacobianJ(2,0) =  jacobianJ(1,0) * jacobianJ(2,1) - jacobianJ(1,1) * jacobianJ(2,0);
                    invJacobianJ(2,1) = -jacobianJ(0,0) * jacobianJ(2,1) + jacobianJ(0,1) * jacobianJ(2,0);
                    invJacobianJ(2,2) =  jacobianJ(0,0) * jacobianJ(1,1) - jacobianJ(0,1) * jacobianJ(1,0);
                    SCAL(1.0 / jacobianDeterminant, invJacobianJ);
                    GEMM(invJacobianJ, shapeDerivatives, shapeDerivativesX);
                    std::for_each(std::execution::unseq, iota_view.begin(), iota_view.end(), [&](I64 tk) {
                        B_L(0, 3 * tk + 0) = shapeDerivativesX(0, tk);
                        B_L(1, 3 * tk + 1) = shapeDerivativesX(1, tk);
                        B_L(2, 3 * tk + 2) = shapeDerivativesX(2, tk);
                        B_L(3, 3 * tk + 0) = shapeDerivativesX(1, tk);
                        B_L(3, 3 * tk + 1) = shapeDerivativesX(0, tk);
                        B_L(4, 3 * tk + 1) = shapeDerivativesX(2, tk);
                        B_L(4, 3 * tk + 2) = shapeDerivativesX(1, tk);
                        B_L(5, 3 * tk + 0) = shapeDerivativesX(2, tk);
                        B_L(5, 3 * tk + 2) = shapeDerivativesX(0, tk);
                    });
                    GEMTM(B_L, elasticity, B_LTE);
                    GEPEMM(trilQuad.weights[tj] * jacobianDeterminant, B_LTE, B_L, K_E);
                }
                //
                I64 elementBase = inputMesh.effePrefSum[ti] * 576;
                for(I64 tj = 0; tj < 8; ++ tj){
                    I64 row_j = 3 * element.cornerNodes[tj];
                    for(I64 tk = 0; tk < 3; ++ tk){
                        I64 row_jk = row_j + tk;
                        for(I64 tm = 0; tm < 8; ++ tm){
                            I64 col_m = 3 * element.cornerNodes[tm];
                            for(I64 tn = 0; tn  < 3; ++ tn){
                                I64 base_n = baseJKMN[tj][tk][tm][tn];
                                I64 col_mn = col_m + tn;
                                I64 ebPbn = elementBase + base_n;
                                coo[ebPbn].row = row_jk;
                                coo[ebPbn].col = col_mn;
                                coo[ebPbn].val = K_E(3 * tj + tk, 3 * tm + tn);
                            }
                        }
                    }
                }
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
    outputMatrix.Coo2Csr(M, N, coo, nnz, threadTask, nestLevel);
    Deallocate<Triplet>(coo, nnz);
}

void CalculateElasticity(
    const Real materialElasticity, 
    const Real materialPoisson, 
    DenseMatrix& outputElasticity){
    //
    assert((outputElasticity.rows == 6 && outputElasticity.cols == 6) 
        && "Elasticity matrix must be 6x6");

    Real lambda = materialElasticity * materialPoisson 
        / (1 + materialPoisson) / (1 - 2 * materialPoisson);
    Real mu = materialElasticity / (2 + 2 * materialPoisson);
    //symmetric
    outputElasticity.data = {
		(Real)2.0 * mu + lambda, lambda, lambda, 0.0, 0.0, 0.0,
		lambda, (Real)2.0 * mu + lambda, lambda, 0.0, 0.0, 0.0,
		lambda, lambda, (Real)2.0 * mu + lambda, 0.0, 0.0, 0.0,
    	0.0, 0.0, 0.0, mu, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, mu, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, mu};
}

}