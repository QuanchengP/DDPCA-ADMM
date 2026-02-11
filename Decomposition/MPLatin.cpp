#include "MPLatin.hpp"

namespace Ddpca{

void MPLatin::SubTransfer(){
    Log("MPLatin::SubTransfer");
    //
    const I64 inteSize = (*interfaces).size();
    contactUserSolver.resize(inteSize);
    accuContProl.resize(inteSize);
    accuContProlT.resize(inteSize);
    I64 tv_master = 0;
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const std::map<I64, I64>& tempNici = (*nodeId2ContactId)[ts][tv_master];
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 tempDomaId = tempInterface.domainIndex[tv_master];
        const SingleDomain& tempDomain = (*domains)[tempDomaId];
        const I64 numbNode = tempDomain.geomMult.nodeLevelPosition.size();
        const I64 maxLevel = tempDomain.mesh.maxLevel;
        //
        I64 M_0 = tempNici.size();
        I64 N_0 = numbNode;
        I64 nnz_0 = M_0;
        Triplet* coo_0 = AlignedAllocate<Triplet>(nnz_0);
        for(I64 ti = 0; ti < numbNode; ++ ti){
            I64 tempNode = tempDomain.geomMult.positionNode[ti];
            auto itNici = tempNici.find(tempNode);
            if(itNici != tempNici.end()){
                coo_0[itNici->second].row = itNici->second;
                coo_0[itNici->second].col = ti;
                coo_0[itNici->second].val = 1.0;
            }
        }
        contactUserSolver[ts].Coo2Csr(M_0, N_0, coo_0, nnz_0, threadManager.interfaceS2M[ts], 1);
        Deallocate<Triplet>(coo_0, nnz_0);
        //
        SparseMatrix tempAccuContProl = contactUserSolver[ts];
        // hanging nodes are not allowed to be contact nodes.
        for(I64 tl = maxLevel; tl >= realDomaLeve[tempDomaId]; -- tl){
            SparseMatrix tempProlong = std::move(tempAccuContProl);
            GEMM(tempProlong, tempDomain.geomMult.scalarProlong[tl], tempAccuContProl, 
                threadManager.interfaceS2M[ts], 1);
        }
        std::vector<I64> usedColumn(tempAccuContProl.N + 1, 0);
        for(I64 ti = 0; ti < tempAccuContProl.M; ++ ti){
            I64 tj_start = tempAccuContProl.row_ptr[ti];
            I64 tj_end = tempAccuContProl.row_ptr[ti + 1];
            for(I64 tj = tj_start; tj < tj_end; ++ tj){
                usedColumn[tempAccuContProl.col_ind[tj]] = 1;
            }
        }
        std::exclusive_scan(
            usedColumn.begin(), usedColumn.end(), usedColumn.begin(), 0);
        //
        if(tempInterface.frictionCoefficient == 0.0){
            I64 M = M_0;
            I64 N = usedColumn.back();
            I64 nnz = tempAccuContProl.nnz;
            Triplet* coo = AlignedAllocate<Triplet>(nnz);
            for(I64 ti = 0; ti < tempAccuContProl.M; ++ ti){
                I64 tj_start = tempAccuContProl.row_ptr[ti];
                I64 tj_end = tempAccuContProl.row_ptr[ti + 1];
                for(I64 tj = tj_start; tj < tj_end; ++ tj){
                    I64 col_tj = tempAccuContProl.col_ind[tj];
                    coo[tj].row = ti;
                    coo[tj].col = usedColumn[col_tj];
                    coo[tj].val = tempAccuContProl.val[tj];
                }
            }
            accuContProl[ts].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
            Deallocate<Triplet>(coo, nnz);
        }
        else{
            I64 M = 3 * M_0;
            I64 N = 3 * usedColumn.back();
            I64 nnz = 3 * tempAccuContProl.nnz;
            Triplet* coo = AlignedAllocate<Triplet>(nnz);
            for(I64 ti = 0; ti < tempAccuContProl.M; ++ ti){
                I64 tj_start = tempAccuContProl.row_ptr[ti];
                I64 tj_end = tempAccuContProl.row_ptr[ti + 1];
                for(I64 tj = tj_start; tj < tj_end; ++ tj){
                    I64 col_tj = tempAccuContProl.col_ind[tj];
                    for(I64 tk = 0; tk < 3; ++ tk){
                        const I64 tjtk = 3 * tj + tk;
                        coo[tjtk].row = 3 * ti + tk;
                        coo[tjtk].col = 3 * usedColumn[col_tj] + tk;
                        coo[tjtk].val = tempAccuContProl.val[tj];
                    }
                }
            }
            accuContProl[ts].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
            Deallocate<Triplet>(coo, nnz);
        }
        TRANSPOSE(accuContProl[ts], accuContProlT[ts]);
    });
}

void MPLatin::SubAccuDomaProl(){
    Log("MPLatin::SubAccuDomaProl");
    //
    const I64 domaSize = (*domains).size();
    accuDomaProl.resize(domaSize);
    accuDomaProlT.resize(domaSize);
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        const SingleDomain& tempDomain = (*domains)[td];
        accuDomaProl[td] = tempDomain.geomMult.rom;
        const I64 maxLevel = tempDomain.mesh.maxLevel;
        for(I64 tl = maxLevel - 1; tl >= realDomaLeve[td]; -- tl){
            SparseMatrix tempMatrix_0 = std::move(accuDomaProl[td]);
            GEMM(tempMatrix_0, tempDomain.conjGrad.pcGM.realProlong[tl], 
                accuDomaProl[td], threadManager.domainS2M[td], 1);
        }
        TRANSPOSE(accuDomaProl[td], accuDomaProlT[td]);
    });
}

void MPLatin::SubDomaAuxi(){
    Log("MPLatin::SubDomaAuxi");
    //
    const I64 inteSize = (*interfaces).size();
    domaAuxi.resize(inteSize);
    I64 tv_master = 0;
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        const std::map<I64, I64>& tempNici_master = (*nodeId2ContactId)[ts][tv_master];
        for(I64 tv = 0; tv < 2; ++ tv){
            const I64 tempDomaId = tempInterface.domainIndex[tv];
            const SingleDomain& tempDomain = (*domains)[tempDomaId];
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = 3 * tempDomain.mesh.coordinate2Node.size();
                I64 N = tempNici_master.size();
                I64 nnz = 48 * inpoSize; // 12 * 4
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(12, 4), N_e(3, 12), normVect(1, 3), M_0e(1, 4), tempMatr(3, 4);
                std::array<Real, 4> M_e;
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.begin());
                    // col-major
                    N_e.Fill(
                        M_e[0], 0.0, 0.0,
                        0.0, M_e[0], 0.0, 
                        0.0, 0.0, M_e[0], 
                        M_e[1], 0.0, 0.0,
                        0.0, M_e[1], 0.0, 
                        0.0, 0.0, M_e[1], 
                        M_e[2], 0.0, 0.0,
                        0.0, M_e[2], 0.0, 
                        0.0, 0.0, M_e[2], 
                        M_e[3], 0.0, 0.0,
                        0.0, M_e[3], 0.0, 
                        0.0, 0.0, M_e[3]);
                    std::copy(std::execution::unseq, 
                        tempInpo.basisVector[0].begin(), tempInpo.basisVector[0].end(), 
                        normVect.data.begin());
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(),  
                        M_0e.data.begin());
                    GEMTM(normVect, M_0e, tempMatr);
                    GEMTM(N_e, tempMatr, elemMatr);
                    SCAL(tempInpo.quadratureWeight * tempInterface.normPenaPara, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
                            for(I64 tm = 0; tm < 4; tm ++){
                                auto itNici_master = tempNici_master.find(tempInpo.node[tv_master][tm]);
                                coo[cooIndex].row = free_tk;
                                coo[cooIndex].col = itNici_master->second;
                                coo[cooIndex].val = elemMatr(3 * tj + tk, tm);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                domaAuxi[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            else{
                I64 M = 3 * tempDomain.mesh.coordinate2Node.size();
                I64 N = 3 * tempNici_master.size();
                I64 nnz = 144 * inpoSize; // 12 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(12, 12), N_e(3, 12), normTang(3, 3), M_0e(1, 4), N_0e(3, 12), 
                    tempPena(3, 3), tempMatr_0(3, 12), tempMatr_1(3, 12), tempMatr_2(3, 12);
                std::array<Real, 4> M_e;
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.begin());
                    // col-major
                    N_e.Fill(
                        M_e[0], 0.0, 0.0,
                        0.0, M_e[0], 0.0, 
                        0.0, 0.0, M_e[0], 
                        M_e[1], 0.0, 0.0,
                        0.0, M_e[1], 0.0, 
                        0.0, 0.0, M_e[1], 
                        M_e[2], 0.0, 0.0,
                        0.0, M_e[2], 0.0, 
                        0.0, 0.0, M_e[2], 
                        M_e[3], 0.0, 0.0,
                        0.0, M_e[3], 0.0, 
                        0.0, 0.0, M_e[3]);
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                        M_0e.data.begin());
                    // col-major
                    N_0e.Fill(
                        M_0e(0,0), 0.0, 0.0,
                        0.0, M_0e(0,0), 0.0, 
                        0.0, 0.0, M_0e(0,0), 
                        M_0e(0,1), 0.0, 0.0,
                        0.0, M_0e(0,1), 0.0, 
                        0.0, 0.0, M_0e(0,1), 
                        M_0e(0,2), 0.0, 0.0,
                        0.0, M_0e(0,2), 0.0, 
                        0.0, 0.0, M_0e(0,2), 
                        M_0e(0,3), 0.0, 0.0,
                        0.0, M_0e(0,3), 0.0, 
                        0.0, 0.0, M_0e(0,3));
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(normTang, N_0e, tempMatr_0);
                    GEMM(tempPena, tempMatr_0, tempMatr_1);
                    GEMTM(normTang, tempMatr_1, tempMatr_2);
                    GEMTM(N_e, tempMatr_2, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
                            for(I64 tm = 0; tm < 4; ++ tm){
                                auto itNici_master = 
                                    tempNici_master.find(tempInpo.node[tv_master][tm]);
                                for(I64 tn = 0; tn < 3; ++ tn){
                                    coo[cooIndex].row = free_tk;
                                    coo[cooIndex].col = 3 * itNici_master->second + tn;
                                    coo[cooIndex].val = elemMatr(3 * tj + tk, 3 * tm + tn);
                                    ++ cooIndex;
                                }
                            }
                        }
                    }
                }
                domaAuxi[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            //
            SparseMatrix tempMatrix_0 = std::move(domaAuxi[ts][tv]);
            GEMM(tempMatrix_0, accuContProl[ts], 
                domaAuxi[ts][tv], threadManager.interfaceS2M[ts], 1);
            SparseMatrix tempMatrix_1 = std::move(domaAuxi[ts][tv]);
            GEMM(accuDomaProlT[tempDomaId], tempMatrix_1, 
                domaAuxi[ts][tv], threadManager.interfaceS2M[ts], 1);
        }
    });
}

void MPLatin::SubAuxiAuxi(){
    //
    Log("MPLatin::SubAuxiAuxi");
    const I64 inteSize = (*interfaces).size();
    auxiAuxi.resize(inteSize);
    I64 tv_master = 0;
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        const std::map<I64, I64>& tempNici_master = (*nodeId2ContactId)[ts][tv_master];
        if(tempInterface.frictionCoefficient == 0.0){
            I64 M = tempNici_master.size();
            I64 N = tempNici_master.size();
            I64 nnz = 16 * inpoSize; // 4 * 4
            Triplet* coo = AlignedAllocate<Triplet>(nnz);
            I64 cooIndex = 0;
            DenseMatrix elemMatr(4, 4);
            DenseMatrix M_0e(1, 4);
            for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                std::copy(std::execution::unseq, 
                    tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                    M_0e.data.begin());
                GEMTM(M_0e, M_0e, elemMatr);
                SCAL(tempInpo.quadratureWeight * tempInterface.normPenaPara, elemMatr);
                for(I64 tj = 0; tj < 4; ++ tj){
                    auto itNici_tj = tempNici_master.find(tempInpo.node[tv_master][tj]);
                    I64 row_tj = itNici_tj->second;
                    for(I64 tk = 0; tk < 4; ++ tk){
                        auto itNici_tk = tempNici_master.find(tempInpo.node[tv_master][tk]);
                        I64 col_tk = itNici_tk->second;
                        coo[cooIndex].row = row_tj;
                        coo[cooIndex].col = col_tk;
                        coo[cooIndex].val = elemMatr(tj, tk);
                        ++ cooIndex;
                    }
                }
            }
            auxiAuxi[ts].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
            Deallocate<Triplet>(coo, nnz);
        }
        else{
            I64 M = 3 * tempNici_master.size();
            I64 N = 3 * tempNici_master.size();
            I64 nnz = 144 * inpoSize; // 12 * 12
            Triplet* coo = AlignedAllocate<Triplet>(nnz);
            I64 cooIndex = 0;
            DenseMatrix elemMatr(12, 12), normTang(3, 3), M_0e(1, 4), N_0e(3, 12), 
                tempPena(3, 3), tempMatr_0(3, 12), tempMatr_1(3, 12), tempMatr_2(3, 12);
            for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                normTang.Fill(
                    tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                    tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                    tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                std::copy(std::execution::unseq, 
                    tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                    M_0e.data.begin());
                // col-major
                N_0e.Fill(
                    M_0e(0,0), 0.0, 0.0,
                    0.0, M_0e(0,0), 0.0, 
                    0.0, 0.0, M_0e(0,0), 
                    M_0e(0,1), 0.0, 0.0,
                    0.0, M_0e(0,1), 0.0, 
                    0.0, 0.0, M_0e(0,1), 
                    M_0e(0,2), 0.0, 0.0,
                    0.0, M_0e(0,2), 0.0, 
                    0.0, 0.0, M_0e(0,2), 
                    M_0e(0,3), 0.0, 0.0,
                    0.0, M_0e(0,3), 0.0, 
                    0.0, 0.0, M_0e(0,3));
                tempPena.Fill(
                    tempInterface.normPenaPara, 0.0, 0.0, 
                    0.0, tempInterface.tangPenaPara, 0.0, 
                    0.0, 0.0, tempInterface.tangPenaPara);
                GEMM(normTang, N_0e, tempMatr_0);
                GEMM(tempPena, tempMatr_0, tempMatr_1);
                GEMTM(normTang, tempMatr_1, tempMatr_2);
                GEMTM(N_0e, tempMatr_2, elemMatr);
                SCAL(tempInpo.quadratureWeight, elemMatr);
                for(I64 tj = 0; tj < 4; ++ tj){
                    auto itNici_tj = tempNici_master.find(tempInpo.node[tv_master][tj]);
                    I64 row_tj = itNici_tj->second;
                    for(I64 tk = 0; tk  < 3; ++ tk){
                        I64 free_tk = 3 * row_tj + tk;
                        for(I64 tm = 0; tm < 4; ++ tm){
                            auto itNici_tm = tempNici_master.find(tempInpo.node[tv_master][tm]);
                            I64 col_tm = itNici_tm->second;
                            for(I64 tn = 0; tn < 3; ++ tn){
                                coo[cooIndex].row = free_tk;
                                coo[cooIndex].col = 3 * col_tm + tn;
                                coo[cooIndex].val = elemMatr(3 * tj + tk, 3 * tm + tn);
                                ++ cooIndex;
                            }
                        }
                    }
                }
            }
            auxiAuxi[ts].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
            Deallocate<Triplet>(coo, nnz);
        }
        SparseMatrix tempMatrix = std::move(auxiAuxi[ts]);
        GEMM(tempMatrix, accuContProl[ts], auxiAuxi[ts], threadManager.interfaceS2M[ts], 1);
        tempMatrix = std::move(auxiAuxi[ts]);
        GEMM(accuContProlT[ts], tempMatrix, auxiAuxi[ts], threadManager.interfaceS2M[ts], 1);
    });
}

void MPLatin::SubGlobalCouple(){
    Log("MPLatin::SubGlobalCouple");
    //
    I64 M = 0;
    I64 nnz = 0;
    const I64 domaSize = (*domains).size();
    const I64 inteSize = (*interfaces).size();
    MAccu.assign(domaSize + inteSize + 1, 0);
    std::vector<I64> nnzAccu(domaSize + inteSize + 1, 0);
    for(I64 td = 0; td < domaSize; ++ td){
        const SingleDomain& tempDomain = (*domains)[td];
        const SparseMatrix& tempStiffness = 
            tempDomain.conjGrad.pcGM.hierarchyStiffness[realDomaLeve[td]];
        M += tempStiffness.M;
        nnz += tempStiffness.nnz;
        MAccu[td + 1] = MAccu[td] + tempStiffness.M;
        nnzAccu[td + 1] = nnzAccu[td] + tempStiffness.nnz;
    }
    for(I64 ts = 0; ts < inteSize; ++ ts){
        I64 tempM = auxiAuxi[ts].M;
        I64 tempNnz = 2 * (domaAuxi[ts][0].nnz + domaAuxi[ts][1].nnz) + auxiAuxi[ts].nnz;
        M += tempM;
        nnz += tempNnz;
        MAccu[domaSize + ts + 1] = MAccu[domaSize + ts] + tempM;
        nnzAccu[domaSize + ts + 1] = nnzAccu[domaSize + ts] + tempNnz;
    }
    //
    I64 N = M;
    Triplet* coo = AlignedAllocate<Triplet>(nnz);
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        const SingleDomain& tempDomain = (*domains)[td];
        const SparseMatrix& tempStiffness = 
            tempDomain.conjGrad.pcGM.hierarchyStiffness[realDomaLeve[td]];
        const I64 tempMAccu = MAccu[td];
        const I64 tempNnzAccu = nnzAccu[td];
        for(I64 ti = 0; ti < tempStiffness.M; ++ ti){
            I64 row = tempMAccu + ti;
            I64 tj_start = tempStiffness.row_ptr[ti];
            I64 tj_end = tempStiffness.row_ptr[ti + 1];
            for(I64 tj = tj_start; tj < tj_end; ++ tj){
                I64 col = tempMAccu + tempStiffness.col_ind[tj];
                const I64 tempNnz = tempNnzAccu + tj;
                coo[tempNnz].row = row;
                coo[tempNnz].col = col;
                coo[tempNnz].val = tempStiffness.val[tj];
            }
        }
    });
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 tempMAccu = MAccu[domaSize + ts];
        const I64 tempNnzAccu = nnzAccu[domaSize + ts];
        I64 cooIndex = 0;
        for(I64 tv = 0; tv < 2; ++ tv){
            const I64 tempDomaId = tempInterface.domainIndex[tv];
            const I64 tempMAccu_r = MAccu[tempDomaId];
            const SparseMatrix& tempDispAuxi = domaAuxi[ts][tv];
            for(I64 ti = 0; ti < tempDispAuxi.M; ++ ti){
                I64 row = tempMAccu_r + ti;
                I64 tj_start = tempDispAuxi.row_ptr[ti];
                I64 tj_end = tempDispAuxi.row_ptr[ti + 1];
                for(I64 tj = tj_start; tj < tj_end; ++ tj){
                    I64 col = tempMAccu + tempDispAuxi.col_ind[tj];
                    I64 realNnz = tempNnzAccu + cooIndex;
                    coo[realNnz].row = row;
                    coo[realNnz].col = col;
                    coo[realNnz].val = - tempDispAuxi.val[tj];
                    ++ cooIndex;
                    realNnz = tempNnzAccu + cooIndex;
                    coo[realNnz].row = col;
                    coo[realNnz].col = row;
                    coo[realNnz].val = - tempDispAuxi.val[tj];
                    ++ cooIndex;
                }
            }
        }
        const SparseMatrix& tempAuxiAuxi = auxiAuxi[ts];
        for(I64 ti = 0; ti < tempAuxiAuxi.M; ++ ti){
            I64 row = tempMAccu + ti;
            I64 tj_start = tempAuxiAuxi.row_ptr[ti];
            I64 tj_end = tempAuxiAuxi.row_ptr[ti + 1];
            for(I64 tj = tj_start; tj < tj_end; ++ tj){
                I64 col = tempMAccu + tempAuxiAuxi.col_ind[tj];
                I64 realNnz = tempNnzAccu + cooIndex;
                coo[realNnz].row = row;
                coo[realNnz].col = col;
                coo[realNnz].val = 2.0 * tempAuxiAuxi.val[tj];
                ++ cooIndex;
            }
        }
    });
    //
    globalCouple.Coo2Csr(M, N, coo, nnz, threadManager.one2oneS2S, 0);
    Deallocate<Triplet>(coo, nnz);
    macroDs.Establish(globalCouple);
}

void MPLatin::SubErrorMA(){
    Log("MPLatin::SubErrorMA");
    //
    const I64 inteSize = (*interfaces).size();
    erroMult.resize(inteSize);
    erroPenaAuxi.resize(inteSize);
    I64 tv_master = 0;
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        const std::map<I64, I64>& tempNici_master = (*nodeId2ContactId)[ts][tv_master];
        for(I64 tv = 0; tv < 2; ++ tv){
            const std::map<I64, I64>& tempNici = (*nodeId2ContactId)[ts][tv];
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = tempNici_master.size();
                I64 N = tempNici.size();
                I64 nnz = 16 * inpoSize; // 4 * 4
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* coo_1 = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(4, 4), M_e(1, 4), M_0e(1, 4);
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq,    
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.data.begin());
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                        M_0e.data.begin());
                    GEMTM(M_0e, M_e, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        auto itNici_master = tempNici_master.find(tempInpo.node[tv_master][tj]);
                        I64 row_tj = itNici_master->second;
                        for(I64 tk = 0; tk < 4; ++ tk){
                            auto itNici = tempNici.find(tempInpo.node[tv][tk]);
                            I64 col_tk = itNici->second;
                            coo[cooIndex].row = row_tj;
                            coo[cooIndex].col = col_tk;
                            coo[cooIndex].val = elemMatr(tj, tk);
                            coo_1[cooIndex].row = row_tj;
                            coo_1[cooIndex].col = col_tk;
                            coo_1[cooIndex].val = tempInterface.normPenaPara * elemMatr(tj, tk);
                            ++ cooIndex;
                        }
                    }
                }
                erroMult[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                erroPenaAuxi[ts][tv].Coo2Csr(M, N, coo_1, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo_1, nnz);
            }
            else{
                I64 M = 3 * tempNici_master.size();
                I64 N = 3 * tempNici.size();
                I64 nnz = 144 * inpoSize; // 12 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* coo_1 = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(12, 12), elemMatr_1(12, 12), M_e(1, 4), N_e(3, 12), 
                    normTang(3, 3), M_0e(1, 4), N_0e(3, 12), tempPena(3, 3), 
                    tempMatr_0(3, 12), tempMatr_1(3, 12), tempMatr_2(3, 12), tempMatr_3(3, 12);
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.data.begin());
                    // col-major
                    N_e.Fill(
                        M_e(0,0), 0.0, 0.0,
                        0.0, M_e(0,0), 0.0, 
                        0.0, 0.0, M_e(0,0), 
                        M_e(0,1), 0.0, 0.0,
                        0.0, M_e(0,1), 0.0, 
                        0.0, 0.0, M_e(0,1), 
                        M_e(0,2), 0.0, 0.0,
                        0.0, M_e(0,2), 0.0, 
                        0.0, 0.0, M_e(0,2), 
                        M_e(0,3), 0.0, 0.0,
                        0.0, M_e(0,3), 0.0, 
                        0.0, 0.0, M_e(0,3));
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                        M_0e.data.begin());
                    // col-major
                    N_0e.Fill(
                        M_0e(0,0), 0.0, 0.0,
                        0.0, M_0e(0,0), 0.0, 
                        0.0, 0.0, M_0e(0,0), 
                        M_0e(0,1), 0.0, 0.0,
                        0.0, M_0e(0,1), 0.0, 
                        0.0, 0.0, M_0e(0,1), 
                        M_0e(0,2), 0.0, 0.0,
                        0.0, M_0e(0,2), 0.0, 
                        0.0, 0.0, M_0e(0,2), 
                        M_0e(0,3), 0.0, 0.0,
                        0.0, M_0e(0,3), 0.0, 
                        0.0, 0.0, M_0e(0,3));
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMTM(normTang, tempMatr_0, tempMatr_1);
                    GEMTM(N_0e, tempMatr_1, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    GEMM(tempPena, tempMatr_0, tempMatr_2);
                    GEMTM(normTang, tempMatr_2, tempMatr_3);
                    GEMTM(N_0e, tempMatr_3, elemMatr_1);
                    SCAL(tempInpo.quadratureWeight, elemMatr_1);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        auto itNici_master = tempNici_master.find(tempInpo.node[tv_master][tj]);
                        I64 row_tj = itNici_master->second;
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * row_tj + tk;
                            for(I64 tm = 0; tm < 4; ++ tm){
                                auto itNici = tempNici.find(tempInpo.node[tv][tm]);
                                I64 col_tm = itNici->second;
                                for(I64 tn = 0; tn < 3; ++ tn){
                                    coo[cooIndex].row = free_tk;
                                    coo[cooIndex].col = 3 * col_tm + tn;
                                    coo[cooIndex].val = elemMatr(3 * tj + tk, 3 * tm + tn);
                                    coo_1[cooIndex].row = free_tk;
                                    coo_1[cooIndex].col = 3 * col_tm + tn;
                                    coo_1[cooIndex].val = elemMatr_1(3 * tj + tk, 3 * tm + tn);
                                    ++ cooIndex;
                                }
                            }
                        }
                    }
                }
                erroMult[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                erroPenaAuxi[ts][tv].Coo2Csr(M, N, coo_1, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo_1, nnz);
            }
            SparseMatrix tempMatr = std::move(erroMult[ts][tv]);
            GEMM(accuContProlT[ts], tempMatr, erroMult[ts][tv], threadManager.interfaceS2M[ts], 1);
            SparseMatrix tempMatr_1 = std::move(erroPenaAuxi[ts][tv]);
            GEMM(accuContProlT[ts], tempMatr_1, erroPenaAuxi[ts][tv], threadManager.interfaceS2M[ts], 1);
        }
    });
}

void MPLatin::SubErroPenaDisp(){
    Log("MPLatin::SubErroPenaDisp");
    //
    const I64 inteSize = (*interfaces).size();
    erroPenaDisp.resize(inteSize);
    I64 tv_master = 0;
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = (*interfaces)[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        const std::map<I64, I64>& tempNici_master = (*nodeId2ContactId)[ts][tv_master];
        for(I64 tv = 0; tv < 2; ++ tv){
            const I64 tempDomaId = tempInterface.domainIndex[tv];
            const SingleDomain& tempDomain = (*domains)[tempDomaId];
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = tempNici_master.size();
                I64 N = 3 * tempDomain.mesh.coordinate2Node.size();
                I64 nnz = 48 * inpoSize; // 4 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(4, 12), M_e(1, 4), N_e(3, 12), normVect(1, 3), 
                    M_0e(1, 4), tempMatr_0(1, 12);
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.data.begin());
                    // col-major
                    N_e.Fill(
                        M_e(0,0), 0.0, 0.0,
                        0.0, M_e(0,0), 0.0, 
                        0.0, 0.0, M_e(0,0), 
                        M_e(0,1), 0.0, 0.0,
                        0.0, M_e(0,1), 0.0, 
                        0.0, 0.0, M_e(0,1), 
                        M_e(0,2), 0.0, 0.0,
                        0.0, M_e(0,2), 0.0, 
                        0.0, 0.0, M_e(0,2), 
                        M_e(0,3), 0.0, 0.0,
                        0.0, M_e(0,3), 0.0, 
                        0.0, 0.0, M_e(0,3));
                    std::copy(std::execution::unseq, 
                        tempInpo.basisVector[0].begin(), tempInpo.basisVector[0].end(), 
                        normVect.data.begin());
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                        M_0e.data.begin());
                    GEMM(normVect, N_e, tempMatr_0);
                    GEMTM(M_0e, tempMatr_0, elemMatr);
                    SCAL(tempInpo.quadratureWeight * tempInterface.normPenaPara, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        auto itNici_master = tempNici_master.find(tempInpo.node[tv_master][tj]);
                        I64 row_tj = itNici_master->second;
                        for(I64 tk = 0; tk < 4; ++ tk){
                            I64 node_tk = tempInpo.node[tv][tk];
                            for(I64 tm = 0; tm < 3; ++ tm){
                                I64 free_tm = 3 * node_tk + tm;
                                coo[cooIndex].row = row_tj;
                                coo[cooIndex].col = free_tm;
                                coo[cooIndex].val = elemMatr(tj, 3 * tk + tm);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                erroPenaDisp[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            else{
                I64 M = 3 * tempNici_master.size();
                I64 N = 3 * tempDomain.mesh.coordinate2Node.size();
                I64 nnz = 144 * inpoSize; // 12 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(12, 12), M_e(1, 4), N_e(3, 12), normTang(3, 3), 
                    M_0e(1, 4), N_0e(3, 12), tempPena(3, 3), tempMatr_0(3, 12), 
                    tempMatr_1(3, 12), tempMatr_2(3, 12);
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.data.begin());
                    // col-major
                    N_e.Fill(
                        M_e(0,0), 0.0, 0.0,
                        0.0, M_e(0,0), 0.0, 
                        0.0, 0.0, M_e(0,0), 
                        M_e(0,1), 0.0, 0.0,
                        0.0, M_e(0,1), 0.0, 
                        0.0, 0.0, M_e(0,1), 
                        M_e(0,2), 0.0, 0.0,
                        0.0, M_e(0,2), 0.0, 
                        0.0, 0.0, M_e(0,2), 
                        M_e(0,3), 0.0, 0.0,
                        0.0, M_e(0,3), 0.0, 
                        0.0, 0.0, M_e(0,3));
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv_master].begin(), tempInpo.shapeFunction[tv_master].end(), 
                        M_0e.data.begin());
                    // col-major
                    N_0e.Fill(
                        M_0e(0,0), 0.0, 0.0,
                        0.0, M_0e(0,0), 0.0, 
                        0.0, 0.0, M_0e(0,0), 
                        M_0e(0,1), 0.0, 0.0,
                        0.0, M_0e(0,1), 0.0, 
                        0.0, 0.0, M_0e(0,1), 
                        M_0e(0,2), 0.0, 0.0,
                        0.0, M_0e(0,2), 0.0, 
                        0.0, 0.0, M_0e(0,2), 
                        M_0e(0,3), 0.0, 0.0,
                        0.0, M_0e(0,3), 0.0, 
                        0.0, 0.0, M_0e(0,3));
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMM(tempPena, tempMatr_0, tempMatr_1);
                    GEMTM(normTang, tempMatr_1, tempMatr_2);
                    GEMTM(N_0e, tempMatr_2, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        auto itNici_master = tempNici_master.find(tempInpo.node[tv_master][tj]);
                        I64 row_tj = itNici_master->second;
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * row_tj + tk;
                            for(I64 tm = 0; tm < 4; ++ tm){
                                I64 node_tm = tempInpo.node[tv][tm];
                                for(I64 tn = 0; tn < 3; ++ tn){
                                    I64 free_tn = 3 * node_tm + tn;
                                    coo[cooIndex].row = free_tk;
                                    coo[cooIndex].col = free_tn;
                                    coo[cooIndex].val = elemMatr(3 * tj + tk, 3 * tm + tn);
                                    ++ cooIndex;
                                }
                            }
                        }
                    }
                }
                erroPenaDisp[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            SparseMatrix tempMatr = std::move(erroPenaDisp[ts][tv]);
            GEMM(accuContProlT[ts], tempMatr, erroPenaDisp[ts][tv], threadManager.interfaceS2M[ts], 1);
        }
    });
}

void MPLatin::Establish(){
    Log("MPLatin::Establish");
    //
    SubTransfer();
    SubAccuDomaProl();
    SubDomaAuxi();
    SubAuxiAuxi();
    SubGlobalCouple();
    SubErrorMA();
    SubErroPenaDisp();
    //
    const I64 domaSize = (*domains).size();
    const I64 inteSize = (*interfaces).size();
    inteForc.resize(inteSize);
    for(I64 ts = 0; ts < inteSize; ++ ts){
        inteForc[ts].resize(MAccu[domaSize + ts + 1] - MAccu[domaSize + ts]);
    }
    coarDisp.resize(domaSize);
    for(I64 td = 0; td < domaSize; ++ td){
        coarDisp[td].resize(MAccu[td + 1] - MAccu[td]);
    }
    globalSolution.resize(MAccu[domaSize + inteSize]);
    globalForce.resize(MAccu[domaSize + inteSize]);
}

void MPLatin::Apply(
    std::vector<AlignedVectorRx>& resuDisp, 
    const std::vector<std::array<AlignedVectorRx, 2>>& inteAuxi, 
    const std::vector<std::array<AlignedVectorRx, 2>>& inteMult){
    //
    Log("        MPLatin::Apply");
    //
    const I64 domaSize = (*domains).size();
    std::fill(std::execution::unseq, globalForce.begin(), globalForce.end(), 0.0);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        std::fill(std::execution::unseq, inteForc[ts].begin(), inteForc[ts].end(), 0.0);
        for(I64 tv = 0; tv < 2; ++ tv){
            const I64 tempDomaId = (*interfaces)[ts].domainIndex[tv];
            PEMV(erroMult[ts][tv], inteMult[ts][tv], inteForc[ts], 
                threadManager.interfaceS2M[ts], 1);
            PEMV(-1.0, erroPenaAuxi[ts][tv], inteAuxi[ts][tv], inteForc[ts], 
                threadManager.interfaceS2M[ts], 1);
            PEMV(erroPenaDisp[ts][tv], resuDisp[tempDomaId], inteForc[ts], 
                threadManager.interfaceS2M[ts], 1);
        }
        std::copy(std::execution::unseq,
            inteForc[ts].begin(), inteForc[ts].end(), 
            globalForce.begin() + MAccu[domaSize + ts]);
    });
    //
    macroDs.Solve(globalForce, globalSolution);
    //
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        const I64 MAccu_td = MAccu[td];
        const I64 M_td = MAccu[td + 1] - MAccu_td;
        std::copy(std::execution::unseq, 
            globalSolution.begin() + MAccu_td, 
            globalSolution.begin() + MAccu_td + M_td, 
            coarDisp[td].begin());
        PEMV(accuDomaProl[td], coarDisp[td], resuDisp[td], 
            threadManager.domainS2M[td], 1);
    });
}

} // namespace Ddpca