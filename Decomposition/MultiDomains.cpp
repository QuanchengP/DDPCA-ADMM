#include "MultiDomains.hpp"

namespace Ddpca {

Real MultiDomains::CalculateCharacteristicLength(){
    //
    const I64 numbDomains = domains.size();
    std::vector<Real> domainVolume(numbDomains);
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        domainVolume[td] = domains[td].mesh.Volume(
            threadManager.domainS2M[td], 1);
    });
    //
    Real totalVolume = std::transform_reduce(std::execution::unseq, 
        domainVolume.begin(), domainVolume.end(), 
        0.0, std::plus<Real>(), 
        [](const auto& val) { return val; });
    Real characteristicLength = std::pow(totalVolume / numbDomains, 0.333333333333333333);
    Log("characteristicLength = " + std::to_string(characteristicLength));
    return characteristicLength;
}

void MultiDomains::CalculatePenaltyParameter(){
    Real characteristicLength = CalculateCharacteristicLength();
    const I64 inteSize = interfaces.size();
    for(I64 ts = 0; ts < inteSize; ++ ts){
        Real averageElasticity = (domains[interfaces[ts].domainIndex[0]].elasticity 
            + domains[interfaces[ts].domainIndex[1]].elasticity) / 2.0;
        interfaces[ts].normPenaPara = NormalPenaltyCoef * averageElasticity / characteristicLength;
        interfaces[ts].tangPenaPara = TangentialPenaltyCoef * averageElasticity / characteristicLength;
    }
}

void MultiDomains::SubDomainPenalty(){
    //
    Log("MultiDomains::SubDomainPenalty");
    //
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        for(I64 tv = 0; tv < 2; ++ tv){
            SingleDomain& tempDomain = domains[tempInterface.domainIndex[tv]];
            I64 M = 3 * tempDomain.mesh.node2Coordinate.size();
            I64 N = M;
            I64 nnz = 144 * inpoSize; // 12 * 12
            Triplet* coo = AlignedAllocate<Triplet>(nnz);
            I64 cooIndex = 0;
            DenseMatrix elemMatr(12, 12), N_e(3, 12);
            std::array<Real, 4> M_e;
            if(tempInterface.frictionCoefficient == 0.0){
                DenseMatrix normVect(1, 3), tempMatr_0(1, 12), tempMatr_1(3, 12);
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
                    GEMM(normVect, N_e, tempMatr_0);
                    GEMTM(normVect, tempMatr_0, tempMatr_1);
                    GEMTM(N_e, tempMatr_1, elemMatr);
                    SCAL(tempInpo.quadratureWeight * tempInterface.normPenaPara, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
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
            }
            else{
                DenseMatrix normTang(3, 3), tempPena(3, 3), tempMatr_0(3, 12);
                DenseMatrix tempMatr_1(3, 12), tempMatr_2(3, 12);
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMM(tempPena, tempMatr_0, tempMatr_1);
                    GEMTM(normTang, tempMatr_1, tempMatr_2);
                    GEMTM(N_e, tempMatr_2, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
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
            }
            assert(cooIndex == nnz && "cooIndex != nnz");
            SparseMatrix domainPenalty;
            domainPenalty.Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
            Deallocate<Triplet>(coo, nnz);
            tempDomain.stiffness.KeepNPAdd(domainPenalty, threadManager.interfaceS2M[ts], 1);
        }
    });
}

void MultiDomains::SubNi2ci(){
    Log("MultiDomains::SubNi2ci");
    const I64 inteSize = interfaces.size();
    nodeId2ContactId.resize(inteSize);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        for(I64 tv = 0; tv < 2; ++ tv){
            std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            for(const auto& tempInpo : interfaces[ts].integralPoints){
                for(I64 tk = 0; tk < 4; tk ++){
                    const I64 tempNode = tempInpo.node[tv][tk];
                    auto itNi2ci = tempNici.find(tempNode);
                    if(itNi2ci == tempNici.end()){
                        const I64 tempSize = tempNici.size();
                        tempNici.emplace(tempNode, tempSize);
                    }
                }
            }
        }
    });
}

void MultiDomains::SubDomaInte(){
    Log("MultiDomains::SubDomaInte");
    //
    const I64 inteSize = interfaces.size();
    domaInte.resize(inteSize);
    domaPenaInte.resize(inteSize);
    domaPenaInte_T.resize(inteSize);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        for(I64 tv = 0; tv < 2; ++ tv){
            SingleDomain& tempDomain = domains[tempInterface.domainIndex[tv]];
            const I64 tempDoaminDof = 3 * tempDomain.mesh.node2Coordinate.size();
            std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = tempDoaminDof;
                I64 N = tempNici.size();
                I64 nnz = 48 * inpoSize; // 12 * 4
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* cooPena = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(12, 4), M_e(1, 4), N_e(3, 12), normVect(1, 3), tempMatr_0(3, 4);
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
                    GEMTM(normVect, M_e, tempMatr_0);
                    GEMTM(N_e, tempMatr_0, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
                            for(I64 tm = 0; tm < 4; ++ tm){
                                auto itNici = tempNici.find(tempInpo.node[tv][tm]);
                                coo[cooIndex].row = free_tk;
                                coo[cooIndex].col = itNici->second;
                                coo[cooIndex].val = elemMatr(3 * tj + tk, tm);
                                cooPena[cooIndex].row = free_tk;
                                cooPena[cooIndex].col = itNici->second;
                                cooPena[cooIndex].val = tempInterface.normPenaPara* elemMatr(3 * tj + tk, tm);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                domaInte[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                domaPenaInte[ts][tv].Coo2Csr(M, N, cooPena, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(cooPena, nnz);
            }
            else{
                I64 M = tempDoaminDof;
                I64 N = 3 * tempNici.size();
                I64 nnz = 144 * inpoSize; // 12 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* cooPena = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr_0(12, 12), elemMatr_1(12, 12);
                std::array<Real, 4> M_e;
                DenseMatrix N_e(3, 12), normTang(3, 3), tempMatr_0(3, 12), 
                    tempMatr_1(3, 12), tempPena(3, 3), tempMatr_2(3, 12);
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMTM(normTang, tempMatr_0, tempMatr_1);
                    GEMTM(N_e, tempMatr_1, elemMatr_0);
                    SCAL(tempInpo.quadratureWeight, elemMatr_0);
                    //
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(tempPena, tempMatr_0, tempMatr_1);
                    GEMTM(normTang, tempMatr_1, tempMatr_2);
                    GEMTM(N_e, tempMatr_2, elemMatr_1);
                    SCAL(tempInpo.quadratureWeight, elemMatr_1);
                    for(I64 tj = 0; tj < 4; ++ tj){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; ++ tk){
                            I64 free_tk = 3 * node_tj + tk;
                            for(I64 tm = 0; tm < 4; ++ tm){
                                auto itNici = tempNici.find(tempInpo.node[tv][tm]);
                                for(I64 tn = 0; tn < 3; ++ tn){
                                    I64 col_tn = 3 * itNici->second + tn;
                                    coo[cooIndex].row = free_tk;
                                    coo[cooIndex].col = col_tn;
                                    coo[cooIndex].val = elemMatr_0(3 * tj + tk, 3 * tm + tn);
                                    cooPena[cooIndex].row = free_tk;
                                    cooPena[cooIndex].col = col_tn;
                                    cooPena[cooIndex].val = elemMatr_1(3 * tj + tk, 3 * tm + tn);
                                    ++ cooIndex;
                                }
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                domaInte[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                domaPenaInte[ts][tv].Coo2Csr(M, N, cooPena, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(cooPena, nnz);
            }
            TRANSPOSE(domaPenaInte[ts][tv], domaPenaInte_T[ts][tv]);
        }
    });
}

void MultiDomains::SubInterfaceMass(){
    Log("MultiDomains::SubInterfaceMass");
    //
    const I64 inteSize = interfaces.size();
    inteInte.resize(inteSize);
    intePenaInte.resize(inteSize);
    inteInteDs.resize(inteSize);
    intePenaInteDs.resize(inteSize);
    //
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        for(I64 tv = 0; tv < 2; ++ tv){
            std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            const I64 niciSize = tempNici.size();
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = niciSize;
                I64 N = niciSize;
                I64 nnz = 16 * inpoSize; // 4 * 4
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* cooPena = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr(4, 4), M_e(1, 4);
                for(const IntegralPoint& tempInpo : tempInterface.integralPoints){
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.data.begin());
                    GEMTM(M_e, M_e, elemMatr);
                    SCAL(tempInpo.quadratureWeight, elemMatr);
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        I64 row_tj = iter_Cj->second;
                        for(I64 tk = 0; tk < 4; tk ++){
                            auto iter_Ck = tempNici.find(tempInpo.node[tv][tk]);
                            I64 col_tk = iter_Ck->second;
                            coo[cooIndex].row = row_tj;
                            coo[cooIndex].col = col_tk;
                            coo[cooIndex].val = elemMatr(tj, tk);
                            cooPena[cooIndex].row = row_tj;
                            cooPena[cooIndex].col = col_tk;
                            cooPena[cooIndex].val = tempInterface.normPenaPara* elemMatr(tj, tk);
                            ++ cooIndex;
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inteInte[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                intePenaInte[ts][tv].Coo2Csr(M, N, cooPena, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(cooPena, nnz);
            }
            else{
                I64 M = 3 * niciSize;
                I64 N = 3 * niciSize;
                I64 nnz = 144 * inpoSize; // 12 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                Triplet* cooPena = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr_0(12, 12), elemMatr_1(12, 12);
                std::array<Real, 4> M_e;
                DenseMatrix N_e(3, 12), normTang(3, 3), tempMatr_0(3, 12), 
                    tempMatr_1(3, 12), tempPena(3, 3), tempMatr_2(3, 12);
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMTM(normTang, tempMatr_0, tempMatr_1);
                    GEMTM(N_e, tempMatr_1, elemMatr_0);
                    SCAL(tempInpo.quadratureWeight, elemMatr_0);
                    //
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(tempPena, tempMatr_0, tempMatr_1);
                    GEMTM(normTang, tempMatr_1, tempMatr_2);
                    GEMTM(N_e, tempMatr_2, elemMatr_1);
                    SCAL(tempInpo.quadratureWeight, elemMatr_1);
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        I64 row_tj = iter_Cj->second;
                        for(I64 tk = 0; tk < 4; tk ++){
                            auto iter_Ck = tempNici.find(tempInpo.node[tv][tk]);
                            I64 col_tk = iter_Ck->second;
                            for(I64 tl = 0; tl < 3; tl ++){
                                for(I64 tm = 0; tm < 3; tm ++){
                                    coo[cooIndex].row = 3 * row_tj + tl;
                                    coo[cooIndex].col = 3 * col_tk + tm;
                                    coo[cooIndex].val = elemMatr_0(3 * tj + tl, 3 * tk + tm);
                                    cooPena[cooIndex].row = 3 * row_tj + tl;
                                    cooPena[cooIndex].col = 3 * col_tk + tm;
                                    cooPena[cooIndex].val = elemMatr_1(3 * tj + tl, 3 * tk + tm);
                                    ++ cooIndex;
                                }
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inteInte[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
                intePenaInte[ts][tv].Coo2Csr(M, N, cooPena, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(cooPena, nnz);
            }
            inteInteDs[ts][tv].Establish(inteInte[ts][tv]);
            intePenaInteDs[ts][tv].Establish(intePenaInte[ts][tv]);
        }
    });
}

void MultiDomains::SubInpoInteLagr(){
    Log("MultiDomains::SubInpoInteLagr");
    //
    const I64 inteSize = interfaces.size();
    inpoInteLagr.resize(inteSize);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        for(I64 tv = 0; tv < 2; ++ tv){
            std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            const I64 niciSize = tempNici.size();
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = inpoSize;
                I64 N = niciSize;
                I64 nnz = 4 * inpoSize;
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.begin());
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        I64 col_tj = iter_Cj->second;
                        coo[cooIndex].row = tp;
                        coo[cooIndex].col = col_tj;
                        coo[cooIndex].val = M_e[tj];
                        ++ cooIndex;
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inpoInteLagr[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            else{
                I64 M = 3 * inpoSize;
                I64 N = 3 * niciSize;
                I64 nnz = 36 * inpoSize; // 3 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr_0(3, 12), N_e(3, 12), normTang(3, 3);
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    GEMM(normTang, N_e, elemMatr_0);
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        for(I64 tk = 0; tk < 3; tk ++){
                            I64 free_tk = 3 * iter_Cj->second + tk;
                            for(I64 tm = 0; tm < 3; tm ++){
                                coo[cooIndex].row = 3 * tp + tm;
                                coo[cooIndex].col = free_tk;
                                coo[cooIndex].val = elemMatr_0(tm, 3 * tj + tk);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inpoInteLagr[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
        }
    });
}

void MultiDomains::SubInpoPenaDomaDisp(){
    Log("MultiDomains::SubInpoPenaDomaDisp");
    //
    const I64 inteSize = interfaces.size();
    inpoPenaDomaDisp.resize(inteSize);
    inpoNgapStre.resize(inteSize);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        AlignedVectorRx& tempInngStre = inpoNgapStre[ts];
        for(I64 tv = 0; tv < 2; ++ tv){
            const I64 tempDomaId = tempInterface.domainIndex[tv];
            SingleDomain& tempDomain = domains[tempDomaId];
            const I64 domaDof = 3 * tempDomain.mesh.node2Coordinate.size();
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = inpoSize;
                I64 N = domaDof;
                I64 nnz = 12 * inpoSize; // 1 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                tempInngStre.resize(inpoSize);
                DenseMatrix elemMatr(1, 12), N_e(3, 12), normVect(1, 3);
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
                    tempInngStre[tp] = tempInterface.normPenaPara * tempInpo.initialNormalGap;
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
                    GEMM(normVect, N_e, elemMatr);
                    SCAL(tempInterface.normPenaPara, elemMatr);
                    for(I64 tj = 0; tj < 4; tj ++){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; tk ++){
                            coo[cooIndex].row = tp;
                            coo[cooIndex].col = 3 * node_tj + tk;
                            coo[cooIndex].val = elemMatr(0, 3 * tj + tk);
                            ++ cooIndex;
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inpoPenaDomaDisp[ts][tv].Coo2Csr(M, N, coo, nnz,  threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            else{
                I64 M = 3 * inpoSize;
                I64 N = domaDof;
                I64 nnz = 36 * inpoSize; // 3 * 12
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                tempInngStre.assign(3 * inpoSize, 0.0);
                DenseMatrix elemMatr_0(3, 12), N_e(3, 12), normTang(3, 3), tempPena(3, 3), tempMatr_0(3, 12);
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
                    tempInngStre[3 * tp + 0] = tempInterface.normPenaPara * tempInpo.initialNormalGap;
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    tempPena.Fill(
                        tempInterface.normPenaPara, 0.0, 0.0, 
                        0.0, tempInterface.tangPenaPara, 0.0, 
                        0.0, 0.0, tempInterface.tangPenaPara);
                    GEMM(normTang, N_e, tempMatr_0);
                    GEMM(tempPena, tempMatr_0, elemMatr_0);
                    for(I64 tj = 0; tj < 4; tj ++){
                        I64 node_tj = tempInpo.node[tv][tj];
                        for(I64 tk = 0; tk < 3; tk ++){
                            for(I64 tm = 0; tm < 3; tm ++){
                                coo[cooIndex].row = 3 * tp + tm;
                                coo[cooIndex].col = 3 * node_tj + tk;
                                coo[cooIndex].val = elemMatr_0(tm, 3 * tj + tk);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inpoPenaDomaDisp[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
        }
    });
}

void MultiDomains::SubInteInpo(){
    Log("MultiDomains::SubInteInpo");
    //
    const I64 inteSize = interfaces.size();
    inteInpo.resize(inteSize);
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        for(I64 tv = 0; tv < 2; ++ tv){
            std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            const I64 niciSize = tempNici.size();
            if(tempInterface.frictionCoefficient == 0.0){
                I64 M = niciSize;
                I64 N = inpoSize;
                I64 nnz = 4 * inpoSize; // 4 * 1
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
                    std::copy(std::execution::unseq, 
                        tempInpo.shapeFunction[tv].begin(), tempInpo.shapeFunction[tv].end(), 
                        M_e.begin());
                    if(tv == 0){
                        SCAL(-1.0, M_e);
                    }
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        coo[cooIndex].row = iter_Cj->second;
                        coo[cooIndex].col = tp;
                        coo[cooIndex].val = tempInpo.quadratureWeight * M_e[tj];
                        ++ cooIndex;
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inteInpo[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
            else{
                I64 M = 3 * niciSize;
                I64 N = 3 * inpoSize;
                I64 nnz = 36 * inpoSize; // 12 * 3
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                I64 cooIndex = 0;
                DenseMatrix elemMatr_0(12, 3), N_e(3, 12), normTang(3, 3);
                std::array<Real, 4> M_e;
                for(I64 tp = 0; tp < inpoSize; ++ tp){
                    const IntegralPoint& tempInpo = tempInterface.integralPoints[tp];
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
                    // col-major
                    normTang.Fill(
                        tempInpo.basisVector[0][0], tempInpo.basisVector[1][0], tempInpo.basisVector[2][0],
                        tempInpo.basisVector[0][1], tempInpo.basisVector[1][1], tempInpo.basisVector[2][1],
                        tempInpo.basisVector[0][2], tempInpo.basisVector[1][2], tempInpo.basisVector[2][2]);
                    GEMTMT(N_e, normTang, elemMatr_0);
                    SCAL(tempInpo.quadratureWeight, elemMatr_0);
                    if(tv == 0){
                        SCAL(-1.0, elemMatr_0);
                    }
                    for(I64 tj = 0; tj < 4; tj ++){
                        auto iter_Cj = tempNici.find(tempInpo.node[tv][tj]);
                        for(I64 tk = 0; tk < 3; tk ++){
                            I64 row_tk = 3 * iter_Cj->second + tk;
                            for(I64 tm = 0; tm < 3; tm ++){
                                coo[cooIndex].row = row_tk;
                                coo[cooIndex].col = 3 * tp + tm;
                                coo[cooIndex].val = elemMatr_0(3 * tj + tk, tm);
                                ++ cooIndex;
                            }
                        }
                    }
                }
                assert(cooIndex == nnz && "cooIndex != nnz");
                inteInpo[ts][tv].Coo2Csr(M, N, coo, nnz, threadManager.interfaceS2M[ts], 1);
                Deallocate<Triplet>(coo, nnz);
            }
        }
    });
}

void MultiDomains::Establish(){
    //
    Log("MultiDomains::Establish");
    //
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        SingleDomain& tempDomain = domains[td];
        tempDomain.conjGrad.preconditioner = PreconditionerMultigrid;
        tempDomain.conjGrad.pcGM.relaxation = RelaxationSSOR;
        //
        tempDomain.geomMult.dataMesh = &(tempDomain.mesh);
        tempDomain.geomMult.dataBoundary = &(tempDomain.boundary);
        tempDomain.geomMult.dataStiffness = &(tempDomain.stiffness);
        tempDomain.geomMult.dataPcg = &(tempDomain.conjGrad.pcGM);
        tempDomain.geomMult.Transfer(threadManager.domainS2M[td], 1);
        //must after Transfer
        FemStiffnessMatrix(tempDomain.mesh, 
            tempDomain.elasticity, 
            tempDomain.poissonRatio, 
            tempDomain.stiffness, 
            threadManager.domainS2M[td], 
            1);
    });
    //
    SubDomainPenalty();
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        SingleDomain& tempDomain = domains[td];
        tempDomain.geomMult.Establish(threadManager.domainS2M[td], 1);
        const I64 maxLevel = tempDomain.conjGrad.pcGM.realProlong.size();
        const SparseMatrix& tempStiff = tempDomain.conjGrad.pcGM.hierarchyStiffness[maxLevel];
        if(tempStiff.M > SingleDomain::dscoLimi){
            tempDomain.conjGrad.Establish(tempStiff, threadManager.domainS2M[td], 1);
        }
        else{
            tempDomain.cholDs.Establish(tempStiff);
            //can not be omitted, the MPLatin will use the following matrices
            tempDomain.conjGrad.pcGM.realProlongT.resize(maxLevel);
            //
            I64 numbPart = threadManager.domainS2M[td].size();
            std::vector<I64> startIndex(numbPart), endIndex(numbPart);
            EvenlyDistribute(maxLevel, numbPart, startIndex, endIndex);
            threadManager.RunTask(1, threadManager.domainS2M[td], 
                [&](I64 tp){
                    I64 start_tp = startIndex[tp];
                    I64 end_tp = endIndex[tp];
                    for(I64 tl = start_tp; tl < end_tp; ++ tl){
                        TRANSPOSE(tempDomain.conjGrad.pcGM.realProlong[tl], 
                            tempDomain.conjGrad.pcGM.realProlongT[tl]);
                    }
                });
            //
            for(I64 tl = maxLevel - 1; tl >= 0; -- tl){
                SparseMatrix tempIK;
                GEMM(tempDomain.conjGrad.pcGM.realProlongT[tl], 
                    tempDomain.conjGrad.pcGM.hierarchyStiffness[tl + 1], 
                    tempIK, 
                    threadManager.domainS2M[td], 1);
                GEMM(tempIK, 
                    tempDomain.conjGrad.pcGM.realProlong[tl], 
                    tempDomain.conjGrad.pcGM.hierarchyStiffness[tl], 
                    threadManager.domainS2M[td], 1);
            }
        }
    });
    //
    SubNi2ci();
    SubDomaInte();
    SubInterfaceMass();
    SubInpoInteLagr();
    SubInpoPenaDomaDisp();
    SubInteInpo();
    //
    mpLatin.domains = &domains;
    mpLatin.interfaces = &interfaces;
    mpLatin.nodeId2ContactId = &nodeId2ContactId;
    mpLatin.Establish();
}

void MultiDomains::OutpAuxiMult(
    const std::vector<AlignedVectorRx>& inpoGamm, 
    const std::string& directoryPath){
    //
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const ContactInterface& tempInterface = interfaces[ts];
        const I64 inpoSize = tempInterface.integralPoints.size();
        std::ofstream tempOfst(directoryPath + "/resuIpga_" + std::to_string(ts) + ".txt", std::ios::out);
        tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
        for(I64 tp = 0; tp < inpoSize; ++ tp){
            const IntegralPoint& tempInpo = interfaces[ts].integralPoints[tp];
            tempOfst << std::setw(30) << tempInpo.contactPoint[0][0] 
                << std::setw(30) << tempInpo.contactPoint[0][1] 
                << std::setw(30) << tempInpo.contactPoint[0][2];
            if(tempInterface.frictionCoefficient == 0.0){
                tempOfst << std::setw(30) << inpoGamm[ts][tp];
            }
            else{
                tempOfst << std::setw(30) << inpoGamm[ts][3 * tp + 0] 
                    << std::setw(30) << inpoGamm[ts][3 * tp + 1] 
                    << std::setw(30) << inpoGamm[ts][3 * tp + 2];
            }
            tempOfst << "\n";
        }
        tempOfst.close();
        for(I64 tv = 0; tv < 2; ++ tv){
            tempOfst.open(directoryPath + "/resuAumu_" + std::to_string(ts) 
                + "_" + std::to_string(tv) + ".txt", std::ios::out);
            tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
            const SingleDomain& tempDomain = domains[tempInterface.domainIndex[tv]];
            const std::map<I64, I64>& tempNici = nodeId2ContactId[ts][tv];
            for(auto iter_Cj = tempNici.begin(); iter_Cj != tempNici.end(); ++ iter_Cj){
                I64 tempNode = iter_Cj->first;
                auto it_Noco = tempDomain.mesh.node2Coordinate.find(tempNode);
                tempOfst << std::setw(30) << it_Noco->second[0] 
                    << std::setw(30) << it_Noco->second[1] 
                    << std::setw(30) << it_Noco->second[2];
                I64 tempCid = iter_Cj->second;
                if(tempInterface.frictionCoefficient == 0.0){
                    tempOfst << std::setw(30) << inteAuxi[ts][tv][tempCid]
                        << std::setw(30) << inteMult[ts][tv][tempCid];
                }
                else{
                    tempOfst << std::setw(30) << inteAuxi[ts][tv][3 * tempCid + 0]
                        << std::setw(30) << inteAuxi[ts][tv][3 * tempCid + 1]
                        << std::setw(30) << inteAuxi[ts][tv][3 * tempCid + 2]
                        << std::setw(30) << inteMult[ts][tv][3 * tempCid + 0]
                        << std::setw(30) << inteMult[ts][tv][3 * tempCid + 1]
                        << std::setw(30) << inteMult[ts][tv][3 * tempCid + 2];
                }
                tempOfst << "\n";
            }
            tempOfst.close();
        }
    });
}

void MultiDomains::ADMM(const std::string& directoryPath){
    Log("MultiDomains::ADMM");
    //
	std::ofstream tempOfst(directoryPath + "/resuMoni.txt", std::ios::out);
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
    //
    const I64 domaSize = domains.size();
    const I64 inteSize = interfaces.size();
    resuDisp.resize(domaSize);
    inteAuxi.resize(inteSize);
    inteMult.resize(inteSize);
    std::vector<AlignedVectorRx> resuDisp_0(domaSize);
    std::vector<std::array<AlignedVectorRx, 2>> inteAuxi_0(inteSize);
    std::vector<std::array<AlignedVectorRx, 2>> inteMult_0(inteSize); // currently not used
    //
    std::vector<AlignedVectorRx> additionalForce(domaSize);
    std::vector<AlignedVectorRx> solverForce(domaSize);
    std::vector<AlignedVectorRx> solverDisp(domaSize);
    std::vector<AlignedVectorRx> inpoGamm(inteSize);
    std::vector<std::array<AlignedVectorRx, 2>> inteForc(inteSize);
    std::vector<std::array<AlignedVectorRx, 2>> inteDeltMult(inteSize);
    std::vector<std::vector<I64>> startIndex(inteSize), endIndex(inteSize);
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        const I64 tempSize = 3 * domains[td].mesh.coordinate2Node.size();
        resuDisp[td].assign(tempSize, 0.0);
        //
        resuDisp_0[td].assign(tempSize, 0.0);
        additionalForce[td].resize(tempSize);
        solverForce[td].resize(domains[td].geomMult.romT.M);
        solverDisp[td].resize(solverForce[td].size());
    });
    threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
        const I64 inpoSize = interfaces[ts].integralPoints.size();
        if(interfaces[ts].frictionCoefficient == 0.0){
            inpoGamm[ts].resize(inpoSize);
        }
        else{
            inpoGamm[ts].resize(3 * inpoSize);
        }
        for(I64 tv = 0; tv < 2; ++ tv){
            if(interfaces[ts].frictionCoefficient == 0.0){
                const I64 tempSize = nodeId2ContactId[ts][tv].size();
                inteAuxi[ts][tv].assign(tempSize, 0.0);
                inteMult[ts][tv].assign(tempSize, 0.0);
                //
                inteAuxi_0[ts][tv].assign(tempSize, 0.0);
                inteMult_0[ts][tv].assign(tempSize, 0.0);
                inteForc[ts][tv].assign(tempSize, 0.0);
                inteDeltMult[ts][tv].assign(tempSize, 0.0);
            }
            else{
                const I64 tempSize = 3 * nodeId2ContactId[ts][tv].size();
                inteAuxi[ts][tv].assign(tempSize, 0.0);
                inteMult[ts][tv].assign(tempSize, 0.0);
                //
                inteAuxi_0[ts][tv].assign(tempSize, 0.0);
                inteMult_0[ts][tv].assign(tempSize, 0.0);
                inteForc[ts][tv].assign(tempSize, 0.0);
                inteDeltMult[ts][tv].assign(tempSize, 0.0);
            }
        }
        I64 numbPart = threadManager.interfaceS2M[ts].size();
        startIndex[ts].resize(numbPart);
        endIndex[ts].resize(numbPart);
        EvenlyDistribute(inpoSize, numbPart, startIndex[ts], endIndex[ts]);
    });
    //
    Real stoppingCriterion = 1.0E-6; // norm but not the square of the norm
    I64 tc = 0;
    const I64 maxIteration = 3000;
	for(tc = 0; tc < maxIteration; ++ tc){
        //
		Log( "        MultiDomains::The " + std::to_string(tc) + "-th iteration...");
        threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
            std::copy(std::execution::unseq, 
                resuDisp[td].begin(), resuDisp[td].end(), resuDisp_0[td].begin());
        });
        threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
            for(I64 tv = 0; tv < 2; ++ tv){
                std::copy(std::execution::unseq, 
                    inteAuxi[ts][tv].begin(), inteAuxi[ts][tv].end(), inteAuxi_0[ts][tv].begin());
                std::copy(std::execution::unseq, 
                    inteMult[ts][tv].begin(), inteMult[ts][tv].end(), inteMult_0[ts][tv].begin());
            }
        });
		//***********************************subdomain equilibrium*********************************
        Log("        MultiDomains::Subdomain equilibrium");
        // for(I64 td = 0; td < domaSize; ++ td){
        threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
            std::fill(additionalForce[td].begin(), additionalForce[td].end(), 0.0);
            for(I64 ts = 0; ts < inteSize; ++ ts){
                for(I64 tv = 0; tv < 2; ++ tv){
                    if(interfaces[ts].domainIndex[tv] != td){
                        continue;
                    }
                    PEMV(domaPenaInte[ts][tv], inteAuxi[ts][tv], 
                        additionalForce[td], threadManager.domainS2M[td], 1);
                    PEMV(-1.0, domaInte[ts][tv], inteMult[ts][tv], 
                        additionalForce[td], threadManager.domainS2M[td], 1);
                }
            }
            MV(1.0, domains[td].geomMult.romT, additionalForce[td], 
                1.0, domains[td].boundary.loadVect, 
                solverForce[td], threadManager.domainS2M[td], 1);
            const I64 maxLevel = domains[td].conjGrad.pcGM.realProlong.size();
            if(solverDisp[td].size() > SingleDomain::dscoLimi){
                std::fill(solverDisp[td].begin(), solverDisp[td].end(), 0.0);
                domains[td].conjGrad.Solve(
                    domains[td].conjGrad.pcGM.hierarchyStiffness[maxLevel], 
                    solverForce[td], solverDisp[td], threadManager.domainS2M[td], 1);
            }
            else{
                domains[td].cholDs.Solve(solverForce[td], solverDisp[td]);
            }
            domains[td].geomMult.DispPost(solverDisp[td], resuDisp[td], 
                threadManager.domainS2M[td], 1);
            #ifndef NDEBUG
                domains[td].mesh.OutputDisplacement(directoryPath, td, resuDisp[td]);
            #endif
        });
        // }
		//********************************Macroscopic problem in LATIN*****************************
        mpLatin.Apply(resuDisp, inteAuxi, inteMult);
        #ifndef NDEBUG
            threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
                (domains[td]).mesh.OutputDisplacement(directoryPath, td, resuDisp[td]);
            });
        #endif
		//***********************************interface equilibrium*********************************
        Log("        MultiDomains::Interface equilibrium");
        threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
            AlignedVectorRx& tempInpoGamm = inpoGamm[ts];
            std::fill(tempInpoGamm.begin(), tempInpoGamm.end(), 0.0);
            PEMV(inpoInteLagr[ts][0], inteMult[ts][0], tempInpoGamm, threadManager.interfaceS2M[ts], 1);
            PEMV(-1.0, inpoInteLagr[ts][1], inteMult[ts][1], tempInpoGamm, threadManager.interfaceS2M[ts], 1);
            PEMV(inpoPenaDomaDisp[ts][0], resuDisp[interfaces[ts].domainIndex[0]], 
                tempInpoGamm, threadManager.interfaceS2M[ts], 1);
            PEMV(-1.0, inpoPenaDomaDisp[ts][1], resuDisp[interfaces[ts].domainIndex[1]], 
                tempInpoGamm, threadManager.interfaceS2M[ts], 1);
            AXPY(-1.0, inpoNgapStre[ts], tempInpoGamm);
            SCAL(0.5, tempInpoGamm);
            //
            bool tempFlagZero = (interfaces[ts].frictionCoefficient == 0.0);
            bool tempFlagGreat = (interfaces[ts].frictionCoefficient > 0.0);
            if(tempFlagZero || tempFlagGreat){
                threadManager.RunTask(1, threadManager.interfaceS2M[ts], [&](I64 tp){
                    I64 start_tp = startIndex[ts][tp];
                    I64 end_tp = endIndex[ts][tp];
                    for(I64 ti = start_tp; ti < end_tp; ++ ti){
                        if(tempFlagZero){
                            tempInpoGamm[ti] = std::max(0.0, tempInpoGamm[ti]);
                        }
                        else{
                            tempInpoGamm[ti * 3 + 0] = std::max(0.0, tempInpoGamm[ti * 3 + 0]);
                        }
                    }
                });
            }
            if(tempFlagGreat){
                threadManager.RunTask(1, threadManager.interfaceS2M[ts], [&](I64 tp){
                    I64 start_tp = startIndex[ts][tp];
                    I64 end_tp = endIndex[ts][tp];
                    for(I64 ti = start_tp; ti < end_tp; ++ ti){
                        if(tempInpoGamm[ti * 3 + 0] > 0.0){
                            Real tempSlide = interfaces[ts].frictionCoefficient * tempInpoGamm[ti * 3 + 0];
                            std::array<Real, 2> tempStick = {tempInpoGamm[ti * 3 + 1], tempInpoGamm[ti * 3 + 2]};
                            Real tempNorm = NRM2(tempStick);
                            if(tempNorm >= tempSlide){
                                SCAL(tempSlide / tempNorm, tempStick);
                                tempInpoGamm[ti * 3 + 1] = tempStick[0];
                                tempInpoGamm[ti * 3 + 2] = tempStick[1];
                            }
                        }
                        else{
                            tempInpoGamm[ti * 3 + 1] = 0.0;
                            tempInpoGamm[ti * 3 + 2] = 0.0;
                        }
                    }
                });
            }
            //
            for(I64 tv = 0; tv < 2; ++ tv){
                AlignedVectorRx& tempInteForc = inteForc[ts][tv];
                MV(domaPenaInte_T[ts][tv], resuDisp[interfaces[ts].domainIndex[tv]], 
                    tempInteForc, threadManager.interfaceS2M[ts], 1);
                PEMV(inteInte[ts][tv], inteMult[ts][tv], 
                    tempInteForc, threadManager.interfaceS2M[ts], 1);
                PEMV(inteInpo[ts][tv], tempInpoGamm, 
                    tempInteForc, threadManager.interfaceS2M[ts], 1);
                intePenaInteDs[ts][tv].Solve(tempInteForc, inteAuxi[ts][tv]);
            }
        });
		//*************************************update multiplier***********************************
        Log("        MultiDomains::Update multiplier");
        threadManager.RunTask(0, threadManager.interfaceS2S, [&](I64 ts){
            const ContactInterface& tempInterface = interfaces[ts];
            for(I64 tv = 0; tv < 2; ++ tv){
                AlignedVectorRx& tempInteForc = inteForc[ts][tv];
                AlignedVectorRx& tempForc = inteDeltMult[ts][tv];
                MV(domaPenaInte_T[ts][tv], resuDisp[tempInterface.domainIndex[tv]], 
                    tempInteForc, threadManager.interfaceS2M[ts], 1);
                PEMV(-1.0, intePenaInte[ts][tv], inteAuxi[ts][tv], tempInteForc, threadManager.interfaceS2M[ts], 1);
                inteInteDs[ts][tv].Solve(tempInteForc, tempForc);
                XPEY(inteMult[ts][tv], tempForc);
            }
        });
        #ifndef NDEBUG
            OutpAuxiMult(inpoGamm, directoryPath);
        #endif
		//************************************stopping criterion***********************************
        bool tempFlag = true;
        for(I64 td = 0; td < domaSize; ++ td){
            Real tempNorm = NRM2(resuDisp[td]);
            AXPY(-1.0, resuDisp[td], resuDisp_0[td]);
            Real diffNorm = NRM2(resuDisp_0[td]);
            tempOfst << std::setw(30) << diffNorm << std::setw(30) << tempNorm;
            if(diffNorm > stoppingCriterion * tempNorm){
                tempFlag = false;
            }
        }
        tempOfst << "\n";
        for(I64 ts = 0; ts < inteSize; ++ ts){
            for(I64 tv = 0; tv < 2; ++ tv){
                Real tempNorm = NRM2(inteAuxi[ts][tv]);
                AXPY(-1.0, inteAuxi[ts][tv], inteAuxi_0[ts][tv]);
                Real diffNorm = NRM2(inteAuxi_0[ts][tv]);
                if(diffNorm > stoppingCriterion * tempNorm){
                    tempFlag = false;
                    break;
                }
            }
        }
        if(tempFlag){
            break;
        }
    }
	tempOfst.close();
    if(tc >= maxIteration){
        Log("        MultiDomains::Nonconvergence in MultiDomains::ADMM");
    }
    else{
        Log("        MultiDomains::Converge after " + std::to_string(tc) + "-th iteration");
    }
    threadManager.RunTask(0, threadManager.domainS2S, [&](I64 td){
        (domains[td]).mesh.OutputDisplacement(directoryPath, td, resuDisp[td]);
    });
    OutpAuxiMult(inpoGamm, directoryPath);
}

} // namespace Ddpca