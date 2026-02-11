#include "GeometricMultigrid.hpp"
#include "../General/AlignedVector.hpp"

#include <numeric>
#include <execution>

namespace Ddpca{

void GeometricMultigrid::Preprocess(
    std::vector<std::map<std::vector<I64>,I64>>& ininTran, 
    std::vector<std::set<I64>>& levelNode_s, 
    bool willMultigrid){
	//
    Log("GeometricMultigrid::Preprocess");
    //the level difference of two adjcent elements <= 1
	//no need to successively from level 0 to level dataMesh->maxiLeve - 1?????
	//
    //the edge center of nodes 0/1 is the node 1 of children element 0
	static const std::vector<std::vector<std::vector<I64>>> elementLine = {{
		{0, 1, 0, 1}, {1, 2, 1, 2}, {2, 3, 3, 3}, {3, 0, 2, 0}, 
		{0, 4, 0, 4}, {1, 5, 1, 5}, {2, 6, 3, 6}, {3, 7, 2, 7}, 
		{4, 5, 4, 5}, {5, 6, 5, 6}, {6, 7, 7, 7}, {7, 4, 6, 4}
	}, {
		{0, 1, 0, 1}, {1, 2, 1, 2}, {2, 3, 3, 3}, {3, 0, 2, 0}, 
		{4, 5, 0, 5}, {5, 6, 1, 6}, {6, 7, 3, 7}, {7, 4, 2, 4}
	}, {
		{0, 3, 0, 3}, {3, 7, 1, 7}, {7, 4, 3, 4}, {4, 0, 2, 0}, 
		{1, 2, 0, 2}, {2, 6, 1, 6}, {6, 5, 3, 5}, {5, 1, 2, 1}
	}, {
		{0, 4, 0, 4}, {4, 5, 1, 5}, {5, 1, 3, 1}, {1, 0, 2, 0}, 
		{3, 7, 0, 7}, {7, 6, 1, 6}, {6, 2, 3, 2}, {2, 3, 2, 3}
	}, {
		{0, 1, 0, 1}, {2, 3, 0, 2}, {4, 5, 0, 5}, {6, 7, 0, 6}
	}, {
		{0, 3, 0, 3}, {1, 2, 0, 2}, {4, 7, 0, 7}, {5, 6, 0, 6}
	}, {
		{0, 4, 0, 4}, {1, 5, 0, 5}, {3, 7, 0, 7}, {2, 6, 0, 6}
	}};
    //the face center of nodes 0/1/2/3 is the node 2 of children element 0
	static const std::vector<std::vector<std::vector<I64>>> elementFace = {{
		{0, 1, 2, 3, 0, 2}, {4, 5, 6, 7, 4, 6}, 
		{0, 3, 7, 4, 0, 7}, {1, 2, 6, 5, 3, 5}, 
		{0, 4, 5, 1, 0, 5}, {3, 7, 6, 2, 3, 7}
	}, {
		{0, 1, 2, 3, 0, 2}, {4, 5, 6, 7, 0, 6}
	}, {
		{0, 3, 7, 4, 0, 7}, {1, 2, 6, 5, 0, 6}
	}, {
		{0, 4, 5, 1, 0, 5}, {3, 7, 6, 2, 0, 6}
	}, {
	}, {
	}, {
	}};
	const I64 maxLevel = dataMesh->maxLevel;
	const std::vector<OctreeElement>& elements = dataMesh->elements;
	const std::map<std::array<I64, 2>, std::set<I64>>& lineUsedByElement 
        = dataMesh->lineUsedByElement;
	const std::map<std::array<I64, 4>, std::set<I64>>& faceUsedByElement 
        = dataMesh->faceUsedByElement;
	//
	I64 elementSize = elements.size();
	for(I64 ti = 0; ti < elementSize; ++ ti){
        I64 tempLevel = elements[ti].level;
		if(tempLevel == 0 && willMultigrid == true){
			for(I64 tj = 0; tj < 8; ++ tj){
				levelNode_s[0].emplace(elements[ti].cornerNodes[tj]);
			}
		}
		if(elements[ti].children.size() != 0){
			int pattern = elements[ti].refinementPattern;
			if(pattern == 0 && willMultigrid == true){
				//
				std::vector<I64> cornerNodes(
                    elements[ti].cornerNodes.begin(), elements[ti].cornerNodes.end());
				std::sort(cornerNodes.begin(), cornerNodes.end());
                I64 sonNode = elements[elements[ti].children[0]].cornerNodes[6];
				ininTran[tempLevel].emplace(cornerNodes, sonNode);
				levelNode_s[tempLevel + 1].emplace(sonNode);
			}
			//
			I64 numbElps = elementLine[pattern].size();
			for(I64 tj = 0; tj < numbElps; ++ tj){
				std::array<I64, 2> tempLine = {
                    elements[ti].cornerNodes[elementLine[pattern][tj][0]], 
					elements[ti].cornerNodes[elementLine[pattern][tj][1]]};
				std::sort(tempLine.begin(), tempLine.end());
				auto iteratorLube = lineUsedByElement.find(tempLine);
				//iteratorLube != lineUsedByElement.end();
				bool tempFlag = false;
				for(const auto& byElement : (iteratorLube->second)){
					if(elements[byElement].children.size() == 0){
						tempFlag = true;
						break;
					}
				}
				I64 tempNode = elements[
						elements[ti].children[elementLine[pattern][tj][2]]
					].cornerNodes[elementLine[pattern][tj][3]];
				if(tempFlag == true){
					ininTran[maxLevel].emplace(
						std::vector<I64>(tempLine.begin(), tempLine.end()), tempNode);
					levelNode_s[maxLevel + 1].emplace(tempNode);
				}
				else if(willMultigrid == true){
					ininTran[tempLevel].emplace(
						std::vector<I64>(tempLine.begin(), tempLine.end()), tempNode);
					levelNode_s[tempLevel + 1].emplace(tempNode);
				}
			}
			//
			I64 numbEfps = elementFace[pattern].size();
			for(I64 tj = 0; tj < numbEfps; ++ tj){
				std::array<I64, 4> tempFace = {
                    elements[ti].cornerNodes[elementFace[pattern][tj][0]], 
					elements[ti].cornerNodes[elementFace[pattern][tj][1]], 
					elements[ti].cornerNodes[elementFace[pattern][tj][2]], 
					elements[ti].cornerNodes[elementFace[pattern][tj][3]]};
				std::sort(tempFace.begin(), tempFace.end());
				auto iteratorFube = faceUsedByElement.find(tempFace);
				//iteratorFube != faceUsedByElement.end();
				bool tempFlag = false;
				for(const auto& byElement : (iteratorFube->second)){
					if(elements[byElement].children.size() == 0){
						tempFlag = true;
						break;
					}
				}
				I64 tempNode = elements[
						elements[ti].children[elementFace[pattern][tj][4]]
					].cornerNodes[elementFace[pattern][tj][5]];
				if(tempFlag == true){
					ininTran[maxLevel].emplace(
						std::vector<I64>(tempFace.begin(), tempFace.end()), tempNode);
					levelNode_s[maxLevel + 1].emplace(tempNode);
				}
				else if(willMultigrid == true){
					ininTran[tempLevel].emplace(
						std::vector<I64>(tempFace.begin(), tempFace.end()), tempNode);
					levelNode_s[tempLevel + 1].emplace(tempNode);
				}
			}
		}
	}
}

void GeometricMultigrid::NestedHangingNode(
    const std::vector<std::map<std::vector<I64>,I64>>& ininTran){
    //
    Log("GeometricMultigrid::NestedHangingNode");
	std::map<I64,Coordinate>& node2Coordinate = dataMesh->node2Coordinate;
	std::map<Coordinate,I64>& coordinate2Node = dataMesh->coordinate2Node;
	//
    I64 maxLevel = ininTran.size() - 1;
	for(const auto& fathersSon : ininTran[maxLevel]){
		I64 outputNode = fathersSon.second;
		Coordinate outputCoordinate(0.0, 0.0, 0.0);
		I64 numbFsfs = (fathersSon.first).size();
		for(I64 ti = 0; ti < numbFsfs; ++ ti){
			I64 inputNode = (fathersSon.first)[ti];
			const Coordinate& inputCoordinate = node2Coordinate[inputNode];
            XPEY(outputCoordinate.data, inputCoordinate.data);
		}
        SCAL(1.0 / numbFsfs, outputCoordinate.data);
		//
		auto iteratorNcfo = node2Coordinate.find(outputNode);
        assert(iteratorNcfo != node2Coordinate.end()
            && "Can not find outputNode in node2Coordinate");
        //Extract the node (without releasing memory)
        I64 nodeIndex = coordinate2Node[iteratorNcfo->second];
        coordinate2Node.erase(iteratorNcfo->second);
        iteratorNcfo->second = outputCoordinate;
        coordinate2Node[outputCoordinate] = nodeIndex;
		//
		hanging2Coarser.emplace(fathersSon.second, fathersSon.first);
	}
	coarser2Hanging = ininTran[maxLevel];
}

void GeometricMultigrid::Transfer(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("GeometricMultigrid::Transfer");
    //
    const I64 maxLevel = dataMesh->maxLevel;
    const I64 numbDmnc = (dataMesh->node2Coordinate).size();
    //
    std::vector<std::map<std::vector<I64>,I64>> ininTran;
	std::vector<std::set<I64>> levelNode_s;
	ininTran.resize(maxLevel + 1);
	levelNode_s.resize(maxLevel + 2);
    Preprocess(ininTran, levelNode_s, true);
    NestedHangingNode(ininTran);
	//
	levelNode.assign(maxLevel + 2, std::vector<I64>());
    nodeLevelPosition.assign(numbDmnc, {0, 0});
	positionNode.assign(numbDmnc, 0);
    //
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel + 2, numbPart, startIndex, endIndex);
    //
    std::function<void(I64)> taskFunction_0 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                levelNode[ti].reserve(levelNode_s[ti].size());
                std::copy(std::execution::unseq, 
                    levelNode_s[ti].begin(), 
                    levelNode_s[ti].end(), 
                    std::back_inserter(levelNode[ti]));
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
    //
    std::vector<I64> leveAccu(maxLevel + 3, 0);
    for(I64 ti = 1; ti <= maxLevel + 2; ++ ti){
        leveAccu[ti] = leveAccu[ti - 1] + levelNode[ti - 1].size();
    }
    std::function<void(I64)> taskFunction_1 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 ti_maxi = levelNode[ti].size();
                for(I64 tj = 0; tj < ti_maxi; ++ tj){
                    I64 node_ij = levelNode[ti][tj];
                    I64 acnu_ij = leveAccu[ti] + tj;
                    nodeLevelPosition[node_ij][0] = ti;
                    nodeLevelPosition[node_ij][1] = acnu_ij;
                    positionNode[acnu_ij] = node_ij;
                }
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
    //
    scalarProlong.resize(maxLevel + 1);
    EvenlyDistribute(maxLevel + 1, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction_2 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                //
                I64 M = leveAccu[ti + 2];
                I64 N = leveAccu[ti + 1];
                I64 ti_max = leveAccu[ti + 1];
                I64 nnz = ti_max;
                for(const auto& fathersSon : ininTran[ti]){
                    nnz += (fathersSon.first).size();
                }
                Triplet* coo = AlignedAllocate<Triplet>(nnz);
                //
                for(I64 tj = 0; tj < ti_max; ++ tj){
                    coo[tj].row = tj;
                    coo[tj].col = tj;
                    coo[tj].val = 1.0;
                }
                I64 tempIndex = ti_max;
                for(const auto& fathersSon : ininTran[ti]){
                    I64 sonLepo = nodeLevelPosition[fathersSon.second][1];
                    I64 numbFsfs = (fathersSon.first).size();
                    Real invNumb = 1.0 / numbFsfs;
                    for(I64 tj = 0; tj < numbFsfs; ++ tj){
                        I64 fatherLepo = nodeLevelPosition[(fathersSon.first)[tj]][1];
                        coo[tempIndex].row = sonLepo;
                        coo[tempIndex].col = fatherLepo;
                        coo[tempIndex].val = invNumb;
                        ++ tempIndex;
                    }
                }
                //
                scalarProlong[ti].Coo2Csr(M, N, coo, nnz, threadManager.levelS2M[ti], nestLevel + 1);
                Deallocate<Triplet>(coo, nnz);
            }
        };
    switch(nestLevel){
        case 0: case 1: 
            threadManager.RunTask(nestLevel, threadTask, taskFunction_2);
            break;
        default:
            taskFunction_2(0);
            break;
    }
}

void GeometricMultigrid::DofRotateOrder(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("GeometricMultigrid::DofRotateOrder");
    //
    const I64 numbDmnc = (dataMesh->node2Coordinate).size();
    const std::map<I64, DenseMatrix>& nodeRotation = dataBoundary->nodeRotation;
    //
    I64 M = dataStiffness->M;
    I64 N = M;
    I64 nnz = M + 6 * nodeRotation.size();
    Triplet* coo = AlignedAllocate<Triplet>(nnz);
    // std::vector<Triplet> coo(nnz);
    //
    I64 tempIndex = 0;
    for(I64 ti = 0; ti < numbDmnc; ++ ti){
        auto iteratorNr = nodeRotation.find(ti);
        if(iteratorNr == nodeRotation.end()){
            for(I64 tj = 0; tj < 3; ++ tj){
                coo[tempIndex].row = 3 * ti + tj;
                coo[tempIndex].col = 3 * ti + tj;
                coo[tempIndex].val = 1.0;
                ++ tempIndex;
            }
        }
        else{
            for(I64 tj = 0; tj < 3; ++ tj){
                for(I64 tk = 0; tk < 3; ++ tk){
                    coo[tempIndex].row = 3 * ti + tj;
                    coo[tempIndex].col = 3 * ti + tk;
                    coo[tempIndex].val = (iteratorNr->second)(tj, tk);
                    ++ tempIndex;
                }
            }
        }
    }
    xyzRtz.Coo2Csr(M, N, coo, nnz, threadTask, nestLevel);
    Deallocate<Triplet>(coo, nnz);
    //
    M = dataStiffness->M;
    N = M;
    nnz = M;
    coo = AlignedAllocate<Triplet>(nnz);
    // coo.resize(nnz);
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(numbDmnc, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                I64 resultPosiiton = 3 * nodeLevelPosition[ti][1];
                for(I64 tj = 0; tj < 3; ++ tj){
                    coo[3 * ti + tj].row = 3 * ti + tj;
                    coo[3 * ti + tj].col = resultPosiiton + tj;
                    coo[3 * ti + tj].val = 1.0;
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
    //
    userSolver.Coo2Csr(M, N, coo, nnz, threadTask, nestLevel);
    Deallocate<Triplet>(coo, nnz);

    //
    GEMM(xyzRtz, userSolver, rotateOrder, threadTask, nestLevel);
}

void GeometricMultigrid::RotateProlong(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("GeometricMultigrid::RotateProlong");
	const I64& maxLevel = dataMesh->maxLevel;
    const std::map<I64, DenseMatrix>& nodeRotation = dataBoundary->nodeRotation;
    //
    rotateProlong.resize(maxLevel + 1);
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel + 1, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tl = start_tp; tl < end_tp; ++ tl){
                std::vector<Triplet> coo;
                I64 numbRserve = 3.3 * scalarProlong[tl].nnz;
                coo.reserve(numbRserve);
                I64 ti_max = 3 * scalarProlong[tl].N;
                for(I64 ti = 0; ti < ti_max; ++ ti){
                    coo.emplace_back(ti, ti, 1.0);
                }
                for(I64 ti = scalarProlong[tl].N; ti < scalarProlong[tl].M; ++ ti){
                    I64 tj_start = scalarProlong[tl].row_ptr[ti];
                    I64 tj_end = scalarProlong[tl].row_ptr[ti + 1];
                    for(I64 tj = tj_start; tj < tj_end; ++ tj){
                        I64 col_tj = scalarProlong[tl].col_ind[tj];
                        Real val_tj = scalarProlong[tl].val[tj];
                        //
                        I64 sonNode = positionNode[ti];
                        auto iteratorSon = nodeRotation.find(sonNode);
                        I64 fatherNode = positionNode[col_tj];
                        auto iteratorFather = nodeRotation.find(fatherNode);
                        if(iteratorSon == nodeRotation.end() && iteratorFather == nodeRotation.end()){
                            for(I64 tk = 0; tk < 3; ++ tk){
                                coo.emplace_back(3 * ti + tk, 3 * col_tj + tk, val_tj);
                            }
                            continue;
                        }
                        //
                        DenseMatrix tempRotate(3,3);
                        tempRotate(0,0) = val_tj;
                        tempRotate(1,1) = val_tj;
                        tempRotate(2,2) = val_tj;
                        //
                        DenseMatrix tempResultS(3,3);
                        if(iteratorSon != nodeRotation.end()){
                            GEMTM((iteratorSon->second), tempRotate, tempResultS);
                        }
                        else{
                            tempResultS = tempRotate;
                        }
                        //
                        DenseMatrix tempResultF(3,3);
                        if(iteratorFather != nodeRotation.end()){
                            GEMM(tempResultS, (iteratorFather->second), tempResultF);
                        }
                        else{
                            tempResultF = tempResultS;
                        }
                        //
                        for(I64 tk = 0; tk < 3; ++ tk){
                            for(I64 tm = 0; tm < 3; ++ tm){
                                coo.emplace_back(3 * ti + tk, 3 * col_tj + tm, tempResultF(tk,tm));
                            }
                        }
                    }
                }
                //
                rotateProlong[tl].Coo2Csr(
                    3 * scalarProlong[tl].M, 
                    3 * scalarProlong[tl].N, 
                    coo, 
                    threadManager.levelS2M[tl], 
                    nestLevel + 1);
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

void GeometricMultigrid::ConstraintProlong(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    Log("GeometricMultigrid::ConstraintProlong");
    //
	const I64 maxLevel = dataMesh->maxLevel;
    const I64 numbDmnc = (dataMesh->node2Coordinate).size();
    const std::map<I64,Real>& constrainedDof = dataBoundary->constrainedDof;
    std::vector<I64>& consFlag = dataBoundary->consFlag;
    std::vector<I64>& consPrefSum = dataBoundary->consPrefSum;
    //
    consFlag.assign(3 * numbDmnc, 1);
    for(const auto& [nodeDof, _] : constrainedDof){
        I64 nodeIndex = nodeDof / 3;
        I64 dofIndex = nodeDof % 3;
        consFlag[3 * nodeLevelPosition[nodeIndex][1] + dofIndex] = 0;
    }
    consPrefSum.resize(consFlag.size());
    std::exclusive_scan(
        std::execution::unseq, 
        consFlag.begin(), consFlag.end(),  
        consPrefSum.begin(), 0);
    //
    constraintProlong.resize(maxLevel + 1);
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel + 1, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tl = start_tp; tl < end_tp; ++ tl){
                I64 ti_max = 3 * scalarProlong[tl].N;
                //
                std::vector<Triplet> coo;
                coo.reserve(ti_max);
                for(I64 ti = 0; ti < ti_max; ++ ti){
                    if(consFlag[ti] == 1){
                        coo.emplace_back(ti, consPrefSum[ti], 1.0);
                    }
                }
                //
                constraintProlong[tl].Coo2Csr(
                    ti_max, 
                    consPrefSum[ti_max - 1] + consFlag[ti_max - 1], 
                    coo, 
                    threadManager.levelS2M[tl], 
                    nestLevel + 1);
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

void GeometricMultigrid::Establish(
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){ 
    //
    Log("GeometricMultigrid::Real");
    //
    // Transfer(numbThreads);
    //V[maxLevel + 1] = DOF_ORDERING * DOF_ROTATE * U[maxLevel + 1]
    //V is original xyz, U is reordered and rotated
    if(nestLevel == 0){
        DofRotateOrder(threadManager.one2oneS2S, nestLevel);
    }
    else{
        DofRotateOrder(threadTask, nestLevel);
    }
    RotateProlong(threadTask, nestLevel);
    ConstraintProlong(threadTask, nestLevel);
    //
	const I64 maxLevel = dataMesh->maxLevel;
    (dataPcg->hierarchyStiffness).resize(maxLevel + 1);
    SparseMatrix& solverMatrix = (dataPcg->hierarchyStiffness)[maxLevel];
    const I64 numbDmnc = (dataMesh->node2Coordinate).size();
    AlignedVectorRx& loadVect = dataBoundary->loadVect;
    std::vector<SparseMatrix>& realProlong = dataPcg->realProlong;
    //
    //for i+1 == maxLevel + 1: V[i+1] = rotaProl[i] * consProl[i] * V[i]
    if(nestLevel == 0){
        GEMM(rotateProlong[maxLevel], constraintProlong[maxLevel], 
            maxProlong, threadManager.one2oneS2S, nestLevel);
        SparseMatrix tempROMTS;
        GEMM(rotateOrder, maxProlong, rom, threadManager.one2oneS2S, nestLevel);
        TRANSPOSE(rom, romT);
        GEMM(romT, *dataStiffness, tempROMTS, threadManager.one2oneS2S, nestLevel);
        GEMM(tempROMTS, rom, solverMatrix, threadManager.one2oneS2S, nestLevel);
    }
    else{
        GEMM(rotateProlong[maxLevel], constraintProlong[maxLevel], 
            maxProlong, threadTask, nestLevel);
        SparseMatrix tempROMTS;
        GEMM(rotateOrder, maxProlong, rom, threadTask, nestLevel);
        TRANSPOSE(rom, romT);
        GEMM(romT, *dataStiffness, tempROMTS, threadTask, nestLevel);
        GEMM(tempROMTS, rom, solverMatrix, threadTask, nestLevel);
    }
    //
    AlignedVectorRx tempLoad(3 * numbDmnc, 0.0);
    for(const auto& [dof, value] : dataBoundary->externalForce){
        //externalForce is already in rtz
        tempLoad[3 * nodeLevelPosition[dof / 3][1] + dof % 3] = value;
    }
    loadVect.assign(maxProlong.N, 0.0);
    SparseMatrix maxProlongT;
    TRANSPOSE(maxProlong, maxProlongT);
    if(nestLevel == 0){
        MV(maxProlongT, tempLoad, loadVect, threadManager.one2oneS2S, nestLevel);
    }
    else{
        MV(maxProlongT, tempLoad, loadVect, threadTask, nestLevel);
    }
    //
    //for i+1 <= maxLevel: V[i+1] = realProlong[i] * V[i] 
    //                            = consProl[i+1]^T * rotaProl[i] * consProl[i] * V[i]
    realProlong.resize(maxLevel);
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tl = start_tp; tl < end_tp; ++ tl){
                SparseMatrix tempCTR, tempCpT;
                TRANSPOSE(constraintProlong[tl + 1], tempCpT);
                GEMM(tempCpT, rotateProlong[tl], tempCTR, threadManager.levelS2M[tl], nestLevel + 1);
                GEMM(tempCTR, constraintProlong[tl], realProlong[tl], threadManager.levelS2M[tl], nestLevel + 1);
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

void GeometricMultigrid::DispPost(
    const AlignedVectorRx& inpuDisp, 
    AlignedVectorRx& outpDisp, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel) const{
    //
    MV(rom, inpuDisp, outpDisp, threadTask, nestLevel);
}

void GeometricMultigrid::Output(
    const std::string& directoryPath, 
    const I64& fileIden, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel) const{
    //
	std::ofstream tempOfst_0(directoryPath + "/xyzRtz_" + std::to_string(fileIden) + ".txt", std::ios::out);
    xyzRtz.Output(tempOfst_0);
    tempOfst_0.close();
    //
    tempOfst_0.open(directoryPath + "/maxProlong_" + std::to_string(fileIden) + ".txt", std::ios::out);
    maxProlong.Output(tempOfst_0);
    tempOfst_0.close();
    //
    const I64 maxLevel = dataMesh->maxLevel;
    const std::vector<SparseMatrix>& realProlong = dataPcg->realProlong;
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(maxLevel + 1, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 tl = start_tp; tl < end_tp; ++ tl){
                //
                std::ofstream tempOfst_1(
                    directoryPath + "/scalarProlong_" + std::to_string(fileIden) 
                        + "_" + std::to_string(tl) + ".txt", 
                    std::ios::out);
                (scalarProlong[tl]).Output(tempOfst_1);
                tempOfst_1.close();
                //
                tempOfst_1.open(
                    directoryPath + "/rotateProlong_" + std::to_string(fileIden) 
                        + "_" + std::to_string(tl) + ".txt", 
                    std::ios::out);
                (rotateProlong[tl]).Output(tempOfst_1);
                tempOfst_1.close();
                //
                tempOfst_1.open(
                    directoryPath + "/constraintProlong_" + std::to_string(fileIden) 
                        + "_" + std::to_string(tl) + ".txt", 
                    std::ios::out);
                (constraintProlong[tl]).Output(tempOfst_1);
                tempOfst_1.close();
                //
                if(tl == maxLevel) continue;
                tempOfst_1.open(
                    directoryPath + "/realProlong_" + std::to_string(fileIden) 
                        + "_" + std::to_string(tl) + ".txt", 
                    std::ios::out);
                (realProlong[tl]).Output(tempOfst_1);
                tempOfst_1.close();
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

} // namespace Ddpca