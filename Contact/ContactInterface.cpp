#include "ContactInterface.hpp"
#include "../Contact/CSQuadSegment.hpp"

#include <mutex>

namespace Ddpca{

void ContactInterface::OutputSegments(const std::string& directoryPath, const I64& fileIden) const{
	//
	Log("        ContactInterface::OutputSegments");
    //
	std::ofstream tempOfst(directoryPath + "/resuSegm_" + std::to_string(fileIden) + "_0.txt", std::ios::out);
	const I64 maseSize = masterSegments.size();
	for(I64 ti = 0; ti < maseSize; ++ ti){
		tempOfst << std::setw(10) << masterSegments[ti][0] 
			<< std::setw(10) << masterSegments[ti][1] 
			<< std::setw(10) << masterSegments[ti][2] 
			<< std::setw(10) << masterSegments[ti][3] << "\n";
	}
	tempOfst.close();
    //
	tempOfst.open(directoryPath + "/resuSegm_" + std::to_string(fileIden) + "_1.txt", std::ios::out);
	const I64 slseSize = slaveSegments.size();
	for(I64 ti = 0; ti < slseSize; ++ ti){
		tempOfst << std::setw(10) << slaveSegments[ti][0] 
			<< std::setw(10) << slaveSegments[ti][1] 
			<< std::setw(10) << slaveSegments[ti][2] 
			<< std::setw(10) << slaveSegments[ti][3] << "\n";
	}
	tempOfst.close();
}

void ContactInterface::OutputIntegralPoints(
    const std::string& directoryPath, const I64& fileIden) const{
    //
    Log("        ContactInterface::OutputIntegralPoints");
    //
    std::ofstream tempOfst(directoryPath + "/resuInpo_" + std::to_string(fileIden) + ".txt", std::ios::out);
    const I64 inpoSize = integralPoints.size();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
    for(I64 ti = 0; ti < inpoSize; ++ ti){
		tempOfst << std::setw(30) << (integralPoints[ti].contactPoint)[0][0] 
			<< std::setw(30) << (integralPoints[ti].contactPoint)[0][1] 
			<< std::setw(30) << (integralPoints[ti].contactPoint)[0][2] 
			<< std::setw(30) << (integralPoints[ti].contactPoint)[1][0] 
			<< std::setw(30) << (integralPoints[ti].contactPoint)[1][1] 
			<< std::setw(30) << (integralPoints[ti].contactPoint)[1][2] 
			<< std::setw(30) << integralPoints[ti].initialNormalGap << "\n";
    }
    tempOfst.close();
}

void ContactInterface::BucketSort(
    const std::array<std::vector<Real>, 2>& masterLocal, std::array<I64, 2> divisionNumber){
    //
    assert((divisionNumber[0] > 0 && divisionNumber[1] > 0) && "divisionNumber is zero");
    //diviNumb: refer to the sparse side, but not dense side
    bucket.resize(divisionNumber[0]);
    for(I64 ti = 0; ti < divisionNumber[0]; ++ ti){
        bucket[ti].resize(divisionNumber[1]);
    }
    //
    for(I64 tl = 0; tl < 2; ++ tl){
        // Real min_tl = 1.0E200, max_tl = -1.0E200;
        // I64 maloSize = masterLocal[tl].size();
        // for(I64 ti = 0; ti < maloSize; ++ ti){
        //     min_tl = std::min(min_tl, masterLocal[tl][ti]);
        //     max_tl = std::max(max_tl, masterLocal[tl][ti]);
        // }
        Real min_tl = * std::min_element(std::execution::unseq, 
            masterLocal[tl].begin(), masterLocal[tl].end());
        Real max_tl = * std::max_element(std::execution::unseq, 
            masterLocal[tl].begin(), masterLocal[tl].end());
        Real increment_tl = (max_tl - min_tl) / divisionNumber[tl];
        increment_tl = std::max(increment_tl, 1.0E-10); // not enough, still needs improvement
        bucketLocal[tl][0] = min_tl - increment_tl;
        bucketLocal[tl][1] = max_tl + increment_tl;
        bucketLocal[tl][2] = (bucketLocal[tl][1] - bucketLocal[tl][0]) / divisionNumber[tl];
    }
    //
    const I64 maseSize = masterSegments.size();
    for(I64 ti = 0; ti < maseSize; ++ ti){
        I64 buck_r = (masterLocal[0][ti] - bucketLocal[0][0]) / bucketLocal[0][2];
        I64 buck_c = (masterLocal[1][ti] - bucketLocal[1][0]) / bucketLocal[1][2];
        bucket[buck_r][buck_c].emplace_back(ti);
    }
}

void ContactInterface::LocalSearch(
    const Mesh& masterMesh, 
    const Mesh& slaveMesh, 
    const std::vector<std::array<Real, 8>>& slaveLocal, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
	const I64 nestLevel, 
    Real maxiDist){
    //
    const I64 buckSize_r = bucket.size();
    const I64 buckSize_c = bucket[0].size();
    const I64 slseSize = slaveSegments.size();
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(slseSize, numbPart, startIndex, endIndex);
    std::mutex mutexTemp;
    std::function<void(I64)> taskFunction_0 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            Real minXi, maxXi, minEta, maxEta;
            I64 buck_rmin, buck_rmax, buck_cmin, buck_cmax, tj_min, tj_max, tk_min, tk_max;
            CSQuadSegment tempCs;
            std::array<I64, 4> tempMasterSegment, tempSlaveSegment;
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                minXi = 1.0E20;
                maxXi = -1.0E20;
                minEta = 1.0E20;
                maxEta = -1.0E20;
                for(I64 tc = 0; tc < 4; ++ tc){
                    minXi = std::min(minXi, slaveLocal[ti][tc * 2 + 0]);
                    maxXi = std::max(maxXi, slaveLocal[ti][tc * 2 + 0]);
                    minEta = std::min(minEta, slaveLocal[ti][tc * 2 + 1]);
                    maxEta = std::max(maxEta, slaveLocal[ti][tc * 2 + 1]);
                }
                buck_rmin = (minXi - bucketLocal[0][0]) / bucketLocal[0][2];
                buck_rmax = (maxXi - bucketLocal[0][0]) / bucketLocal[0][2];
                buck_cmin = (minEta - bucketLocal[1][0]) / bucketLocal[1][2];
                buck_cmax = (maxEta - bucketLocal[1][0]) / bucketLocal[1][2];
                tj_min = std::max(buck_rmin - 1, (I64)0);
                tj_max = std::min(buck_rmax + 1, (I64)(buckSize_r - 1));
                tk_min = std::max(buck_cmin - 1, (I64)0);
                tk_max = std::min(buck_cmax + 1, (I64)(buckSize_c - 1));
                for(I64 tj = tj_min; tj <= tj_max; tj ++){
                    for(I64 tk = tk_min; tk <= tk_max; tk ++){
                        for(const auto &tempMaster : bucket[tj][tk]){
                            //
                            tempCs.integralPoints.clear();
                            tempMasterSegment = masterSegments[tempMaster];
                            tempSlaveSegment = slaveSegments[ti];
                            const std::map<I64, Coordinate>& masterNode2Coordinate = masterMesh.node2Coordinate;
                            const std::map<I64, Coordinate>& slaveNode2Coordinate = slaveMesh.node2Coordinate;
                            for(I64 ti = 0; ti < 4; ti ++){
                                auto iterNoco = masterNode2Coordinate.find(tempMasterSegment[ti]);
                                tempCs.masterCorners[ti] = {(iterNoco->second)[0], 
                                    (iterNoco->second)[1], (iterNoco->second)[2]};
                                iterNoco = slaveNode2Coordinate.find(tempSlaveSegment[ti]);
                                tempCs.slaveCorners[ti] = {(iterNoco->second)[0], 
                                    (iterNoco->second)[1], (iterNoco->second)[2]};
                            }
                            tempCs.Search(tempMasterSegment, tempSlaveSegment);
                            //
                            bool tempFlag = false;
                            for(const auto &iterInpo : tempCs.integralPoints){
                                if(iterInpo.initialNormalGap <= maxiDist){
                                    tempFlag = true;
                                    break;
                                }
                            }
                            if(tempFlag){
                                std::lock_guard<std::mutex> lock(mutexTemp);
                                integralPoints.insert(integralPoints.end(), 
                                    tempCs.integralPoints.begin(), tempCs.integralPoints.end());
                            }
                        }
                    }
                }
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
    assert(!integralPoints.empty() && "ContactInterface::LocalSearch empty");
    // basisVector[0~2] must be perpendicular to each other
    const I64 inpoSize = integralPoints.size();
    EvenlyDistribute(inpoSize, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction_1 = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            for(I64 ti = start_tp; ti < end_tp; ++ ti){
                IntegralPoint& tempInpo = integralPoints[ti];
                Real tempProj = DOT(tempInpo.basisVector[1], tempInpo.basisVector[2]);
                AXPY(- tempProj, tempInpo.basisVector[1], tempInpo.basisVector[2]);
                SCAL(1.0 / NRM2(tempInpo.basisVector[2]), tempInpo.basisVector[2]);
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
    auto iteratorMax = std::max_element(std::execution::unseq, 
        integralPoints.begin(), integralPoints.end(), 
        [](const IntegralPoint &a, const IntegralPoint &b){
            return a.initialNormalGap < b.initialNormalGap;
        });
    auto iteratorMin = std::min_element(std::execution::unseq, 
        integralPoints.begin(), integralPoints.end(), 
        [](const IntegralPoint &a, const IntegralPoint &b){
            return a.initialNormalGap < b.initialNormalGap;
        });
    Log("        ContactInterface::LocalSearch gap: minimum = " 
        + Double2String(iteratorMin->initialNormalGap)
        + ", maximum = " + Double2String(iteratorMax->initialNormalGap));
    // if(frictionCoefficient == 0.0){
    //     InterpolateBasisGap(masterMesh, threadTask, nestLevel);
    // }
}

void ContactInterface::InterpolateBasisGap(const Mesh& masterMesh, 
    const std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const I64 nestLevel){
    //
    std::map<I64,std::array<Real, 3>> nodeBasis;
    const BilinearQuadrature<2>& biliQuad = GetBilinearQuadrature();
    //
    const I64 maseSize = masterSegments.size();
    DenseMatrix factMatr(3, 4);
    std::array<std::array<Real, 3>, 2> PrmaPxie;
    std::array<Real, 3> tempNorm;
    for(I64 ti = 0; ti < maseSize; ++ ti){
        const std::array<I64, 4>& tiMasterSegment = masterSegments[ti];
        factMatr.Fill(0.0);
        for(I64 tj = 0; tj < 4; ++ tj){
            auto iterNoco = masterMesh.node2Coordinate.find(tiMasterSegment[tj]);
            for(I64 tk = 0; tk < 3; ++ tk){
                Real tkCoor = (iterNoco->second)[tk];
                factMatr(tk, 0) += tkCoor / 4.0;
                factMatr(tk, 1) += tkCoor * biliQuad.cornerNodes[tj][0] / 4.0;
                factMatr(tk, 2) += tkCoor * biliQuad.cornerNodes[tj][1] / 4.0;
                factMatr(tk, 3) += tkCoor 
                    * biliQuad.cornerNodes[tj][0] * biliQuad.cornerNodes[tj][1] / 4.0;
            }
        }
        for(I64 tj = 0; tj < 4; ++ tj){
            for(I64 tk = 0; tk < 3; ++ tk){
                PrmaPxie[0][tk] = factMatr(tk, 1) + factMatr(tk, 3) * biliQuad.cornerNodes[tj][1];
                PrmaPxie[1][tk] = factMatr(tk, 2) + factMatr(tk, 3) * biliQuad.cornerNodes[tj][0];
            }
            tempNorm = Cross(PrmaPxie[0], PrmaPxie[1]);
            SCAL(1.0 / NRM2(tempNorm), tempNorm);
            I64 tjNode = tiMasterSegment[tj];
            auto iterNoba = nodeBasis.find(tjNode);
            if(iterNoba == nodeBasis.end()){
                nodeBasis.emplace(tjNode, tempNorm);
            }
            else{
                XPEY(iterNoba->second, tempNorm);
            }
        }
    }
    for(auto& iterNoba : nodeBasis){
        NORMALIZE(iterNoba.second);
    }
    //
    const I64 inpoSize = integralPoints.size();
    I64 numbPart = (nestLevel <  threadManager.maxNestLevel) ? threadTask.size() : 1;
    std::vector<I64> startIndex(numbPart), endIndex(numbPart);
    EvenlyDistribute(inpoSize, numbPart, startIndex, endIndex);
    std::function<void(I64)> taskFunction = 
        [&](I64 tp){
            I64 start_tp = startIndex[tp];
            I64 end_tp = endIndex[tp];
            std::array<Real, 3> tempNorm, tempCross, tempDiff;
            for(I64 tt = start_tp; tt < end_tp; ++ tt){
                //
                IntegralPoint& iterInpo = integralPoints[tt];
                tempNorm.fill(0.0);
                for(I64 ti = 0; ti < 4; ++ ti){
                    I64 tjNode = iterInpo.node[0][ti];
                    AXPY(iterInpo.shapeFunction[0][ti], nodeBasis[tjNode], tempNorm);
                }
                NORMALIZE(tempNorm);
                iterInpo.basisVector[0] = tempNorm;
                //
                iterInpo.basisVector[1].fill(0.0);
                iterInpo.basisVector[2].fill(0.0);
                Real d1x = tempNorm[0];
                Real d1y = tempNorm[1];
                Real d1z = tempNorm[2];
                Real EPSI = 1.0E-14;
                if(std::abs(d1y) < EPSI){
                    iterInpo.basisVector[1][1] = 1.0;
                    if(std::abs(d1x) < EPSI){
                        iterInpo.basisVector[2][0] = 1.0;
                    }
                    else{
                        if(std::abs(d1z) < EPSI){
                            iterInpo.basisVector[2][2] = 1.0;
                        }
                        else{
                            iterInpo.basisVector[2][0] = d1z / std::sqrt(std::pow(d1z, 2.0) + std::pow(d1x, 2.0));
                            iterInpo.basisVector[2][2] = - d1x / std::sqrt(std::pow(d1z, 2.0) + std::pow(d1x, 2.0));
                        }
                    }
                }
                else if(std::abs(d1z) < EPSI){
                    iterInpo.basisVector[1][2] = 1.0;
                    if(std::abs(d1x) < EPSI){
                        iterInpo.basisVector[2][0] = 1.0;
                    }
                    else{
                        iterInpo.basisVector[2][0] = d1y / std::sqrt(std::pow(d1y, 2.0) + std::pow(d1x, 2.0));
                        iterInpo.basisVector[2][1] = - d1x / std::sqrt(std::pow(d1y, 2.0) + std::pow(d1x, 2.0));
                    }
                }
                else{
                    Real b2c2 = std::pow(d1y, 2.0) + std::pow(d1z, 2.0);
                    iterInpo.basisVector[1][1] = d1z / std::sqrt(b2c2);
                    iterInpo.basisVector[1][2] = - d1y / std::sqrt(b2c2);
                    Real abc = std::sqrt(std::pow(b2c2, 2.0) + std::pow(d1x * d1y, 2.0) + std::pow(d1x * d1z, 2.0));
                    iterInpo.basisVector[2][0] = b2c2 / abc;
                    iterInpo.basisVector[2][1] = (- d1x * d1y) / abc;
                    iterInpo.basisVector[2][2] = (- d1x * d1z) / abc;
                }
                tempCross = Cross(iterInpo.basisVector[0], iterInpo.basisVector[1]);
                Real tempDot = DOT(tempCross, iterInpo.basisVector[2]);
                if(tempDot < 0.0){
                    SCAL(-1.0, iterInpo.basisVector[2]);
                }
                //
                tempDiff = iterInpo.contactPoint[1];
                AXPY(-1.0, iterInpo.contactPoint[0], tempDiff);
                iterInpo.initialNormalGap = DOT(iterInpo.basisVector[0], tempDiff);
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

void ContactInterface::AdaptiveRefine(Mesh& masterMesh, Mesh& slaveMesh, bool &isRefi, 
    CurvedSurface& mastSurf, CurvedSurface& slavSurf, I64 tempLeve, 
    Real distCrit, std::array<I64, 2> divisionNumber, 
    std::function<void(Coordinate, Real &, Real &)> CART_CURV){
    //
    masterSegments.clear();
    mastSurf.Initialize();
    while(mastSurf.Increment(masterMesh)){
		//NO NEED: if((* masterMesh).elemVect[iterEfsu_0.eid].level == tempLeve){
        masterSegments.emplace_back(mastSurf.currentFace);
    }
    slaveSegments.clear();
    slavSurf.Initialize();
	while(slavSurf.Increment(slaveMesh)){
		slaveSegments.emplace_back(slavSurf.currentFace);
	}
	//
    Ddpca::I64 maseSize = masterSegments.size();
    std::array<std::vector<Ddpca::Real>, 2> masterLocal;
    masterLocal[0].resize(maseSize);
    masterLocal[1].resize(maseSize);
    Ddpca::I64 node_tj;
    Ddpca::Real tempXi, tempEta;
	for(long ti = 0; ti < maseSize; ti ++){
		tempXi = 0.0;
		tempEta = 0.0;
		for(long tj = 0; tj < 4; tj ++){
			node_tj = masterSegments[ti][tj];
			auto iterNoco = masterMesh.node2Coordinate.find(node_tj);
			CART_CURV(iterNoco->second, tempXi, tempEta);
		}
        masterLocal[0][ti] = tempXi / 4.0;
        masterLocal[1][ti] = tempEta / 4.0;
	}
    BucketSort(masterLocal, divisionNumber);
    Ddpca::I64 slseSize = slaveSegments.size();
    std::vector<std::array<Ddpca::Real, 8>> slaveLocal(slseSize);
	for(long ti = 0; ti < slseSize; ti ++){
		for(long tj = 0; tj < 4; tj ++){
            node_tj = slaveSegments[ti][tj];
            auto iterNoco = slaveMesh.node2Coordinate.find(node_tj);
            tempXi = 0.0;
            tempEta = 0.0;
			CART_CURV(iterNoco->second, tempXi, tempEta);
            slaveLocal[ti][tj * 2 + 0] = tempXi;
            slaveLocal[ti][tj * 2 + 1] = tempEta;
		}
	}
	//
	std::array<std::set<long>,2> spliNode;
    const I64 buckSize_r = bucket.size();
    const I64 buckSize_c = bucket[0].size();
    Real minXi, maxXi, minEta, maxEta;
    I64 buck_rmin, buck_rmax, buck_cmin, buck_cmax, tj_min, tj_max, tk_min, tk_max;
    CSQuadSegment tempCs;
    std::array<I64, 4> tempMasterSegment, tempSlaveSegment;
    const std::map<I64, Coordinate>& masterNode2Coordinate = masterMesh.node2Coordinate;
    const std::map<I64, Coordinate>& slaveNode2Coordinate = slaveMesh.node2Coordinate;
    for(I64 ti = 0; ti < slseSize; ++ ti){
        minXi = 1.0E20;
        maxXi = -1.0E20;
        minEta = 1.0E20;
        maxEta = -1.0E20;
        for(I64 tc = 0; tc < 4; ++ tc){
            minXi = std::min(minXi, slaveLocal[ti][tc * 2 + 0]);
            maxXi = std::max(maxXi, slaveLocal[ti][tc * 2 + 0]);
            minEta = std::min(minEta, slaveLocal[ti][tc * 2 + 1]);
            maxEta = std::max(maxEta, slaveLocal[ti][tc * 2 + 1]);
        }
        buck_rmin = (minXi - bucketLocal[0][0]) / bucketLocal[0][2];
        buck_rmax = (maxXi - bucketLocal[0][0]) / bucketLocal[0][2];
        buck_cmin = (minEta - bucketLocal[1][0]) / bucketLocal[1][2];
        buck_cmax = (maxEta - bucketLocal[1][0]) / bucketLocal[1][2];
        tj_min = std::max(buck_rmin - 1, (I64)0);
        tj_max = std::min(buck_rmax + 1, (I64)(buckSize_r - 1));
        tk_min = std::max(buck_cmin - 1, (I64)0);
        tk_max = std::min(buck_cmax + 1, (I64)(buckSize_c - 1));
        for(I64 tj = tj_min; tj <= tj_max; tj ++){
            for(I64 tk = tk_min; tk <= tk_max; tk ++){
                for(const auto &tempMaster : bucket[tj][tk]){
                    //
                    tempCs.integralPoints.clear();
                    tempMasterSegment = masterSegments[tempMaster];
                    tempSlaveSegment = slaveSegments[ti];
                    for(I64 ti = 0; ti < 4; ti ++){
                        auto iterNoco = masterNode2Coordinate.find(tempMasterSegment[ti]);
                        tempCs.masterCorners[ti] = {(iterNoco->second)[0], 
                            (iterNoco->second)[1], (iterNoco->second)[2]};
                        iterNoco = slaveNode2Coordinate.find(tempSlaveSegment[ti]);
                        tempCs.slaveCorners[ti] = {(iterNoco->second)[0], 
                            (iterNoco->second)[1], (iterNoco->second)[2]};
                    }
                    tempCs.Search(tempMasterSegment, tempSlaveSegment);
                    //
                    bool tempFlag = false;
                    for(const auto &iterInpo : tempCs.integralPoints){
                        if(iterInpo.initialNormalGap <= distCrit){
                            tempFlag = true;
                            break;
                        }
                    }
                    if(tempFlag){
						for(const auto &iterInpo : tempCs.integralPoints){
							for(long tm = 0; tm < 4; tm ++){
								spliNode[0].emplace(iterInpo.node[0][tm]);
								spliNode[1].emplace(iterInpo.node[1][tm]);
							}
						}
                    }
                }
            }
        }
    }
	if(spliNode[0].size() == 0 && spliNode[1].size() == 0){
		isRefi = false;
		return;
	}
	else{
		isRefi = true;
	}
	//
    std::set<I64> elementsToSplit;
    std::map<I64, std::set<I64>> subElements;
    std::map<std::vector<I64>, Coordinate> curvInte;
	for(long tv = 0; tv < 2; tv ++){
		Mesh& tempMesh = (tv == 0) ? masterMesh : slaveMesh;
		CurvedSurface& tempSurf = (tv == 0) ? mastSurf : slavSurf;
		//
		elementsToSplit.clear();
        const I64 teelSize = tempMesh.elements.size();
		for(long ti = 0; ti < teelSize; ti ++){
			auto &tempElem = tempMesh.elements[ti];
			if(tempElem.children.size() > 0 || tempElem.level != tempLeve){
				continue;
			}
			bool tempFlag = false;
			for(long tj = 0; tj < 8; tj ++){
				auto iterSpno = spliNode[tv].find(tempElem.cornerNodes[tj]);
				if(iterSpno != spliNode[tv].end()){
					tempFlag = true;
					break;
				}
			}
			if(tempFlag == true){
				elementsToSplit.insert(ti);
				tempElem.refinementPattern = OctreeElement::REFINEMENT_FULL;
			}
		}
		curvInte.clear();
        tempSurf.Refine(tempMesh, elementsToSplit, curvInte);
        tempMesh.Refine(elementsToSplit, subElements, curvInte);
	}
}

} // namespace Ddpca