
#ifndef _CURVEDS_H
#define _CURVEDS_H

#include "PREP.h"

//curved surface
class CURVEDS{
public:
	//indice to point
	std::vector<std::vector<COOR>> indiPoin;
	//point to indice
	std::map<COOR,std::vector<long>> poinIndi;
	long INSERT(long ti,long tj, COOR inpuCoor);
	//output: 1 - on the curved surface, outpCoor is valid
	//otherwise - not on the curved surface, outpCoor is not valid
	long REFINE_SEARCH(std::vector<COOR> inpuCoor, COOR &outpCoor);
	//if one face of an element is on CURVEDS, this element is gonna be refined
	//(put into planSurf)
	long REFINE(const std::set<long> &spliElem, 
		const std::vector<TREE_ELEM> &elemVect, 
		const std::map<long,COOR> &nodeCoor, 
		std::map<std::vector<long>,COOR> &planSurf
	);
	//rotation and translation of rigid body
	long RIGI_ROTR(Eigen::Matrix3d rotaMatr, Eigen::Vector3d tranVect);
};

long CURVEDS::INSERT(long ti,long tj, COOR inpuCoor){
	//not considering duplicated key
	indiPoin[ti][tj] = inpuCoor;
	std::vector<long> tempVect = {ti, tj};
	poinIndi.emplace(inpuCoor, tempVect);
	return 1;
}

long CURVEDS::REFINE_SEARCH(std::vector<COOR> inpuCoor, COOR &outpCoor){
	bool tempFlag = true;
	std::vector<long> resuIndi(2, 0);
	for(long ti = 0; ti < inpuCoor.size(); ti ++){
		auto iterPoin = poinIndi.find(inpuCoor[ti]);
		if(iterPoin == poinIndi.end()){
			tempFlag = false;
			break;
		}
		resuIndi[0] += (iterPoin->second)[0];
		resuIndi[1] += (iterPoin->second)[1];
	}
	if(tempFlag == false){
		return -1;
	}
	resuIndi[0] /= inpuCoor.size();
	resuIndi[1] /= inpuCoor.size();
	outpCoor = indiPoin[resuIndi[0]][resuIndi[1]];
	return 1;
}

long CURVEDS::REFINE(const std::set<long> &spliElem, 
	const std::vector<TREE_ELEM> &elemVect, 
	const std::map<long,COOR> &nodeCoor, 
	std::map<std::vector<long>,COOR> &planSurf){
	for(const auto &iterSpel : spliElem){
		//
		for(long tj = 0; tj < hexaLine.size(); tj ++){
			std::vector<long> inpuNode(hexaLine[tj].size());
			std::vector<COOR> inpuCoor(hexaLine[tj].size());
			for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
				inpuNode[tk] = elemVect[iterSpel].cornNode[hexaLine[tj][tk]];
				auto iterNoco = nodeCoor.find(inpuNode[tk]);
				inpuCoor[tk] = iterNoco->second;
			}
			COOR outpCoor;
			long searResu = REFINE_SEARCH(inpuCoor, outpCoor);
			if(searResu == 1){
				sort(inpuNode.begin(), inpuNode.end());
				planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
					inpuNode, outpCoor
				));
			}
		}
		//
		for(long tj = 0; tj < hexaFace.size(); tj ++){
			std::vector<long> inpuNode(hexaFace[tj].size());
			std::vector<COOR> inpuCoor(hexaFace[tj].size());
			for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
				inpuNode[tk] = elemVect[iterSpel].cornNode[hexaFace[tj][tk]];
				auto iterNoco = nodeCoor.find(inpuNode[tk]);
				inpuCoor[tk] = iterNoco->second;
			}
			COOR outpCoor;
			long searResu = REFINE_SEARCH(inpuCoor, outpCoor);
			if(searResu == 1){
				sort(inpuNode.begin(), inpuNode.end());
				planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
					inpuNode, outpCoor
				));
			}
		}
	}
	return 1;
}

long CURVEDS::RIGI_ROTR(Eigen::Matrix3d rotaMatr, Eigen::Vector3d tranVect){
	poinIndi.clear();
	for(long ti = 0; ti < indiPoin.size(); ti ++){
		for(long tj = 0; tj < indiPoin[ti].size(); tj ++){
			if(indiPoin[ti][tj].size() == 0){
				continue;
			}
			Eigen::Vector3d tempCoor;
			tempCoor << indiPoin[ti][tj][0], indiPoin[ti][tj][1], indiPoin[ti][tj][2];
			tempCoor = rotaMatr * tempCoor + tranVect;
			indiPoin[ti][tj] = {tempCoor(0), tempCoor(1), tempCoor(2)};
			std::vector<long> tempVect = {ti, tj};
			poinIndi.emplace(indiPoin[ti][tj], tempVect);
		}
	}
	return 1;
}

#endif
