
#ifndef _CSEARCH_H
#define _CSEARCH_H

#include "CURVEDS.h"
#include "MULTIGRID.h"

/*************************************************************************************************/
//Projecting 3d faces onto 2d plane: Cartesian xyz TO curvilinear xi,eta.
//Two faces in 2d are intersecting, the minimum value of the intersecting area.
//If the area of a triangle is less than miniArea, the triangle is degenerated to a line-segment.
const double miniArea = 1.0E-12;

//Quadrature over master side, one quadrature point = one INTEGRAL_POINT.
//Project slave face onto master face (curvilinear plane), a polygon is obtained,
//divide the polygon into sub triangles by its centroid,
//2*2 Guass quadrature over one sub triangle,
//one sub triangle = four quadrature points = four INTEGRAL_POINT
class INTEGRAL_POINT{
public:
	VECTOR2L node;//four nodes: one element face
	VECTOR2D shapFunc;//four shape/basis functions
	Eigen::VectorXd dual;//four dual basis functions on non-mortar side
	std::vector<Eigen::Vector3d> contPoin;
	//basis vector, 0 - normal of master side, 1/2 - tangential
	std::vector<Eigen::Vector3d> basiVect;
	double initNgap;//= normVect * (contPoin 1 - contPoin 0)
	double quadWeig;//quadrature weight
	INTEGRAL_POINT() : node(2,std::vector<long>(4)), shapFunc(2,std::vector<double>(4)), 
		contPoin(2), basiVect(3){
	}
};

//element face of *targGrid locates on surface *targSurf
class EFACE_SURFACE{
public:
	const MULTIGRID *targGrid;
	const CURVEDS *targSurf;
	long ti;
	long tj;
	std::vector<long> currNode;
	long eid;
	EFACE_SURFACE(const MULTIGRID *inpuGrid, const CURVEDS *inpuSurf){
		targGrid = inpuGrid;
		targSurf = inpuSurf;
		ti = 0;
		tj = 0;
		currNode.resize(4);
	}
	long INCREMENT(){
		for(; ti < (*targGrid).elemVect.size(); ti ++){
			if((*targGrid).elemVect[ti].children.size() > 0){
				tj = 0;
				continue;
			}
			for(; tj < hexaFace.size(); tj ++){
				bool tempFlag = true;
				for(long tk = 0; tk < 4; tk ++){
					currNode[tk] = (*targGrid).elemVect[ti].cornNode[hexaFace[tj][tk]];
					auto iterNoco = (*targGrid).nodeCoor.find(currNode[tk]);
					auto iterPoin = (*targSurf).poinIndi.find(iterNoco->second);
					if(iterPoin == (*targSurf).poinIndi.end()){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == true){
					eid = ti;
					tj ++;
					if(tj == hexaFace.size()){
						tj = 0;
						ti ++;
					}
					return 1;
				}
			}
			tj = 0;
		}
		return -1;
	}
};

/*************************************************************************************************/
//contact search
class CSEARCH{
public:
	/*********************************************************************************************/
	const MULTIGRID *mastGrid;//two-body contact
	const MULTIGRID *slavGrid;
	//n*4, n segments/element faces, 4 nodes per segment
	VECTOR2L mastSegm;
	VECTOR2L slavSegm;
	std::vector<INTEGRAL_POINT> intePoin;//the main output
	long COPY(const CSEARCH &targCsea);//copy from targCsea
	long OUTPUT_COSE(long contIden);//output mastSegm and slavSegm
	/*********************************************************************************************/
	//global contact search, bucket search
	//2d local coordinate, 0 - infimum, 1 - supremum, 2 - increment
	VECTOR2D buckLoco;
	//sort master segments into different buckets
	long BUCKET_SORT(const VECTOR2D &mastCoor, std::vector<long> diviNumb);
	VECTOR3L buck;//bucket
	/*********************************************************************************************/
	//project master face point mastXiet/mastPoin onto slave face slavCorn
	long PROJECT_MTS(const std::vector<Eigen::Vector3d> &mastCorn, 
		const std::vector<Eigen::Vector3d> &slavCorn, 
		const Eigen::Vector2d mastXiet, const Eigen::Vector3d mastPoin, 
		Eigen::Matrix<double,3,2> &PrmaPxie, 
		Eigen::Matrix<double,2,1> &slavXiet, double &ngap
	);
	//subroutine of PROJECT_STM
	long PROJECT_STM_SUB(const std::vector<Eigen::Vector3d> &mastCorn, 
		const Eigen::Vector3d &slavPoin, 
		Eigen::Matrix<double,3,2> &PrmaPxie, Eigen::Matrix<double,2,1> &mastXiet, double &ngap
	);
	//project slave point slavPoin onto master face mastCorn
	long PROJECT_STM(const std::vector<Eigen::Vector3d> &mastCorn, 
		const Eigen::Vector3d &slavPoin, 
		Eigen::Matrix<double,3,2> &PrmaPxie, Eigen::Matrix<double,2,1> &mastXiet, double &ngap
	);
	//project master face point mastXiet onto slave face slavCorn: maslPoin[0] to maslPoin[1]
	long MAST_TO_SLAV(const std::vector<Eigen::Vector3d> &mastCorn, 
		const std::vector<Eigen::Vector3d> &slavCorn, const Eigen::Vector2d mastXiet, 
		VECTOR2D &maslShap, std::vector<Eigen::Vector3d> &maslPoin, 
		std::vector<Eigen::Vector3d> &basiVect, double &weigFact
	);
	//area of 2d triangle
	double TRIANGLE_AREA_2D(Eigen::Vector2d tempPoin_0, 
		Eigen::Vector2d tempPoin_1, Eigen::Vector2d tempPoin_2
	);
	//Guass quadrature over a triangle (tempXieta_0~2)
	long TRIANGLE_QUADRATURE(Eigen::Vector2d tempXiet_0, 
		Eigen::Vector2d tempXiet_1, Eigen::Vector2d tempXiet_2, 
		std::vector<Eigen::Vector2d> &listXiet, std::vector<double> &listWeig
	);
	//sort Eigen::Vector2d by the tempIden-th component
	long SORT_BY_2D(Eigen::Vector2d &tempPoin_0, Eigen::Vector2d &tempPoin_1, long tempIden);
	//whether line tempPoin_0~1 intersect with line tempPoin_2~3
	long IS_CROSS_2D(Eigen::Vector2d tempPoin_0, Eigen::Vector2d tempPoin_1, 
		Eigen::Vector2d tempPoin_2, Eigen::Vector2d tempPoin_3
	);
	//the intersection of two lines
	long LINE_INTERSECT_2D(Eigen::Vector2d tempPoin_0, Eigen::Vector2d tempPoin_1, 
		Eigen::Vector2d tempPoin_2, Eigen::Vector2d tempPoin_3, 
		std::vector<Eigen::Vector2d> &resuInte
	);
	//whether tempPoin in 2d face tempCorn
	long IN_CQUAD_2D(Eigen::Vector2d tempPoin, const std::vector<Eigen::Vector2d> tempCorn);
	//subroutine of SEGMENT_INTERSECT
	long SI_SUB(const std::vector<Eigen::Vector3d> &mastCorn, 
		const std::vector<Eigen::Vector3d> &slavCorn, 
		std::vector<Eigen::Vector2d> &listXiet, std::vector<double> &listWeig
	);
	//intersection of segment tempMast and segment tempSlav
	long SEGMENT_INTERSECT(long tempMast, long tempSlav, std::vector<INTEGRAL_POINT> &inpuInpo);
	long CONTACT_SEARCH(const VECTOR2D &slavCoor, double maxiDist = 1.0E12);
	//output intePoin
	long OUTPUT_INPO(long contIden);
	/*********************************************************************************************/
	//refine faces by distance distCrit
	long ADAPTIVE_REFINE(MULTIGRID *tempMast, MULTIGRID *tempSlav, bool &isnoRefi, 
		CURVEDS *mastSurf, CURVEDS *slavSurf, long tempLeve, 
		double distCrit, std::vector<long> buckNumb, 
		std::function<void(COOR, double &, double &)> CART_CURV
	);
};

long CSEARCH::COPY(const CSEARCH &targCsea){
	//NO: mastGrid, slavGrid
	mastSegm = targCsea.mastSegm;
	slavSegm = targCsea.slavSegm;
	intePoin = targCsea.intePoin;
	buckLoco = targCsea.buckLoco;
	buck = targCsea.buck;
	return 1;
}

long CSEARCH::OUTPUT_COSE(long contIden){
	std::stringstream tempStre;
	tempStre << "resuSegm_" << contIden << "_0.txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	for(long ti = 0; ti < mastSegm.size(); ti ++){
		tempOfst << std::setw(10) << mastSegm[ti][0] 
			<< std::setw(10) << mastSegm[ti][1] 
			<< std::setw(10) << mastSegm[ti][2] 
			<< std::setw(10) << mastSegm[ti][3] << std::endl;
	}
	tempOfst.close();
	tempStre << "resuSegm_" << contIden << "_1.txt";
	tempOfst.open(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	for(long ti = 0; ti < slavSegm.size(); ti ++){
		tempOfst << std::setw(10) << slavSegm[ti][0] 
			<< std::setw(10) << slavSegm[ti][1] 
			<< std::setw(10) << slavSegm[ti][2] 
			<< std::setw(10) << slavSegm[ti][3] << std::endl;
	}
	tempOfst.close();
	return 1;
}

long CSEARCH::BUCKET_SORT(const VECTOR2D &mastCoor, std::vector<long> diviNumb){
	//diviNumb: refer to the sparse side, but not dense side
	buckLoco.resize(2);
	buckLoco[0].resize(3);
	buckLoco[1].resize(3);
	VECT_RESI(buck, diviNumb[0], diviNumb[1]);
	//
	for(long ti = 0; ti < 2; ti ++){
		auto mini_ti = std::min_element(mastCoor[ti].begin(), mastCoor[ti].end());
		auto maxi_ti = std::max_element(mastCoor[ti].begin(), mastCoor[ti].end());
		double incr_ti;
		incr_ti = (*maxi_ti - *mini_ti) / diviNumb[ti];
		if(abs(incr_ti) < 1.0E-10){
			incr_ti = 1.0E-10;//not enough, still needs improvement
		}
		buckLoco[ti][0] = *mini_ti - incr_ti;
		buckLoco[ti][1] = *maxi_ti + incr_ti;
		buckLoco[ti][2] = (buckLoco[ti][1] - buckLoco[ti][0]) / diviNumb[ti];
	}
	for(long ti = 0; ti < mastSegm.size(); ti ++){
		long buck_r = (mastCoor[0][ti] - buckLoco[0][0]) / buckLoco[0][2];
		long buck_c = (mastCoor[1][ti] - buckLoco[1][0]) / buckLoco[1][2];
		buck[buck_r][buck_c].emplace_back(ti);
	}
	return 1;
}

long CSEARCH::PROJECT_MTS(const std::vector<Eigen::Vector3d> &mastCorn, 
	const std::vector<Eigen::Vector3d> &slavCorn, 
	const Eigen::Vector2d mastXiet, const Eigen::Vector3d mastPoin, 
	Eigen::Matrix<double,3,2> &PrmaPxie, 
	Eigen::Matrix<double,2,1> &slavXiet, double &ngap){
	PrmaPxie = Eigen::MatrixXd::Zero(3,2);
	for(long ti = 0; ti < 3; ti ++){
		for(long tj = 0; tj < 4; tj ++){
			PrmaPxie(ti,0) += mastCorn[tj](ti) * (biliQuad.nacoCorn[tj][0] / 4.0 
				+ biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] * mastXiet(1) / 4.0);
			PrmaPxie(ti,1) += mastCorn[tj](ti) * (biliQuad.nacoCorn[tj][1] / 4.0 
				+ biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] * mastXiet(0) / 4.0);
		}
	}
	Eigen::Matrix<double,3,4> factMatr = Eigen::MatrixXd::Zero(3,4);
	for(long ti = 0; ti < 3; ti ++){
		for(long tj = 0; tj < 4; tj ++){
			factMatr(ti,0) += slavCorn[tj](ti) / 4.0;
			factMatr(ti,1) += slavCorn[tj](ti) * biliQuad.nacoCorn[tj][0] / 4.0;
			factMatr(ti,2) += slavCorn[tj](ti) * biliQuad.nacoCorn[tj][1] / 4.0;
			factMatr(ti,3) += slavCorn[tj](ti) 
				* biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] / 4.0;
		}
	}
	Eigen::Matrix<double,2,4> equaFact;
	equaFact(0,0) = factMatr(0,0) * PrmaPxie(0,0) + factMatr(1,0) * PrmaPxie(1,0) 
		+ factMatr(2,0) * PrmaPxie(2,0) - mastPoin(0) * PrmaPxie(0,0) 
		- mastPoin(1) * PrmaPxie(1,0) - mastPoin(2) * PrmaPxie(2,0);
	equaFact(0,1) = factMatr(0,1) * PrmaPxie(0,0) + factMatr(1,1) * PrmaPxie(1,0) 
		+ factMatr(2,1) * PrmaPxie(2,0);
	equaFact(0,2) = factMatr(0,2) * PrmaPxie(0,0) + factMatr(1,2) * PrmaPxie(1,0) 
		+ factMatr(2,2) * PrmaPxie(2,0);
	equaFact(0,3) = factMatr(0,3) * PrmaPxie(0,0) + factMatr(1,3) * PrmaPxie(1,0) 
		+ factMatr(2,3) * PrmaPxie(2,0);
	equaFact(1,0) = factMatr(0,0) * PrmaPxie(0,1) + factMatr(1,0) * PrmaPxie(1,1) 
		+ factMatr(2,0) * PrmaPxie(2,1) - mastPoin(0) * PrmaPxie(0,1) 
		- mastPoin(1) * PrmaPxie(1,1) - mastPoin(2) * PrmaPxie(2,1);
	equaFact(1,1) = factMatr(0,1) * PrmaPxie(0,1) + factMatr(1,1) * PrmaPxie(1,1) 
		+ factMatr(2,1) * PrmaPxie(2,1);
	equaFact(1,2) = factMatr(0,2) * PrmaPxie(0,1) + factMatr(1,2) * PrmaPxie(1,1) 
		+ factMatr(2,2) * PrmaPxie(2,1);
	equaFact(1,3) = factMatr(0,3) * PrmaPxie(0,1) + factMatr(1,3) * PrmaPxie(1,1) 
		+ factMatr(2,3) * PrmaPxie(2,1);
	slavXiet = Eigen::MatrixXd::Zero(2,1);
	long tj;
	for(tj = 0; tj < 60; tj ++){
		Eigen::Matrix<double,2,1> resiResu;
		resiResu(0) = equaFact(0,0) + equaFact(0,1) * slavXiet(0) + equaFact(0,2) * slavXiet(1) 
			+ equaFact(0,3) * slavXiet(0) * slavXiet(1);
		resiResu(1) = equaFact(1,0) + equaFact(1,1) * slavXiet(0) + equaFact(1,2) * slavXiet(1) 
			+ equaFact(1,3) * slavXiet(0) * slavXiet(1);
		Eigen::Matrix<double,2,2> hessMatr;
		hessMatr(0,0) = equaFact(0,1) + equaFact(0,3) * slavXiet(1);
		hessMatr(0,1) = equaFact(0,2) + equaFact(0,3) * slavXiet(0);
		hessMatr(1,0) = equaFact(1,1) + equaFact(1,3) * slavXiet(1);
		hessMatr(1,1) = equaFact(1,2) + equaFact(1,3) * slavXiet(0);
		Eigen::Matrix<double,2,1> DELTxiet = - hessMatr.fullPivLu().solve(resiResu);
		if(sqrt(DELTxiet(0) * DELTxiet(0) + DELTxiet(1) * DELTxiet(1)) < 1.0E-14){
			if(sqrt(resiResu(0) * resiResu(0) + resiResu(1) * resiResu(1)) < 1.0E-15){
				break;
			}
		}
		slavXiet = slavXiet + DELTxiet;
	}
	if(tj >= 60){
		std::cout << "Error 1 in CONTACT::PROJECT_MTS!" << std::endl;
	}
	Eigen::Vector3d tang_1 = PrmaPxie.block(0,0,3,1);
	Eigen::Vector3d tang_2 = PrmaPxie.block(0,1,3,1);
	Eigen::Vector3d tempNorm = (tang_1.cross(tang_2)).normalized();
	Eigen::Matrix<double,4,1> tempVect;
	tempVect << 1.0, slavXiet(0), slavXiet(1), slavXiet(0) * slavXiet(1);
	Eigen::Vector3d slavPoin = factMatr * tempVect;
	ngap = tempNorm.dot(slavPoin - mastPoin);
	return 1;
}

long CSEARCH::PROJECT_STM_SUB(const std::vector<Eigen::Vector3d> &mastCorn, 
	const Eigen::Vector3d &slavPoin, 
	Eigen::Matrix<double,3,2> &PrmaPxie, Eigen::Matrix<double,2,1> &mastXiet, double &ngap){
	Eigen::Matrix<double,3,4> factMatr = Eigen::MatrixXd::Zero(3,4);
	for(long ti = 0; ti < 3; ti ++){
		for(long tj = 0; tj < 4; tj ++){
			factMatr(ti,0) += mastCorn[tj](ti) / 4.0;
			factMatr(ti,1) += mastCorn[tj](ti) * biliQuad.nacoCorn[tj][0] / 4.0;
			factMatr(ti,2) += mastCorn[tj](ti) * biliQuad.nacoCorn[tj][1] / 4.0;
			factMatr(ti,3) += mastCorn[tj](ti) 
				* biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] / 4.0;
		}
	}
	Eigen::Matrix<double,2,6> equaFact;
	equaFact(0,0) = factMatr(0,0) * factMatr(0,1) - factMatr(0,1) * slavPoin(0) 
		+ factMatr(1,0) * factMatr(1,1) - factMatr(1,1) * slavPoin(1) 
		+ factMatr(2,0) * factMatr(2,1) - factMatr(2,1) * slavPoin(2);
	equaFact(0,1) = factMatr(0,1) * factMatr(0,1) 
		+ factMatr(1,1) * factMatr(1,1) + factMatr(2,1) * factMatr(2,1);
	equaFact(0,2) = factMatr(0,1) * factMatr(0,2) +factMatr(0,0) * factMatr(0,3) 
		- factMatr(0,3) * slavPoin(0) 
		+ factMatr(1,1) * factMatr(1,2) +factMatr(1,0) * factMatr(1,3) 
		- factMatr(1,3) * slavPoin(1) 
		+ factMatr(2,1) * factMatr(2,2) +factMatr(2,0) * factMatr(2,3) 
		- factMatr(2,3) * slavPoin(2);
	equaFact(0,3) = 2.0 * factMatr(0,1) * factMatr(0,3) 
		+ 2.0 * factMatr(1,1) * factMatr(1,3) 
		+ 2.0 * factMatr(2,1) * factMatr(2,3);
	equaFact(0,4) = factMatr(0,2) * factMatr(0,3) 
		+ factMatr(1,2) * factMatr(1,3) + factMatr(2,2) * factMatr(2,3);
	equaFact(0,5) = factMatr(0,3) * factMatr(0,3) 
		+ factMatr(1,3) * factMatr(1,3) + factMatr(2,3) * factMatr(2,3);
	equaFact(1,0) = factMatr(0,0) * factMatr(0,2) - factMatr(0,2) * slavPoin(0) 
		+ factMatr(1,0) * factMatr(1,2) - factMatr(1,2) * slavPoin(1) 
		+ factMatr(2,0) * factMatr(2,2) - factMatr(2,2) * slavPoin(2);
	equaFact(1,1) = factMatr(0,2) * factMatr(0,2) 
		+ factMatr(1,2) * factMatr(1,2) + factMatr(2,2) * factMatr(2,2);
	equaFact(1,2) = equaFact(0,2);
	equaFact(1,3) = 2.0 * factMatr(0,2) * factMatr(0,3) 
		+ 2.0 * factMatr(1,2) * factMatr(1,3) 
		+ 2.0 * factMatr(2,2) * factMatr(2,3);
	equaFact(1,4) = factMatr(0,1) * factMatr(0,3) 
		+ factMatr(1,1) * factMatr(1,3) + factMatr(2,1) * factMatr(2,3);
	equaFact(1,5) = equaFact(0,5);
	// mastXiet = Eigen::MatrixXd::Zero(2,1);
	//
	long tj;
	for(tj = 0; tj < 60; tj ++){
		Eigen::Matrix<double,2,1> resiResu;
		resiResu(0) = equaFact(0,0) + equaFact(0,1) * mastXiet(0) + equaFact(0,2) * mastXiet(1) 
			+ equaFact(0,3) * mastXiet(0) * mastXiet(1) 
			+ equaFact(0,4) * mastXiet(1) * mastXiet(1) 
			+ equaFact(0,5) * mastXiet(0) * mastXiet(1) * mastXiet(1);
		resiResu(1) = equaFact(1,0) + equaFact(1,1) * mastXiet(1) + equaFact(1,2) * mastXiet(0) 
			+ equaFact(1,3) * mastXiet(0) * mastXiet(1) 
			+ equaFact(1,4) * mastXiet(0) * mastXiet(0) 
			+ equaFact(1,5) * mastXiet(0) * mastXiet(0) * mastXiet(1);
		Eigen::Matrix<double,2,2> hessMatr;
		hessMatr(0,0) = equaFact(0,1) + equaFact(0,3) * mastXiet(1) 
			+ equaFact(0,5) * mastXiet(1) * mastXiet(1);
		hessMatr(0,1) = equaFact(0,2) + equaFact(0,3) * mastXiet(0) 
			+ equaFact(0,4) * 2.0 * mastXiet(1) + equaFact(0,5) * mastXiet(0) * 2.0 * mastXiet(1);
		hessMatr(1,0) = equaFact(1,2) + equaFact(1,3) * mastXiet(1) 
			+ equaFact(1,4) * 2.0 * mastXiet(0) + equaFact(1,5) * 2.0 * mastXiet(0) * mastXiet(1);
		hessMatr(1,1) = equaFact(1,1) + equaFact(1,3) * mastXiet(0) 
			+ equaFact(1,5) * mastXiet(0) * mastXiet(0);
		Eigen::Matrix<double,2,1> DELTxiet = - hessMatr.fullPivLu().solve(resiResu);
		if(sqrt(DELTxiet(0) * DELTxiet(0) + DELTxiet(1) * DELTxiet(1)) < 1.0E-12){
			if(sqrt(resiResu(0) * resiResu(0) + resiResu(1) * resiResu(1)) < 1.0E-15){
				break;
			}
		}
		mastXiet = mastXiet + DELTxiet;
	}
	if(tj >= 60 && debuMode >= 1){
		std::cout << "Warning 1 in CONTACT::PROJECT_STM_SUB!" << std::endl;
	}
	for(long ti = 0; ti < 3; ti ++){
		PrmaPxie(ti,0) = factMatr(ti,1) + factMatr(ti,3) * mastXiet(1);
		PrmaPxie(ti,1) = factMatr(ti,2) + factMatr(ti,3) * mastXiet(0);
	}
	Eigen::Vector3d tang_1 = PrmaPxie.block(0,0,3,1);
	Eigen::Vector3d tang_2 = PrmaPxie.block(0,1,3,1);
	Eigen::Vector3d tempNorm = (tang_1.cross(tang_2)).normalized();
	Eigen::Matrix<double,4,1> tempVect;
	tempVect << 1.0, mastXiet(0), mastXiet(1), mastXiet(0) * mastXiet(1);
	Eigen::Vector3d mastPoin = factMatr * tempVect;
	ngap = tempNorm.dot(slavPoin - mastPoin);
	return 1;
}

//the slavePoin maybe very far away from mastCorn
long CSEARCH::PROJECT_STM(const std::vector<Eigen::Vector3d> &mastCorn, 
	const Eigen::Vector3d &slavPoin, 
	Eigen::Matrix<double,3,2> &PrmaPxie, Eigen::Matrix<double,2,1> &mastXiet, double &ngap){
	double elemDiam = 1.0E15;
	Eigen::Vector3d mastCent = Eigen::Vector3d::Zero(3);
	for(long ti = 0; ti < 4; ti ++){
		for(long tj = ti + 1; tj < 4; tj ++){
			elemDiam = std::min(elemDiam, (mastCorn[ti] - mastCorn[tj]).norm());
		}
		mastCent += 0.25 * mastCorn[ti];
	}
	mastXiet = Eigen::MatrixXd::Zero(2,1);
	if((slavPoin - mastCent).norm() <= elemDiam){
		PROJECT_STM_SUB(mastCorn, slavPoin, PrmaPxie, mastXiet, ngap);
	}
	else{
		double realTime = (slavPoin - mastCent).norm() / elemDiam;
		Eigen::Vector3d slavDire = (slavPoin - mastCent) / (slavPoin - mastCent).norm();
		Eigen::Matrix<double,3,2> PrmaPxie_1;
		double ngap_1;
		for(long ti = 0; ti < realTime; ti ++){
			Eigen::Vector3d slavPoin_1 = mastCent + ti * elemDiam * slavDire;
			PROJECT_STM_SUB(mastCorn, slavPoin_1, PrmaPxie_1, mastXiet, ngap_1);
		}
		PROJECT_STM_SUB(mastCorn, slavPoin, PrmaPxie, mastXiet, ngap);
	}
	return 1;
}

long CSEARCH::MAST_TO_SLAV(const std::vector<Eigen::Vector3d> &mastCorn, 
	const std::vector<Eigen::Vector3d> &slavCorn, const Eigen::Vector2d mastXiet, 
	VECTOR2D &maslShap, std::vector<Eigen::Vector3d> &maslPoin, 
	std::vector<Eigen::Vector3d> &basiVect, double &weigFact){
	Eigen::Matrix<double,3,2> PrmaPxie = Eigen::MatrixXd::Zero(3,2);
	Eigen::Matrix<double,2,1> slavXiet = Eigen::MatrixXd::Zero(2,1);
	double ngap;
	PROJECT_MTS(mastCorn, slavCorn, mastXiet, maslPoin[0], PrmaPxie, slavXiet, ngap);
	Eigen::Vector3d tang_1 = PrmaPxie.block(0,0,3,1);
	Eigen::Vector3d tang_2 = PrmaPxie.block(0,1,3,1);
	basiVect.resize(3);
	basiVect[0] = (tang_1.cross(tang_2)).normalized();
	basiVect[1] = tang_1.normalized();
	basiVect[2] = tang_2.normalized();
	weigFact = sqrt(pow(PrmaPxie(1,0) * PrmaPxie(2,1) - PrmaPxie(1,1) * PrmaPxie(2,0), 2.0) 
		+ pow(PrmaPxie(2,0) * PrmaPxie(0,1) - PrmaPxie(2,1) * PrmaPxie(0,0), 2.0) 
		+ pow(PrmaPxie(0,0) * PrmaPxie(1,1) - PrmaPxie(0,1) * PrmaPxie(1,0), 2.0));
	if(debuMode >= 1 && (slavXiet(0) < -1.0 - 1.0E-6 || slavXiet(0) > 1.0 + 1.0E-6 
		|| slavXiet(1) < -1.0 - 1.0E-6 || slavXiet(1) > 1.0 + 1.0E-6)){
		std::cout << "Warning 1 in CONTACT::MAST_TO_SLAV!" << std::endl;
	}
	maslPoin[1] = Eigen::Vector3d::Zero();
	for(long ti = 0; ti < 4; ti ++){
		maslShap[1].resize(4);
		maslShap[1][ti] = (1.0 + biliQuad.nacoCorn[ti][0] * slavXiet(0)) 
			* (1.0 + biliQuad.nacoCorn[ti][1] * slavXiet(1)) / 4.0;
		maslPoin[1] = maslPoin[1] + maslShap[1][ti] * slavCorn[ti];
	}
	return 1;
}

double CSEARCH::TRIANGLE_AREA_2D(Eigen::Vector2d tempPoin_0, 
	Eigen::Vector2d tempPoin_1, Eigen::Vector2d tempPoin_2){
	Eigen::Vector2d tempVect_0 = tempPoin_1 - tempPoin_0;
	Eigen::Vector2d tempVect_1 = tempPoin_2 - tempPoin_0;
	return abs(tempVect_0(0) * tempVect_1(1) - tempVect_0(1) * tempVect_1(0)) / 2.0;
}

long CSEARCH::TRIANGLE_QUADRATURE(Eigen::Vector2d tempXiet_0, 
	Eigen::Vector2d tempXiet_1, Eigen::Vector2d tempXiet_2, 
	std::vector<Eigen::Vector2d> &listXiet, std::vector<double> &listWeig){
	double area = TRIANGLE_AREA_2D(tempXiet_0, tempXiet_1, tempXiet_2);
	for(long ti = 0; ti < triaQuad.numbNgip; ti ++){
		Eigen::Vector2d tempXiet = triaQuad.nacoNgip[ti][0] * tempXiet_0 
			+ triaQuad.nacoNgip[ti][1] * tempXiet_1 
			+ triaQuad.nacoNgip[ti][2] * tempXiet_2;
		listXiet.emplace_back(tempXiet);
		//integral variable transformation: from Cartesian coordinate to shape function
		//x:xmin~xmax, y:ymin(x)~ymax(x)
		//A1:0~1, A2:0~1-A1, A3=1-A0-A1
		listWeig.emplace_back(2.0 * area * triaQuad.niwfNgip[ti]);
	}
	return 1;
}

long CSEARCH::SORT_BY_2D(Eigen::Vector2d &tempPoin_0, Eigen::Vector2d &tempPoin_1, long tempIden){
	if(tempPoin_0(tempIden) > tempPoin_1(tempIden)){
		Eigen::Vector2d tempPoin;
		tempPoin = tempPoin_0;
		tempPoin_0 = tempPoin_1;
		tempPoin_1 = tempPoin;
	}
	return 1;
}

long CSEARCH::IS_CROSS_2D(Eigen::Vector2d tempPoin_0, Eigen::Vector2d tempPoin_1, 
		Eigen::Vector2d tempPoin_2, Eigen::Vector2d tempPoin_3){
	if(std::max(tempPoin_0(0), tempPoin_1(0)) < std::min(tempPoin_2(0), tempPoin_3(0))
		|| std::max(tempPoin_0(1), tempPoin_1(1)) < std::min(tempPoin_2(1), tempPoin_3(1))
		|| std::min(tempPoin_0(0), tempPoin_1(0)) > std::max(tempPoin_2(0), tempPoin_3(0))
		|| std::min(tempPoin_0(1), tempPoin_1(1)) > std::max(tempPoin_2(1), tempPoin_3(1))){
		return -1;
	}
	if(((tempPoin_2(0) - tempPoin_0(0)) * (tempPoin_2(1) - tempPoin_3(1)) 
		- (tempPoin_2(1) - tempPoin_0(1)) * (tempPoin_2(0) - tempPoin_3(0))) 
		* ((tempPoin_2(0) - tempPoin_1(0)) * (tempPoin_2(1) - tempPoin_3(1)) 
		- (tempPoin_2(1) - tempPoin_1(1)) * (tempPoin_2(0) - tempPoin_3(0))) <= 0
		&& ((tempPoin_0(0) - tempPoin_2(0)) * (tempPoin_0(1) - tempPoin_1(1)) 
		- (tempPoin_0(1) - tempPoin_2(1)) * (tempPoin_0(0) - tempPoin_1(0))) 
		* ((tempPoin_0(0) - tempPoin_3(0)) * (tempPoin_0(1) - tempPoin_1(1)) 
		- (tempPoin_0(1) - tempPoin_3(1)) * (tempPoin_0(0) - tempPoin_1(0))) <= 0){
		return 1;
	}
	else{
		return -1;
	}
}

long CSEARCH::LINE_INTERSECT_2D(Eigen::Vector2d tempPoin_0, Eigen::Vector2d tempPoin_1, 
	Eigen::Vector2d tempPoin_2, Eigen::Vector2d tempPoin_3, 
	std::vector<Eigen::Vector2d> &resuInte){
	if(IS_CROSS_2D(tempPoin_0, tempPoin_1, tempPoin_2, tempPoin_3) == -1){
		return 1;
	}
	double area_2 = TRIANGLE_AREA_2D(tempPoin_2, tempPoin_0, tempPoin_1);
	double area_3 = TRIANGLE_AREA_2D(tempPoin_3, tempPoin_0, tempPoin_1);
	double miniArea = 1.0E-12;// * (tempPoin_0 - tempPoin_1).norm();//?????????
	if(abs(area_2) < miniArea && abs(area_3) < miniArea){//collinear
		if(abs(tempPoin_0(0) - tempPoin_1(0)) > 1.0E-10){
			SORT_BY_2D(tempPoin_0, tempPoin_1, 0);
			SORT_BY_2D(tempPoin_2, tempPoin_3, 0);
			double from_x = tempPoin_0(0);
			double from_y = tempPoin_0(1);
			if(tempPoin_0(0) < tempPoin_2(0)){
				from_x = tempPoin_2(0);
				from_y = tempPoin_2(1);
			}
			double to_x = tempPoin_1(0);
			double to_y = tempPoin_1(1);
			if(tempPoin_1(0) > tempPoin_3(0)){
				to_x = tempPoin_3(0);
				to_y = tempPoin_3(1);
			}
			Eigen::Vector2d tempPoin;
			if(abs(from_x - to_x) < 1.0E-10){
				tempPoin << from_x, from_y;
				resuInte.emplace_back(tempPoin);
			}
			else{
				tempPoin << from_x, from_y;
				resuInte.emplace_back(tempPoin);
				tempPoin << to_x, to_y;
				resuInte.emplace_back(tempPoin);
			}
		}
		else{
			SORT_BY_2D(tempPoin_0, tempPoin_1, 1);
			SORT_BY_2D(tempPoin_2, tempPoin_3, 1);
			double from_x = tempPoin_0(0);
			double from_y = tempPoin_0(1);
			if(tempPoin_0(1) < tempPoin_2(1)){
				from_x = tempPoin_2(0);
				from_y = tempPoin_2(1);
			}
			double to_x = tempPoin_1(0);
			double to_y = tempPoin_1(1);
			if(tempPoin_1(1) > tempPoin_3(1)){
				to_x = tempPoin_3(0);
				to_y = tempPoin_3(1);
			}
			Eigen::Vector2d tempPoin;
			if(abs(from_x - to_x) < 1.0E-10){
				tempPoin << from_x, from_y;
				resuInte.emplace_back(tempPoin);
			}
			else{
				tempPoin << from_x, from_y;
				resuInte.emplace_back(tempPoin);
				tempPoin << to_x, to_y;
				resuInte.emplace_back(tempPoin);
			}
		}
	}
	else if(abs(area_2) < miniArea){//one endpoint lies on the another line-segment
		resuInte.emplace_back(tempPoin_2);
	}
	else if(abs(area_3) < miniArea){//one endpoint lies on the another line-segment
		resuInte.emplace_back(tempPoin_3);
	}
	else{//true intersect
		double tempFact = area_2 / area_3;
		Eigen::Vector2d tempPoin;
		tempPoin << (tempPoin_2(0) + tempFact * tempPoin_3(0)) / (1.0 + tempFact), 
			(tempPoin_2(1) + tempFact * tempPoin_3(1)) / (1.0 + tempFact);
		resuInte.emplace_back(tempPoin);
	}
	return 1;
}

long CSEARCH::IN_CQUAD_2D(Eigen::Vector2d tempPoin, const std::vector<Eigen::Vector2d> tempCorn){
	double subqArea = 0.0;
	for(long ti = 0; ti < 4; ti ++){
		subqArea += TRIANGLE_AREA_2D(tempPoin, tempCorn[ti], tempCorn[(ti + 1) % 4]);
	}
	double areaTota = TRIANGLE_AREA_2D(tempCorn[0], tempCorn[1], tempCorn[2]) + 
		TRIANGLE_AREA_2D(tempCorn[2], tempCorn[3], tempCorn[0]);
	if(subqArea <= (1.0 + 1.0E-12) * areaTota){
		return 1;
	}
	else{
		return -1;
	}
}

long CSEARCH::SI_SUB(const std::vector<Eigen::Vector3d> &mastCorn, 
	const std::vector<Eigen::Vector3d> &slavCorn, 
	std::vector<Eigen::Vector2d> &listXiet, std::vector<double> &listWeig){
	//projection
	std::vector<Eigen::Vector2d> mastProj(4), slavProj(4);
	mastProj[0] << -1.0, -1.0;
	mastProj[1] << 1.0, - 1.0;
	mastProj[2] << 1.0, 1.0;
	mastProj[3] << -1.0, 1.0;
	for(long ti = 0; ti < 4; ti ++){
		Eigen::Matrix<double,3,2> PrmaPxie;
		double ngap;
		PROJECT_STM(mastCorn, slavCorn[ti], PrmaPxie, slavProj[ti], ngap);
	}
	//intersection
	std::vector<Eigen::Vector2d> tempInte_0, tempInte_1, tempInte_2;
	for(long ti = 0; ti < 4; ti ++){
		if(IN_CQUAD_2D(slavProj[ti], mastProj) == 1){
			tempInte_0.push_back(slavProj[ti]);
		}
		if(IN_CQUAD_2D(mastProj[ti], slavProj) == 1){
			tempInte_0.push_back(mastProj[ti]);
		}
	}
	for(long ti = 0; ti < 4; ti ++){
		for(long tj = 0; tj < 4; tj ++){
			LINE_INTERSECT_2D(mastProj[ti], mastProj[(ti + 1) % 4], 
				slavProj[tj], slavProj[(tj + 1) % 4], tempInte_0
			);
		}
	}
	if(tempInte_0.size() < 3){
		return -1;
	}
	//no repeat
	Eigen::VectorXi tempInde = Eigen::VectorXi::LinSpaced(
		tempInte_0.size(), 0, tempInte_0.size() - 1
	);
	auto rule = [tempInte_0](int ti, int tj)->bool{
		if(tempInte_0[ti](0) < tempInte_0[tj](0) - 1.0E-10){
			return true;
		}
		else if(tempInte_0[ti](0) <= tempInte_0[tj](0) + 1.0E-10){
			if(tempInte_0[ti](1) < tempInte_0[tj](1) - 1.0E-10){
				return true;
			}
			else if(tempInte_0[ti](1) <= tempInte_0[tj](1) + 1.0E-10){
				return false;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	};
	std::sort(tempInde.data(), tempInde.data() + tempInde.size(), rule);
	tempInte_1.emplace_back(tempInte_0[tempInde(0)]);
	for(long ti = 1; ti < tempInte_0.size(); ti ++){
		if(abs(tempInte_0[tempInde(ti)](0) - tempInte_0[tempInde(ti-1)](0)) > 1.0E-10
			|| abs(tempInte_0[tempInde(ti)](1) - tempInte_0[tempInde(ti-1)](1)) > 1.0E-10){
			tempInte_1.emplace_back(tempInte_0[tempInde(ti)]);
		}
	}
	//sort
	Eigen::Vector2d tempCent = Eigen::Vector2d::Zero(2);
	for(const auto& iterTein : tempInte_1){
		tempCent = tempCent + iterTein;
	}
	tempCent = tempCent / tempInte_1.size();
	Eigen::VectorXd tempAngl;
	tempAngl.resize(tempInte_1.size());
	long tempNumb = 0;
	for(const auto& iterTein : tempInte_1){
		Eigen::Vector2d tempVect = iterTein - tempCent;
		tempAngl(tempNumb) = atan2(tempVect(1), tempVect(0));
		tempNumb ++;
	}
	tempInde = Eigen::VectorXi::LinSpaced(tempInte_1.size(), 0, tempInte_1.size() - 1);
	auto rule_1 = [tempAngl](int ti, int tj)->bool{
		if(tempAngl(ti) < tempAngl(tj)){
			return true;
		}
		else{
			return false;
		}
	};
	std::sort(tempInde.data(), tempInde.data() + tempInde.size(), rule_1);
	for(long ti = 0; ti < tempInte_1.size(); ti ++){
		tempInte_2.emplace_back(tempInte_1[tempInde(ti)]);
	}
	//
	//centroid and area
	//R. Nurnberg. Calculating the area and centroid of a polygon in 2d.
	//https://paulbourke.net/geometry/polygonmesh/centroid.pdf
	double area = 0.0;
	tempCent = Eigen::VectorXd::Zero(2);
	long tempSize = tempInte_2.size();
	for(long ti = 0; ti < tempInte_2.size(); ti ++){
		area += tempInte_2[ti](0) * tempInte_2[(ti + 1) % tempSize](1) 
			- tempInte_2[(ti + 1) % tempSize](0) * tempInte_2[ti](1);
		tempCent(0) += (tempInte_2[ti](0) + tempInte_2[(ti + 1) % tempSize](0)) 
			* (tempInte_2[ti](0) * tempInte_2[(ti + 1) % tempSize](1) 
			- tempInte_2[(ti + 1) % tempSize](0) * tempInte_2[ti](1));
		tempCent(1) += (tempInte_2[ti](1) + tempInte_2[(ti + 1) % tempSize](1)) 
			* (tempInte_2[ti](0) * tempInte_2[(ti + 1) % tempSize](1) 
			- tempInte_2[(ti + 1) % tempSize](0) * tempInte_2[ti](1));
	}
	area /= 2.0;
	if(abs(area) <= miniArea){
		return -1;
	}
	tempCent = tempCent / 6.0 / area;
	for(long ti = 0; ti < tempInte_2.size(); ti ++){
		TRIANGLE_QUADRATURE(tempCent, tempInte_2[ti], 
			tempInte_2[(ti + 1) % tempInte_2.size()], listXiet, listWeig);
	}
	return 1;
}

long CSEARCH::SEGMENT_INTERSECT(long tempMast, long tempSlav, 
	std::vector<INTEGRAL_POINT> &inpuInpo){
	std::vector<Eigen::Vector3d> mastCorn(4), slavCorn(4);
	for(long ti = 0; ti < 4; ti ++){
		auto iterNoco = (*mastGrid).nodeCoor.find(mastSegm[tempMast][ti]);
		mastCorn[ti] << (iterNoco->second)[0], (iterNoco->second)[1], (iterNoco->second)[2];
		iterNoco = (*slavGrid).nodeCoor.find(slavSegm[tempSlav][ti]);
		slavCorn[ti] << (iterNoco->second)[0], (iterNoco->second)[1], (iterNoco->second)[2];
	}
	std::vector<Eigen::Vector2d> listXiet;
	std::vector<double> listWeig;
	SI_SUB(mastCorn, slavCorn, listXiet, listWeig);
	//shape function, initial gap
	for(long ti = 0; ti < listXiet.size(); ti ++){
		//
		VECTOR2D maslShap(2, std::vector<double>(4));
		std::vector<Eigen::Vector3d> maslPoin(2, Eigen::Vector3d::Zero());
		for(long tj = 0; tj < 4; tj ++){
			maslShap[0][tj] = (1.0 + biliQuad.nacoCorn[tj][0] * listXiet[ti](0)) 
				* (1.0 + biliQuad.nacoCorn[tj][1] * listXiet[ti](1)) / 4.0;
			maslPoin[0] = maslPoin[0] + maslShap[0][tj] * mastCorn[tj];
		}
		std::vector<Eigen::Vector3d> basiVect(3);
		double weigFact;
		MAST_TO_SLAV(mastCorn, slavCorn, listXiet[ti], maslShap, maslPoin, basiVect, weigFact);
		Eigen::Vector3d mastSlav = maslPoin[1] - maslPoin[0];
		//
		INTEGRAL_POINT tempInpo;
		for(long tj = 0; tj < 4; tj ++){
			tempInpo.node[0][tj] = mastSegm[tempMast][tj];
			tempInpo.node[1][tj] = slavSegm[tempSlav][tj];
		}
		tempInpo.shapFunc = maslShap;
		tempInpo.contPoin = maslPoin;
		tempInpo.basiVect = basiVect;
		tempInpo.initNgap = mastSlav.dot(basiVect[0]);
		tempInpo.quadWeig = listWeig[ti] * weigFact;
		inpuInpo.emplace_back(tempInpo);
	}
	return 1;
}

long CSEARCH::CONTACT_SEARCH(const VECTOR2D &slavCoor, double maxiDist){
	for(long ti = 0; ti < slavSegm.size(); ti ++){
		if(ti%1000 == 0){
			std::cout << "The " << ti << "/" << slavSegm.size() << "-th segment." << std::endl;
		}
		long buck_r = (slavCoor[0][ti] - buckLoco[0][0]) / buckLoco[0][2];
		long buck_c = (slavCoor[1][ti] - buckLoco[1][0]) / buckLoco[1][2];
		if(buck_r < 0 || buck_r > buck.size() - 1 || buck_c < 0 || buck_c > buck[0].size() - 1){
			continue;
		}
		for(long tj = std::max(buck_r - 1, (long)0); 
			tj <= std::min(buck_r + 1, (long)(buck.size() - 1)); tj ++){
			for(long tk = std::max(buck_c - 1, (long)0); 
				tk <= std::min(buck_c + 1, (long)(buck[0].size() - 1)); tk ++){
				for(const auto &iter_0 : buck[tj][tk]){
					std::vector<INTEGRAL_POINT> intePoin_e;
					SEGMENT_INTERSECT(iter_0, ti, intePoin_e);
					bool tempFlag = false;
					for(const auto &iterInpo : intePoin_e){
						if(iterInpo.initNgap <= maxiDist){
							tempFlag = true;
							break;
						}
					}
					if(tempFlag == true){
						intePoin.insert(intePoin.end(), intePoin_e.begin(), intePoin_e.end());
					}
				}
			}
		}
	}
	double miniNgap = 1.0E12;
	double maxiNgap = -1.0E12;
	for(long ti = 0; ti < intePoin.size(); ti ++){
		miniNgap = std::min(miniNgap, intePoin[ti].initNgap);
		maxiNgap = std::max(maxiNgap, intePoin[ti].initNgap);
	}
	std::cout << "Initial normal gap: minimum = " << miniNgap 
		<< ", maximum = " << maxiNgap << std::endl;
	return 1;
}

long CSEARCH::OUTPUT_INPO(long contIden){
	std::stringstream tempStre;
	tempStre << "resuInpo_" << contIden << ".txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long ti = 0; ti < intePoin.size(); ti ++){
		tempOfst << std::setw(30) << intePoin[ti].contPoin[0](0) 
			<< std::setw(30) << intePoin[ti].contPoin[0](1) 
			<< std::setw(30) << intePoin[ti].contPoin[0](2) 
			<< std::setw(30) << intePoin[ti].contPoin[1](0) 
			<< std::setw(30) << intePoin[ti].contPoin[1](1) 
			<< std::setw(30) << intePoin[ti].contPoin[1](2) 
			<< std::setw(30) << intePoin[ti].initNgap << std::endl;
	}
	tempOfst.close();
	return 1;
}

long CSEARCH::ADAPTIVE_REFINE(MULTIGRID *tempMast, MULTIGRID *tempSlav, bool &isnoRefi, 
	CURVEDS *mastSurf, CURVEDS *slavSurf, long tempLeve, 
	double distCrit, std::vector<long> buckNumb, 
	std::function<void(COOR, double &, double &)> CART_CURV){
	mastGrid = tempMast;
	slavGrid = tempSlav;
	EFACE_SURFACE iterEfsu_0(tempMast, mastSurf);
	while(iterEfsu_0.INCREMENT() == 1){
		//NO NEED: if((* tempMast).elemVect[iterEfsu_0.eid].level == tempLeve){
		mastSegm.emplace_back(iterEfsu_0.currNode);
	}
	EFACE_SURFACE iterEfsu_1(tempSlav, slavSurf);
	while(iterEfsu_1.INCREMENT() == 1){
		slavSegm.emplace_back(iterEfsu_1.currNode);
	}
	if(debuMode >= 1){
		OUTPUT_COSE(0);
	}
	//
	VECTOR2D mastCoor(2);
	for(long ti = 0; ti < mastSegm.size(); ti ++){
		double tempXico = 0.0;
		double tempEtac = 0.0;
		for(long tj = 0; tj < 4; tj ++){
			long node_tj = mastSegm[ti][tj];
			auto iterNoco = (*tempMast).nodeCoor.find(node_tj);
			CART_CURV(iterNoco->second, tempXico, tempEtac);
		}
		mastCoor[0].emplace_back(tempXico / 4.0);
		mastCoor[1].emplace_back(tempEtac / 4.0);
	}
	BUCKET_SORT(mastCoor, buckNumb);
	VECTOR2D slavCoor(2);
	for(long ti = 0; ti < slavSegm.size(); ti ++){
		double tempXico = 0.0;
		double tempEtac = 0.0;
		for(long tj = 0; tj < 4; tj ++){
			long node_tj = slavSegm[ti][tj];
			auto iterNoco = (*tempSlav).nodeCoor.find(node_tj);
			CART_CURV(iterNoco->second, tempXico, tempEtac);
		}
		slavCoor[0].emplace_back(tempXico / 4.0);
		slavCoor[1].emplace_back(tempEtac / 4.0);
	}
	//
	std::vector<std::set<long>> spliNode(2);
	for(long ti = 0; ti < slavSegm.size(); ti ++){
		if(ti%1000 == 0){
			std::cout << "The " << ti << "/" << slavSegm.size() << "-th segment." << std::endl;
		}
		long buck_r = (slavCoor[0][ti] - buckLoco[0][0]) / buckLoco[0][2];
		long buck_c = (slavCoor[1][ti] - buckLoco[1][0]) / buckLoco[1][2];
		if(buck_r < 0 || buck_r > buck.size() - 1 || buck_c < 0 || buck_c > buck[0].size() - 1){
			continue;
		}
		for(long tj = std::max(buck_r - 1, (long)0); 
			tj <= std::min(buck_r + 1, (long)(buck.size() - 1)); tj ++){
			for(long tk = std::max(buck_c - 1, (long)0); 
				tk <= std::min(buck_c + 1, (long)(buck[0].size() - 1)); tk ++){
				for(const auto &iterBuck : buck[tj][tk]){
					intePoin.clear();
					SEGMENT_INTERSECT(iterBuck, ti, intePoin);
					double miniNgap = 1.0E12;
					for(const auto &iterInpo : intePoin){
						miniNgap = std::min(miniNgap, iterInpo.initNgap);
					}
					if(miniNgap <= distCrit){
						for(const auto &iterInpo : intePoin){
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
		isnoRefi = false;
		return 1;
	}
	else{
		isnoRefi = true;
	}
	//
	for(long tv = 0; tv < 2; tv ++){
		MULTIGRID *tempGrid = (tv == 0) ? tempMast : tempSlav;
		CURVEDS *tempSurf = (tv == 0) ? mastSurf : slavSurf;
		//
		std::set<long> spliElem;
		std::map<long,std::set<long>> spliFlag;
		std::map<std::vector<long>,COOR> planSurf;
		spliElem.clear();
		for(long ti = 0; ti < (*tempGrid).elemVect.size(); ti ++){
			const auto &tempElem = (*tempGrid).elemVect[ti];
			if(tempElem.children.size() > 0 || tempElem.level != tempLeve){
				continue;
			}
			bool tempFlag = false;
			for(long tj = 0; tj < 8; tj ++){
				auto iterSpno = spliNode[tv].find(tempElem.cornNode[tj]);
				if(iterSpno != spliNode[tv].end()){
					tempFlag = true;
					break;
				}
			}
			if(tempFlag == true){
				spliElem.insert(ti);
				(*tempGrid).elemVect[ti].refiPatt = 0;
			}
		}
		planSurf.clear();
		(*tempSurf).REFINE(spliElem, (*tempGrid).elemVect, (*tempGrid).nodeCoor, planSurf);
		(*tempGrid).REFINE(spliElem, spliFlag, planSurf);
	}
	return 1;
}

#endif
