
#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include "MGPIS.h"

//map: ordered vector to long
typedef std::map<std::vector<long>,long> MOVL;

class MULTIGRID{
public:
	/*********************************************************************************************/
	//from node numbering to coordinate
	std::map<long,COOR> nodeCoor;
	//from coordinate to node numbering
	std::map<COOR,long> coorNode;
	//line has been used by which elements
	std::map<std::vector<long>,std::set<long>> lineUsed;
	//face has been used by which elements
	std::map<std::vector<long>,std::set<long>> faceUsed;
	//
	std::vector<TREE_ELEM> elemVect;
	//constructor
	MULTIGRID();
	long COPY(const MULTIGRID &targGrid);
	//an example to show the refinement scheme
	long SCHEME(long scheSett);
	long SCHEME_0();
	long SCHEME_1();
	//try to add one node, the added node may already exist
	long TRY_ADD_NODE(COOR tempCoor);
	//add one element
	long ADD_ELEMENT(TREE_ELEM elem);
	//refine the elements in spliElem
	//spliFlag: contains specified sub elements after refinement
	//planSurf: Cartesian - curvilinear - Cartesian
	long REFINE(std::set<long> &spliElem, 
		const std::map<long,std::set<long>> &spliFlag, 
		const std::map<std::vector<long>,COOR> &planSurf
	);
	//gradual level check
	long GRLE_CHECK(std::set<long> &spliElem);
	//output element information
	long OUTPUT_ELEMENT(long fileIden);
	//rotation and translation of rigid body
	long RIGI_ROTR(Eigen::Matrix3d rotaMatr, Eigen::Vector3d tranVect);
	/*********************************************************************************************/
	//nodes of different level: 0 - coarsest
	VECTOR2L leveNode;
	//node is locating at level/position
	VECTOR2L nodeLepo;
	std::vector<long> posiNode;//from position to node
	std::set<long> coupNode;//coupling nodes
	long coupReps;//coupling nodes representation
	//prolongation operator
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> scalProl;
	std::map<long,std::vector<long>> finoCono;//finer node to coarser node
	std::map<std::vector<long>,long> conoFino;//coarser node to finer node
	//local refinement: pass the patch test
	long PATCH(const MOVL &tempMovl);
	//transfer operator
	long TRANSFER();
	//
	double mateElas;//material elasticity
	double matePois;//material Poisson ratio
	Eigen::Matrix<double,6,6> elatMatr;//elastic matrix
	//earliest transformation, nodeCoor : nodeLepo
	Eigen::SparseMatrix<double,Eigen::RowMajor> earlTran;
	Eigen::SparseMatrix<double,Eigen::RowMajor> scalEarl;
	long STIF_MATR();//construct original stiffness matrix
	/*********************************************************************************************/
	//prolongation operator
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> prolOper;
	//original stiffness matrix: without constraint
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> origStif;
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> consOper;
	std::map<long,double> consDofv;//DOF's displacement is given
	std::map<long,double> exteForc;//DOF's external force is given
	Eigen::VectorXi consFlag;//0 - constrained, 1 - not constrained
	Eigen::VectorXd dispForc;//value of constrained DOF
	Eigen::VectorXd consForc;//force vector after displacement constraint
	long LOAD_ACCU(long free_ti, double inpuLoad);//accumulate inpuLoad to DOF free_ti
	std::map<long,Eigen::Matrix3d> nodeRota;//node rotation matrix
	//calculate origStif, mgpi.consStif considering consDofv, exteForc
	long CONSTRAINT(long geomMult);
	long ADDITIONAL_FORCE(Eigen::VectorXd &addiForc);
	/*********************************************************************************************/
	long OUTP_SUB1(const Eigen::VectorXd &resuDisp, Eigen::VectorXd &outpDisp);
	long OUTP_SUB2(const Eigen::VectorXd &outpDisp, long fileIden);
	long OUTP_SUB2(const Eigen::VectorXd &outpDisp, std::string fileIden);
	Eigen::VectorXd OUTP_DISP(Eigen::VectorXd resuDisp, long fileIden);//output displacement
	long STRESS_RECOVERY(Eigen::VectorXd outpDisp, long fileIden);
	MGPIS mgpi;//multigrid preconditioned solver
};

MULTIGRID::MULTIGRID(){
	//
	mateElas = 210.0E9;
	matePois = 0.3;
	//
	mgpi.maxiLeve = -1;
}

long MULTIGRID::COPY(const MULTIGRID &targGrid){
	nodeCoor = targGrid.nodeCoor;
	coorNode = targGrid.coorNode;
	lineUsed = targGrid.lineUsed;
	faceUsed = targGrid.faceUsed;
	elemVect = targGrid.elemVect;
	leveNode = targGrid.leveNode;
	nodeLepo = targGrid.nodeLepo;
	posiNode = targGrid.posiNode;
	coupNode = targGrid.coupNode;
	coupReps = targGrid.coupReps;
	scalProl = targGrid.scalProl;
	finoCono = targGrid.finoCono;
	conoFino = targGrid.conoFino;
	mateElas = targGrid.mateElas;
	matePois = targGrid.matePois;
	elatMatr = targGrid.elatMatr;
	earlTran = targGrid.earlTran;
	scalEarl = targGrid.scalEarl;
	prolOper = targGrid.prolOper;
	origStif = targGrid.origStif;
	consOper = targGrid.consOper;
	consDofv = targGrid.consDofv;
	exteForc = targGrid.exteForc;
	consFlag = targGrid.consFlag;
	dispForc = targGrid.dispForc;
	consForc = targGrid.consForc;
	nodeRota = targGrid.nodeRota;
	mgpi.maxiLeve = targGrid.mgpi.maxiLeve;
	mgpi.realProl = targGrid.mgpi.realProl;
	mgpi.consStif = targGrid.mgpi.consStif;
	mgpi.consLowe = targGrid.mgpi.consLowe;
	mgpi.consDiag = targGrid.mgpi.consDiag;
	mgpi.consUppe = targGrid.mgpi.consUppe;
	return 1;
}

long MULTIGRID::SCHEME(long scheSett){
	//
	mkdir("RScheme", 0755);
	outpDire = "RScheme/";
	if(scheSett == 0){
		SCHEME_0();
	}
	else if(scheSett == 1){
		SCHEME_1();
	}
	return 1;
}

long MULTIGRID::SCHEME_0(){
	//
	double leng = 1.0;
	long diviNumb = 2;
	//nodes
	long tempNode[diviNumb + 1][diviNumb + 1][diviNumb + 1];
	for(long ti = 0; ti <= diviNumb; ti ++){
		for(long tj = 0; tj <= diviNumb; tj ++){
			for(long tk = 0; tk <= diviNumb; tk ++){
				double tempX = - leng / 2.0 + leng / diviNumb * (double)ti;
				double tempY = - leng / 2.0 + leng / diviNumb * (double)tj;
				double tempZ = 0.0 + leng / diviNumb * (double)tk;
				COOR tempCoor(tempX, tempY, tempZ);
				tempNode[ti][tj][tk] = TRY_ADD_NODE(tempCoor);
			}
		}
	}
	//elements
	for(long ti = 0; ti < diviNumb; ti ++){
		for(long tj = 0; tj < diviNumb; tj ++){
			for(long tk = 0; tk < diviNumb; tk ++){
				TREE_ELEM tempElem;
				tempElem.parent = -1;
				tempElem.cornNode[0] = tempNode[tk][tj][ti];
				tempElem.cornNode[1] = tempNode[tk + 1][tj][ti];
				tempElem.cornNode[2] = tempNode[tk + 1][tj + 1][ti];
				tempElem.cornNode[3] = tempNode[tk][tj + 1][ti];
				tempElem.cornNode[4] = tempNode[tk][tj][ti + 1];
				tempElem.cornNode[5] = tempNode[tk + 1][tj][ti + 1];
				tempElem.cornNode[6] = tempNode[tk + 1][tj + 1][ti + 1];
				tempElem.cornNode[7] = tempNode[tk][tj + 1][ti + 1];
				tempElem.level = 0;
				tempElem.refiPatt = 7;
				tempElem.children.resize(0);
				long elemNumb = ADD_ELEMENT(tempElem);
			}
		}
	}
	//refinement level 0
	std::set<long> spliElem;
	std::map<long,std::set<long>> spliFlag;
	std::map<std::vector<long>,COOR> planSurf;
	spliElem.clear();//
	spliElem.insert(6);
	elemVect[6].refiPatt = 0;
	spliElem.insert(7);
	elemVect[7].refiPatt = 0;
	spliFlag.clear();//
	spliFlag[6] = {5, 7};
	spliFlag[7] = {4, 6};
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//refinement level 1
	for(auto const &iterSpel : spliElem){
		elemVect[iterSpel].refiPatt = 0;
	}
	spliFlag.clear();//
	spliFlag[13] = {5, 7};
	spliFlag[15] = {5, 7};
	spliFlag[20] = {4, 6};
	spliFlag[22] = {4, 6};
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//refinement level 2
	for(auto const &iterSpel : spliElem){
		elemVect[iterSpel].refiPatt = 0;
	}
	spliFlag.clear();//
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//
	OUTPUT_ELEMENT(0);
	TRANSFER();
	return 1;
}

long MULTIGRID::SCHEME_1(){
	//
	double leng = 1.0;
	long diviNumb = 1;
	//nodes
	long tempNode[diviNumb + 1][diviNumb + 1][diviNumb + 1];
	for(long ti = 0; ti <= diviNumb; ti ++){
		for(long tj = 0; tj <= diviNumb; tj ++){
			for(long tk = 0; tk <= diviNumb; tk ++){
				double tempX = - leng / 2.0 + leng / diviNumb * (double)ti;
				double tempY = - leng / 2.0 + leng / diviNumb * (double)tj;
				double tempZ = 0.0 + leng / diviNumb * (double)tk;
				COOR tempCoor(tempX, tempY, tempZ);
				tempNode[ti][tj][tk] = TRY_ADD_NODE(tempCoor);
			}
		}
	}
	//elements
	for(long ti = 0; ti < diviNumb; ti ++){
		for(long tj = 0; tj < diviNumb; tj ++){
			for(long tk = 0; tk < diviNumb; tk ++){
				TREE_ELEM tempElem;
				tempElem.parent = -1;
				tempElem.cornNode[0] = tempNode[tk][tj][ti];
				tempElem.cornNode[1] = tempNode[tk + 1][tj][ti];
				tempElem.cornNode[2] = tempNode[tk + 1][tj + 1][ti];
				tempElem.cornNode[3] = tempNode[tk][tj + 1][ti];
				tempElem.cornNode[4] = tempNode[tk][tj][ti + 1];
				tempElem.cornNode[5] = tempNode[tk + 1][tj][ti + 1];
				tempElem.cornNode[6] = tempNode[tk + 1][tj + 1][ti + 1];
				tempElem.cornNode[7] = tempNode[tk][tj + 1][ti + 1];
				tempElem.level = 0;
				tempElem.refiPatt = 7;
				tempElem.children.resize(0);
				long elemNumb = ADD_ELEMENT(tempElem);
			}
		}
	}
	//refinement level 0
	std::set<long> spliElem;
	std::map<long,std::set<long>> spliFlag;
	std::map<std::vector<long>,COOR> planSurf;
	spliElem.clear();//
	spliElem.insert(0);
	elemVect[0].refiPatt = 6;
	spliFlag.clear();//
	spliFlag[0] = {0, 1};
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//refinement level 1
	spliElem.clear();//
	for(long ti = 0; ti < elemVect.size(); ti ++){
		if(elemVect[ti].children.size() > 0){
			continue;
		}
		spliElem.insert(ti);
		elemVect[ti].refiPatt = 3;
	}
	spliFlag.clear();//
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//refinement level 2
	spliElem.clear();//
	for(long ti = 0; ti < elemVect.size(); ti ++){
		if(elemVect[ti].children.size() > 0){
			continue;
		}
		spliElem.insert(ti);
		elemVect[ti].refiPatt = 2;
	}
	spliFlag.clear();//
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//refinement level 3
	spliElem.clear();//
	for(long ti = 0; ti < elemVect.size(); ti ++){
		if(elemVect[ti].children.size() > 0){
			continue;
		}
		spliElem.insert(ti);
		elemVect[ti].refiPatt = 0;
	}
	spliFlag.clear();//
	planSurf.clear();//
	REFINE(spliElem, spliFlag, planSurf);
	//
	OUTPUT_ELEMENT(0);
	TRANSFER();
	return 1;
}

long MULTIGRID::TRY_ADD_NODE(COOR tempCoor){
	auto iterCono = coorNode.find(tempCoor);
	if(iterCono == coorNode.end()){
		long tempSize = nodeCoor.size();
		nodeCoor.emplace(tempSize, tempCoor);
		coorNode.emplace(tempCoor, tempSize);
		return tempSize;
	}
	else{
		return iterCono->second;
	}
}

long MULTIGRID::ADD_ELEMENT(TREE_ELEM tempElem){
	elemVect.push_back(tempElem);
	long elemNumb = elemVect.size() - 1;
	//lines
	for(long ti = 0; ti < hexaLine.size(); ti ++){
		std::vector<long> line = {tempElem.cornNode[hexaLine[ti][0]], 
			tempElem.cornNode[hexaLine[ti][1]]};
		sort(line.begin(), line.end());
		auto iterLius = lineUsed.find(line);
		if(iterLius == lineUsed.end()){
			std::set<long> tempList;
			tempList.insert(elemNumb);
			lineUsed.emplace(line, tempList);
		}
		else{
			(iterLius->second).insert(elemNumb);
		}
	}
	//faces
	for(long ti = 0; ti < hexaFace.size(); ti ++){
		std::vector<long> face = {tempElem.cornNode[hexaFace[ti][0]], 
			tempElem.cornNode[hexaFace[ti][1]], 
			tempElem.cornNode[hexaFace[ti][2]], 
			tempElem.cornNode[hexaFace[ti][3]]};
		sort(face.begin(), face.end());
		auto iterFaus = faceUsed.find(face);
		if(iterFaus == faceUsed.end()){
			std::set<long> tempList;
			tempList.insert(elemNumb);
			faceUsed.emplace(face, tempList);
		}
		else{
			(iterFaus->second).insert(elemNumb);
		}
	}
	//
	mgpi.maxiLeve = std::max(mgpi.maxiLeve, tempElem.level);
	return elemNumb;
}

long MULTIGRID::REFINE(std::set<long> &spliElem, 
	const std::map<long,std::set<long>> &spliFlag, 
	const std::map<std::vector<long>,COOR> &planSurf){
	//
	GRLE_CHECK(spliElem);
	//
	//refinement template
	VECTOR3L refiTemp_1 = {{
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5}, 
		{2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, 
		{0, 3, 7, 4}, {1, 2, 6, 5}, {0, 4, 5, 1}, {3, 7, 6, 2}, {0, 1, 2, 3}, {4, 5, 6, 7}, 
		{0, 1, 2, 3, 4, 5, 6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, 
		{4, 5}, {5, 6}, {6, 7}, {7, 4}, 
		{0, 1, 2, 3}, {4, 5, 6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 3}, {3, 7}, {7, 4}, {4, 0}, 
		{1, 2}, {2, 6}, {6, 5}, {5, 1}, 
		{0, 3, 7, 4}, {1, 2, 6, 5}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 4}, {4, 5}, {5, 1}, {1, 0}, 
		{3, 7}, {7, 6}, {6, 2}, {2, 3}, 
		{0, 4, 5, 1}, {3, 7, 6, 2}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 1}, {2, 3}, 
		{4, 5}, {6, 7}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 3}, {1, 2}, 
		{4, 7}, {5, 6}
	}, {
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, 
		{0, 4}, {1, 5}, 
		{3, 7}, {2, 6}
	}};
	VECTOR3L refiTemp_2 = {{
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, 
		{0, 0, 1}, {2, 0, 1}, {2, 2, 1}, {0, 2, 1}, {1, 0, 2}, {2, 1, 2}, 
		{1, 2, 2}, {0, 1, 2}, {0, 1, 1}, {2, 1, 1}, {1, 0, 1}, {1, 2, 1}, 
		{1, 1, 0}, {1, 1, 2}, {1, 1, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0}, 
		{1, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 1, 2}, {1, 1, 0}, {1, 1, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 1, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 1}, 
		{2, 1, 0}, {2, 2, 1}, {2, 1, 2}, {2, 0, 1}, {0, 1, 1}, {2, 1, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 0, 1}, {1, 0, 2}, {2, 0, 1}, {1, 0, 0}, 
		{0, 2, 1}, {1, 2, 2}, {2, 2, 1}, {1, 2, 0}, {1, 0, 1}, {1, 2, 1}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {1, 0, 0}, {1, 2, 0}, {1, 0, 2}, {1, 2, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 1, 0}, {2, 1, 0}, {0, 1, 2}, {2, 1, 2}
	}, {
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0}, {0, 0, 2}, {2, 0, 2}, 
		{2, 2, 2}, {0, 2, 2}, {0, 0, 1}, {2, 0, 1}, {0, 2, 1}, {2, 2, 1}
	}};
	//refined element 8 sub-elements, 8 nodes, 3 coordinates
	VECTOR4L refiElem_1 = {{
		{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}}, 
		{{1,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,0,1}, {2,0,1}, {2,1,1}, {1,1,1}}, 
		{{0,1,0}, {1,1,0}, {1,2,0}, {0,2,0}, {0,1,1}, {1,1,1}, {1,2,1}, {0,2,1}}, 
		{{1,1,0}, {2,1,0}, {2,2,0}, {1,2,0}, {1,1,1}, {2,1,1}, {2,2,1}, {1,2,1}}, 
		{{0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {0,0,2}, {1,0,2}, {1,1,2}, {0,1,2}}, 
		{{1,0,1}, {2,0,1}, {2,1,1}, {1,1,1}, {1,0,2}, {2,0,2}, {2,1,2}, {1,1,2}}, 
		{{0,1,1}, {1,1,1}, {1,2,1}, {0,2,1}, {0,1,2}, {1,1,2}, {1,2,2}, {0,2,2}}, 
		{{1,1,1}, {2,1,1}, {2,2,1}, {1,2,1}, {1,1,2}, {2,1,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,2}, {1,0,2}, {1,1,2}, {0,1,2}}, 
		{{1,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,0,2}, {2,0,2}, {2,1,2}, {1,1,2}}, 
		{{0,1,0}, {1,1,0}, {1,2,0}, {0,2,0}, {0,1,2}, {1,1,2}, {1,2,2}, {0,2,2}}, 
		{{1,1,0}, {2,1,0}, {2,2,0}, {1,2,0}, {1,1,2}, {2,1,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,1,0}, {0,1,0}, {0,0,1}, {2,0,1}, {2,1,1}, {0,1,1}}, 
		{{0,1,0}, {2,1,0}, {2,2,0}, {0,2,0}, {0,1,1}, {2,1,1}, {2,2,1}, {0,2,1}}, 
		{{0,0,1}, {2,0,1}, {2,1,1}, {0,1,1}, {0,0,2}, {2,0,2}, {2,1,2}, {0,1,2}}, 
		{{0,1,1}, {2,1,1}, {2,2,1}, {0,2,1}, {0,1,2}, {2,1,2}, {2,2,2}, {0,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,2,0}, {0,2,0}, {0,0,1}, {1,0,1}, {1,2,1}, {0,2,1}}, 
		{{0,0,1}, {1,0,1}, {1,2,1}, {0,2,1}, {0,0,2}, {1,0,2}, {1,2,2}, {0,2,2}}, 
		{{1,0,0}, {2,0,0}, {2,2,0}, {1,2,0}, {1,0,1}, {2,0,1}, {2,2,1}, {1,2,1}}, 
		{{1,0,1}, {2,0,1}, {2,2,1}, {1,2,1}, {1,0,2}, {2,0,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {1,0,0}, {1,2,0}, {0,2,0}, {0,0,2}, {1,0,2}, {1,2,2}, {0,2,2}}, 
		{{1,0,0}, {2,0,0}, {2,2,0}, {1,2,0}, {1,0,2}, {2,0,2}, {2,2,2}, {1,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,1,0}, {0,1,0}, {0,0,2}, {2,0,2}, {2,1,2}, {0,1,2}}, 
		{{0,1,0}, {2,1,0}, {2,2,0}, {0,2,0}, {0,1,2}, {2,1,2}, {2,2,2}, {0,2,2}}
	}, {
		{{0,0,0}, {2,0,0}, {2,2,0}, {0,2,0}, {0,0,1}, {2,0,1}, {2,2,1}, {0,2,1}}, 
		{{0,0,1}, {2,0,1}, {2,2,1}, {0,2,1}, {0,0,2}, {2,0,2}, {2,2,2}, {0,2,2}}
	}};
	std::set<long> spliElem_0;
	for(const auto& iterSpel : spliElem){
		//
		std::vector<long> befoCorn = elemVect[iterSpel].cornNode;
		long s = elemVect[iterSpel].refiPatt;
		//
		VECTOR3L AfteSpli(3, VECTOR2L(3, std::vector<long>(3)));
		for(long tk = 0; tk < refiTemp_1[s].size(); tk ++){
			//
			if(refiTemp_1[s][tk].size() == 1){
				AfteSpli[refiTemp_2[s][tk][0]][refiTemp_2[s][tk][1]][refiTemp_2[s][tk][2]] 
					= befoCorn[refiTemp_1[s][tk][0]];
				continue;
			}
			//
			std::vector<long> origNode(refiTemp_1[s][tk].size());
			for(long ti = 0; ti < refiTemp_1[s][tk].size(); ti ++){
				origNode[ti] = befoCorn[refiTemp_1[s][tk][ti]];
			}
			sort(origNode.begin(), origNode.end());
			auto iterPlsu = planSurf.find(origNode);
			if(iterPlsu != planSurf.end()){
				AfteSpli[refiTemp_2[s][tk][0]][refiTemp_2[s][tk][1]][refiTemp_2[s][tk][2]] 
					= TRY_ADD_NODE(iterPlsu->second);
			}
			else{
				COOR tempCoor(0.0, 0.0, 0.0);
				for(long ti = 0; ti < origNode.size(); ti ++){
					tempCoor = tempCoor + nodeCoor[origNode[ti]];
				}
				tempCoor = tempCoor / origNode.size();
				AfteSpli[refiTemp_2[s][tk][0]][refiTemp_2[s][tk][1]][refiTemp_2[s][tk][2]] 
					= TRY_ADD_NODE(tempCoor);
			}
		}
		//
		elemVect[iterSpel].children.resize(8);
		for(long ti = 0; ti < refiElem_1[s].size(); ti ++){
			TREE_ELEM tempElem;
			tempElem.parent = iterSpel;
			tempElem.cornNode = {
			AfteSpli[refiElem_1[s][ti][0][0]][refiElem_1[s][ti][0][1]][refiElem_1[s][ti][0][2]], 
			AfteSpli[refiElem_1[s][ti][1][0]][refiElem_1[s][ti][1][1]][refiElem_1[s][ti][1][2]], 
			AfteSpli[refiElem_1[s][ti][2][0]][refiElem_1[s][ti][2][1]][refiElem_1[s][ti][2][2]], 
			AfteSpli[refiElem_1[s][ti][3][0]][refiElem_1[s][ti][3][1]][refiElem_1[s][ti][3][2]], 
			AfteSpli[refiElem_1[s][ti][4][0]][refiElem_1[s][ti][4][1]][refiElem_1[s][ti][4][2]], 
			AfteSpli[refiElem_1[s][ti][5][0]][refiElem_1[s][ti][5][1]][refiElem_1[s][ti][5][2]], 
			AfteSpli[refiElem_1[s][ti][6][0]][refiElem_1[s][ti][6][1]][refiElem_1[s][ti][6][2]], 
			AfteSpli[refiElem_1[s][ti][7][0]][refiElem_1[s][ti][7][1]][refiElem_1[s][ti][7][2]]
			};
			tempElem.level = elemVect[iterSpel].level + 1;
			tempElem.refiPatt = 7;
			tempElem.children.resize(0);
			elemVect[iterSpel].children[ti] = ADD_ELEMENT(tempElem);
		}
		//
		auto iterSpfl = spliFlag.find(iterSpel);
		if(iterSpfl == spliFlag.end()){
			continue;
		}
		for(const auto& iterCoel : (iterSpfl->second)){
			spliElem_0.insert(elemVect[iterSpel].children[iterCoel]);
		}
	}
	spliElem = spliElem_0;
	return 1;
}

long MULTIGRID::GRLE_CHECK(std::set<long> &spliElem){
	//when refine a line L of spliElem: spliElem's parent does not have this line L;
	//denote spliElem's parent as P, L is a sub segment of one line in the 12 lines of P;
	//P's neighbor should be refined in the same way of as P;
	VECTOR4L pareLine = {{
		{{0, 1}, {0, 3}, {0, 4}}, {{1, 0}, {1, 2}, {1, 5}}, 
		{{3, 0}, {3, 2}, {3, 7}}, {{2, 1}, {2, 3}, {2, 6}}, 
		{{4, 0}, {4, 5}, {4, 7}}, {{5, 1}, {5, 4}, {5, 6}}, 
		{{7, 3}, {7, 6}, {7, 4}}, {{6, 2}, {6, 5}, {6, 7}}
	}, {
		{{0, 1}, {0, 3}, {4, 5}, {4, 7}}, 
		{{1, 0}, {1, 2}, {5, 4}, {5, 6}}, 
		{{3, 0}, {3, 2}, {7, 4}, {7, 6}}, 
		{{2, 1}, {2, 3}, {6, 5}, {6, 7}}
	}, {
		{{0, 3}, {0, 4}, {1, 2}, {1, 5}}, 
		{{3, 0}, {3, 7}, {2, 1}, {2, 6}}, 
		{{4, 0}, {4, 7}, {5, 1}, {5, 6}}, 
		{{7, 3}, {7, 4}, {6, 2}, {6, 5}}
	}, {
		{{0, 1}, {0, 4}, {3, 2}, {3, 7}}, 
		{{4, 0}, {4, 5}, {7, 3}, {7, 6}}, 
		{{1, 0}, {1, 5}, {2, 3}, {2, 6}}, 
		{{5, 1}, {5, 4}, {6, 2}, {6, 7}}
	}, {
		{{0, 1}, {4, 5}, {2, 3}, {7, 6}}, 
		{{0, 1}, {4, 5}, {2, 3}, {7, 6}}
	}, {
		{{0, 3}, {1, 2}, {4, 7}, {5, 6}}, 
		{{0, 3}, {1, 2}, {4, 7}, {5, 6}}
	}, {
		{{0, 4}, {1, 5}, {3, 7}, {2, 6}}, 
		{{0, 4}, {1, 5}, {3, 7}, {2, 6}}
	}};
	VECTOR4L pareFace = {{
		{{0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1}}, 
		{{1, 2, 3, 0}, {1, 2, 6, 5}, {1, 0, 4, 5}}, 
		{{3, 0, 1, 2}, {3, 7, 4, 0}, {3, 7, 6, 2}}, 
		{{2, 3, 0, 1}, {2, 6, 5, 1}, {2, 3, 7, 6}}, 
		{{4, 5, 6, 7}, {4, 0, 3, 7}, {4, 5, 1, 0}}, 
		{{5, 6, 7, 4}, {5, 1, 2, 6}, {5, 1, 0, 4}}, 
		{{7, 4, 5, 6}, {7, 4, 0, 3}, {7, 6, 2, 3}}, 
		{{6, 7, 4, 5}, {6, 5, 1, 2}, {6, 2, 3, 7}}
	}, {
		{{0, 1, 2, 3}, {0, 3, 7, 4}, {0, 4, 5, 1}, {4, 5, 6, 7}}, 
		{{1, 2, 3, 0}, {1, 2, 6, 5}, {1, 0, 4, 5}, {5, 6, 7, 4}}, 
		{{3, 0, 1, 2}, {3, 7, 4, 0}, {3, 7, 6, 2}, {7, 4, 5, 6}}, 
		{{2, 3, 0, 1}, {2, 6, 5, 1}, {2, 3, 7, 6}, {6, 7, 4, 5}}
	}, {
		{{0, 3, 7, 4}, {1, 2, 6, 5}, {0, 4, 5, 1}, {0, 1, 2, 3}}, 
		{{3, 7, 4, 0}, {2, 6, 5, 1}, {3, 7, 6, 2}, {3, 0, 1, 2}}, 
		{{4, 0, 3, 7}, {5, 1, 2, 6}, {4, 5, 6, 7}, {4, 5, 1, 0}}, 
		{{7, 4, 0, 3}, {6, 5, 1, 2}, {7, 6, 2, 3}, {7, 4, 5, 6}}
	}, {
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {0, 1, 2, 3}}, 
		{{4, 5, 1, 0}, {7, 6, 2, 3}, {4, 5, 6, 7}, {4, 0, 3, 7}}, 
		{{1, 0, 4, 5}, {2, 3, 7, 6}, {1, 2, 3, 0}, {1, 2, 6, 5}}, 
		{{5, 1, 0, 4}, {6, 2, 3, 7}, {5, 6, 7, 4}, {5, 1, 2, 6}}
	}, {
		{{0, 1, 2, 3}, {0, 4, 5, 1}, {4, 5, 6, 7}, {3, 7, 6, 2}}, 
		{{0, 1, 2, 3}, {0, 4, 5, 1}, {4, 5, 6, 7}, {3, 7, 6, 2}}
	}, {
		{{0, 3, 7, 4}, {4, 5, 6, 7}, {1, 2, 6, 5}, {0, 1, 2, 3}}, 
		{{0, 3, 7, 4}, {4, 5, 6, 7}, {1, 2, 6, 5}, {0, 1, 2, 3}}
	}, {
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {1, 2, 6, 5}}, 
		{{0, 4, 5, 1}, {3, 7, 6, 2}, {0, 3, 7, 4}, {1, 2, 6, 5}}
	}};
	std::set<long> spliElem_0 = spliElem;
	while(true){
		std::set<long> spliElem_A;
		for(const auto& iterSpvo : spliElem_0){
			TREE_ELEM tempElem = elemVect[iterSpvo];
			if(tempElem.parent != -1){
				TREE_ELEM tempPare = elemVect[tempElem.parent];
				long p = tempPare.refiPatt;
				long isub;
				for(isub = 0; isub < tempPare.children.size(); isub ++){
					if(tempPare.children[isub] == iterSpvo){
						break;
					}
				}
				//
				for(long tj = 0; tj < pareLine[p][isub].size(); tj ++){
					std::vector<long> tempLine = 
						{tempPare.cornNode[pareLine[p][isub][tj][0]], 
						 tempPare.cornNode[pareLine[p][isub][tj][1]]};
					sort(tempLine.begin(), tempLine.end());
					auto iterLius = lineUsed.find(tempLine);
					//iterLius != lineUsed.end()
					for(const auto& iterAdel : (iterLius->second)){
						if(elemVect[iterAdel].children.size() == 0){
							auto iterSpel = spliElem.find(iterAdel);
							if(iterSpel == spliElem.end()){
								spliElem_A.insert(iterAdel);
								elemVect[iterAdel].refiPatt = 0;
							}
						}
					}
				}
				//
				for(long tj = 0; tj < pareFace[p][isub].size(); tj ++){
					std::vector<long> tempFace = 
						{tempPare.cornNode[pareFace[p][isub][tj][0]], 
						 tempPare.cornNode[pareFace[p][isub][tj][1]], 
						 tempPare.cornNode[pareFace[p][isub][tj][2]], 
						 tempPare.cornNode[pareFace[p][isub][tj][3]]};
					sort(tempFace.begin(), tempFace.end());
					auto iterFaus = faceUsed.find(tempFace);
					//iterFaus != faceUsed.end()
					for(const auto& iterAdel : (iterFaus->second)){
						if(elemVect[iterAdel].children.size() == 0){
							auto iterSpel = spliElem.find(iterAdel);
							if(iterSpel == spliElem.end()){
								spliElem_A.insert(iterAdel);
								elemVect[iterAdel].refiPatt = 0;
							}
						}
					}
				}
			}
		}
		if(spliElem_A.size() == 0){
			break;
		}
		else{
			spliElem_0 = spliElem_A;
			spliElem.insert(spliElem_A.begin(), spliElem_A.end());
		}
	}
	return 1;
}

long MULTIGRID::OUTPUT_ELEMENT(long fileIden){
	OUTPUT_TIME("MULTIGRID::OUTPUT_ELEMENT");
	std::stringstream tempStre;
	tempStre << "resuNode_" << fileIden << ".txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(const auto& iterNoco : nodeCoor){
		tempOfst << std::setw(30) << (iterNoco.second)[0]
			<< std::setw(30) << (iterNoco.second)[1]
			<< std::setw(30) << (iterNoco.second)[2] << std::endl;
	}
	tempOfst.close();
	tempStre << "resuElem_" << fileIden << ".txt";
	tempOfst.open(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	for(long ti = 0; ti < elemVect.size(); ti ++){
		if(elemVect[ti].children.size() == 0){
			for(long tj = 0; tj < elemVect[ti].cornNode.size(); tj ++){
				tempOfst << std::setw(10) << elemVect[ti].cornNode[tj];
			}
			tempOfst << std::endl;
		}
	}
	tempOfst.close();
	return 1;
}

long MULTIGRID::RIGI_ROTR(Eigen::Matrix3d rotaMatr, Eigen::Vector3d tranVect){
	coorNode.clear();
	for(auto& iterNoco : nodeCoor){
		Eigen::Vector3d tempCoor;
		tempCoor << (iterNoco.second)[0], (iterNoco.second)[1], (iterNoco.second)[2];
		tempCoor = rotaMatr * tempCoor + tranVect;
		iterNoco.second = {tempCoor(0), tempCoor(1), tempCoor(2)};
		coorNode.emplace(iterNoco.second, iterNoco.first);
	}
	return 1;
}

long MULTIGRID::PATCH(const MOVL &tempMovl){
	//the level difference of two adjcent elements <= 1
	//no need to successively from level 0 to level mgpi.maxiLeve - 1?????
	for(const auto& iterIntr : tempMovl){
		long outpNode = iterIntr.second;
		COOR outpCoor = {0.0, 0.0, 0.0};
		for(long ti = 0; ti < (iterIntr.first).size(); ti ++){
			long inpuNode = (iterIntr.first)[ti];
			COOR inpuCoor = nodeCoor[inpuNode];
			outpCoor = outpCoor + inpuCoor;
		}
		outpCoor = outpCoor / (iterIntr.first).size();
		//
		auto iterNoco = nodeCoor.find(outpNode);
		if(iterNoco == nodeCoor.end()){
			std::cout << "Error 1 in MULTIGRID::PATCH" << std::endl;
		}
		COOR origCoor = iterNoco->second;
		nodeCoor.erase(iterNoco);
		nodeCoor.emplace(outpNode, outpCoor);
		//
		auto iterCono = coorNode.find(origCoor);
		if(iterCono == coorNode.end()){
			std::cout << "Error 2 in MULTIGRID::PATCH" << std::endl;
		}
		coorNode.erase(iterCono);
		coorNode.emplace(outpCoor, outpNode);
		//
		finoCono.emplace(iterIntr.second, iterIntr.first);
	}
	conoFino = tempMovl;
	return 1;
}

long MULTIGRID::TRANSFER(){
	OUTPUT_TIME("MULTIGRID::TRANSFER");
	//
	VECTOR3L elemLine = {{
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
	VECTOR3L elemFace = {{
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
	//inner transfer in the same level
	//inter-transfer in different levels
	std::vector<MOVL> ininTran;
	ininTran.resize(mgpi.maxiLeve + 1);
	//avoid duplicate nodes
	std::vector<std::set<long>> leveNode_s;
	leveNode_s.resize(mgpi.maxiLeve + 2);
	for(long ti = 0; ti < elemVect.size(); ti ++){
		if(elemVect[ti].level == 0){
			for(long tj = 0; tj < 8; tj ++){
				leveNode_s[0].insert(elemVect[ti].cornNode[tj]);
			}
		}
		if(elemVect[ti].children.size() != 0){
			long s = elemVect[ti].refiPatt;
			if(s == 0){
				//
				std::vector<long> tempCorn = {
					elemVect[ti].cornNode[0], elemVect[ti].cornNode[1], 
					elemVect[ti].cornNode[2], elemVect[ti].cornNode[3], 
					elemVect[ti].cornNode[4], elemVect[ti].cornNode[5], 
					elemVect[ti].cornNode[6], elemVect[ti].cornNode[7]
				};
				sort(tempCorn.begin(), tempCorn.end());
				ininTran[elemVect[ti].level].insert(MOVL::value_type(
					tempCorn, 
					elemVect[elemVect[ti].children[0]].cornNode[6])
				);
				leveNode_s[elemVect[ti].level + 1].insert(
					elemVect[elemVect[ti].children[0]].cornNode[6]
				);
			}
			//
			for(long tj = 0; tj < elemLine[s].size(); tj ++){
				std::vector<long> tempLine = {elemVect[ti].cornNode[elemLine[s][tj][0]], 
					elemVect[ti].cornNode[elemLine[s][tj][1]]};
				sort(tempLine.begin(), tempLine.end());
				auto iterLius = lineUsed.find(tempLine);
				//iterLius != lineUsed.end();
				bool tempFlag = false;
				for(const auto& iterCoel : (iterLius->second)){
					if(elemVect[iterCoel].children.size() == 0){
						tempFlag = true;
						break;
					}
				}
				long tempNode = elemVect[
						elemVect[ti].children[elemLine[s][tj][2]]
					].cornNode[elemLine[s][tj][3]];
				if(tempFlag == true){
					ininTran[mgpi.maxiLeve].insert(MOVL::value_type(tempLine, tempNode));
					leveNode_s[mgpi.maxiLeve + 1].insert(tempNode);
				}
				else{
					ininTran[elemVect[ti].level].insert(MOVL::value_type(tempLine, tempNode));
					leveNode_s[elemVect[ti].level + 1].insert(tempNode);
				}
			}
			//
			for(long tj = 0; tj < elemFace[s].size(); tj ++){
				std::vector<long> tempFace = {elemVect[ti].cornNode[elemFace[s][tj][0]], 
					elemVect[ti].cornNode[elemFace[s][tj][1]], 
					elemVect[ti].cornNode[elemFace[s][tj][2]], 
					elemVect[ti].cornNode[elemFace[s][tj][3]]};
				sort(tempFace.begin(), tempFace.end());
				auto iterFaus = faceUsed.find(tempFace);
				//iterFaus != faceUsed.end();
				bool tempFlag = false;
				for(const auto& iterCoel : (iterFaus->second)){
					if(elemVect[iterCoel].children.size() == 0){
						tempFlag = true;
						break;
					}
				}
				long tempNode = elemVect[
						elemVect[ti].children[elemFace[s][tj][4]]
					].cornNode[elemFace[s][tj][5]];
				if(tempFlag == true){
					ininTran[mgpi.maxiLeve].insert(MOVL::value_type(tempFace, tempNode));
					leveNode_s[mgpi.maxiLeve + 1].insert(tempNode);
				}
				else{
					ininTran[elemVect[ti].level].insert(MOVL::value_type(tempFace, tempNode));
					leveNode_s[elemVect[ti].level + 1].insert(tempNode);
				}
			}
		}
	}
	//
	PATCH(ininTran[mgpi.maxiLeve]);
	//
	leveNode.resize(mgpi.maxiLeve + 2);
	VECT_RESI(nodeLepo, nodeCoor.size(), 2);
	posiNode.resize(nodeCoor.size());
	scalProl.resize(mgpi.maxiLeve + 1);
	for(long ti = 0; ti <= mgpi.maxiLeve + 1; ti ++){
		for(const auto& iterLeno : leveNode_s[ti]){
			auto iterCono = coupNode.find(iterLeno);
			if(iterLeno == coupReps){
				leveNode[0].push_back(iterLeno);
			}
			else if(iterCono != coupNode.end()){
				leveNode[mgpi.maxiLeve + 1].push_back(iterLeno);
			}
			else{
				leveNode[ti].push_back(iterLeno);
			}
		}
	}
	long accuNumb = 0;
	for(long ti = 0; ti <= mgpi.maxiLeve + 1; ti ++){
		for(long tj = 0; tj < leveNode[ti].size(); tj ++){
			nodeLepo[leveNode[ti][tj]][0] = ti;
			nodeLepo[leveNode[ti][tj]][1] = accuNumb + tj;
			posiNode[accuNumb + tj] = leveNode[ti][tj];
		}
		accuNumb += leveNode[ti].size();
	}
	accuNumb = 0;
	for(long ti = 0; ti <= mgpi.maxiLeve; ti ++){
		std::vector<Eigen::Triplet<double>> tempList;
		for(long tj = 0; tj < accuNumb + leveNode[ti].size(); tj ++){
			tempList.emplace_back(tj, tj, 1.0);
		}
		for(const auto& iterIntr : ininTran[ti]){
			//
			long resuLepo = nodeLepo[iterIntr.second][1];
			auto iterCono = coupNode.find(iterIntr.second);
			if(iterCono != coupNode.end()){
				continue;
			}
			for(long tj = 0; tj < (iterIntr.first).size(); tj ++){
				long prolLepo;
				auto iterCono = coupNode.find((iterIntr.first)[tj]);
				if(iterCono != coupNode.end()){
					prolLepo = nodeLepo[coupReps][1];
				}
				else{
					prolLepo = nodeLepo[(iterIntr.first)[tj]][1];
				}
				tempList.emplace_back(resuLepo, prolLepo, 1.0 / (iterIntr.first).size());
			}
		}
		if(ti == mgpi.maxiLeve){
			for(const auto &iterCono : coupNode){
				tempList.emplace_back(nodeLepo[iterCono][1], nodeLepo[coupReps][1], 1.0);
			}
		}
		scalProl[ti].resize(accuNumb + leveNode[ti].size() + leveNode[ti + 1].size(), 
			accuNumb + leveNode[ti].size());
		scalProl[ti].setFromTriplets(tempList.begin(), tempList.end());
		tempList.clear();
		accuNumb += leveNode[ti].size();
	}
	return 1;
}

long MULTIGRID::STIF_MATR(){
	OUTPUT_TIME("MULTIGRID::STIF_MATR");
	//
	double mateLame[2];
	mateLame[0] = mateElas * matePois / (1.0 + matePois) / (1.0 - 2.0 * matePois);
	mateLame[1] = mateElas / 2.0 / (1.0 + matePois);
	elatMatr << 2.0 * mateLame[1] + mateLame[0], mateLame[0], mateLame[0], 0.0, 0.0, 0.0,
		mateLame[0], 2.0 * mateLame[1] + mateLame[0], mateLame[0], 0.0, 0.0, 0.0,
		mateLame[0], mateLame[0], 2.0 * mateLame[1] + mateLame[0], 0.0, 0.0, 0.0,
    	0.0, 0.0, 0.0, mateLame[1], 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, mateLame[1], 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, mateLame[1];
	//
	long partSize = 50000;
	long partNumb = elemVect.size() / partSize;
	if(partNumb * partSize < elemVect.size()){
		partNumb ++;
	}
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> partStif(partNumb);
	#pragma omp parallel for
	for(long tp = 0; tp < partNumb; tp ++){
		long star_tp = tp * partSize;
		long endi_tp = star_tp + partSize;
		if(endi_tp > elemVect.size()){
			endi_tp = elemVect.size();
		}
		if(star_tp >= endi_tp){
			continue;
		}
		std::vector<Eigen::Triplet<double>> stifList;
		stifList.reserve((endi_tp - star_tp) * 24 * 24);
		for(long ti = star_tp; ti < endi_tp; ti ++){
			if(ti % 10000 == 0){
				std::cout << "Element stiffness #" << ti << std::endl;
			}
			if(elemVect[ti].children.size() > 0){
				continue;
			}
			Eigen::MatrixXd exyz = Eigen::MatrixXd::Zero(8, 3);
			for(long tj = 0; tj < 8; tj ++){
				auto iterNoco = nodeCoor.find(elemVect[ti].cornNode[tj]);
				for(long tk = 0; tk < 3; tk ++){
					exyz(tj,tk) = (iterNoco->second)[tk];
				}
			}
			//
			Eigen::MatrixXd K_E = Eigen::MatrixXd::Zero(24, 24);
			for(long tj = 0; tj < trilQuad.numbNgip; tj ++){
				Eigen::MatrixXd ngipJact = trilQuad.pnpxNgip[tj] * exyz;
				Eigen::MatrixXd B_L = Eigen::MatrixXd::Zero(6, 24);
				for(long tk = 0; tk < 8; tk ++){
					Eigen::Vector3d pnpX = ngipJact.fullPivLu().solve(
						trilQuad.pnpxNgip[tj].block(0,tk,3,1)
					);
					B_L.block(0, 3 * tk, 6, 3) << 
						pnpX(0), 0.0, 0.0, 
						0.0, pnpX(1), 0.0, 
						0.0, 0.0, pnpX(2), 
						pnpX(1), pnpX(0), 0.0, 
						0.0, pnpX(2), pnpX(1), 
						pnpX(2), 0.0, pnpX(0);
				}
				K_E += trilQuad.niwfNgip[tj] * ngipJact.determinant() 
					* (B_L.transpose() * elatMatr * B_L);
			}
			for(long tj = 0; tj < 8; tj ++){
				for(long tk = 0; tk  < 3; tk ++){
					long row_jk = 3 * elemVect[ti].cornNode[tj] + tk;
					for(long tm = 0; tm < 8; tm ++){
						for(long tn = 0; tn  < 3; tn ++){
							long col_mn = 3 * elemVect[ti].cornNode[tm] + tn;
							stifList.emplace_back(
								row_jk, col_mn, K_E(3 * tj + tk , 3 * tm + tn)
							);
						}
					}
				}
			}
		}
		partStif[tp].resize(3 * nodeCoor.size(), 3 * nodeCoor.size());
		partStif[tp].setFromTriplets(stifList.begin(), stifList.end());
		stifList.clear();
	}
	origStif.resize(mgpi.maxiLeve + 2);
	origStif[mgpi.maxiLeve + 1] = partStif[0];
	for(long tp = 1; tp < partNumb; tp ++){
		origStif[mgpi.maxiLeve + 1] += partStif[tp];
	}
	return 1;
}

long MULTIGRID::LOAD_ACCU(long free_ti, double inpuLoad){
	auto iterCodo = consDofv.find(free_ti);
	if(iterCodo == consDofv.end()){
		auto iterExfo = exteForc.find(free_ti);
		if(iterExfo == exteForc.end()){
			exteForc.emplace(free_ti, inpuLoad);
		}
		else{
			iterExfo->second = iterExfo->second + inpuLoad;
		}
	}
	else if(debuMode >= 1){
		std::cout << "Error 1 in MULTIGRID::LOAD_ACCU!" << std::endl;
		return 0;
	}
	return 1;
}

long MULTIGRID::CONSTRAINT(long geomMult){
	OUTPUT_TIME("MULTIGRID::CONSTRAINT");
	//
	std::vector<Eigen::Triplet<double>> rotaList;
	for(long ti = 0; ti < nodeCoor.size(); ti ++){
		auto iterNoro = nodeRota.find(ti);
		if(iterNoro == nodeRota.end()){
			for(long tk = 0; tk < 3; tk ++){
				rotaList.emplace_back(3 * ti + tk, 3 * ti + tk, 1.0);
			}
		}
		else{
			for(long tj = 0; tj < 3; tj ++){
				for(long tk = 0; tk < 3; tk ++){
					rotaList.emplace_back(3 * ti + tj, 3 * ti + tk, (iterNoro->second)(tj,tk));
				}
			}
		}
	}
	Eigen::SparseMatrix<double,Eigen::RowMajor> rotaMatr;
	rotaMatr.resize(3 * nodeCoor.size(), 3 * nodeCoor.size());
	rotaMatr.setFromTriplets(rotaList.begin(), rotaList.end());
	origStif[mgpi.maxiLeve + 1] = rotaMatr.transpose() * origStif[mgpi.maxiLeve + 1] * rotaMatr;
	//
	std::vector<Eigen::Triplet<double>> tranList, tranList_s;
	for(long ti = 0; ti < nodeCoor.size(); ti ++){
		for(long tj = 0; tj < 3; tj ++){
			tranList.emplace_back(3 * ti + tj, 3 * nodeLepo[ti][1] + tj, 1.0);
		}
		tranList_s.emplace_back(ti, nodeLepo[ti][1], 1.0);
	}
	earlTran.resize(3 * nodeCoor.size(), 3 * nodeCoor.size());
	earlTran.setFromTriplets(tranList.begin(), tranList.end());
	scalEarl.resize(nodeCoor.size(), nodeCoor.size());
	scalEarl.setFromTriplets(tranList_s.begin(), tranList_s.end());
	tranList.clear();
	tranList_s.clear();
	origStif[mgpi.maxiLeve + 1] = earlTran.transpose() * origStif[mgpi.maxiLeve + 1] * earlTran;
	//
	prolOper.resize(mgpi.maxiLeve + 1);
	for(long tl = ((geomMult == 1) ? 0 : mgpi.maxiLeve); tl <= mgpi.maxiLeve; tl ++){
		std::vector<Eigen::Triplet<double>> tempList;
		for(long ti = 0; ti < 3 * scalProl[tl].cols(); ti ++){
			tempList.emplace_back(ti, ti, 1.0);
		}
		for(long ti = scalProl[tl].cols(); ti < scalProl[tl].rows(); ti ++){
			for(RSPA_INNE iterScpr(scalProl[tl], ti); iterScpr; ++ iterScpr){
				long offsNode = posiNode[ti];
				long famoNode = posiNode[iterScpr.col()];
				auto iterOffs = nodeRota.find(offsNode);//offspring
				auto iterFamo = nodeRota.find(famoNode);//father, mother
				Eigen::Matrix3d tempRota = iterScpr.value() * Eigen::Matrix3d::Identity(3, 3);
				auto iterCono = coupNode.find(offsNode);
				if((iterOffs != nodeRota.end() || iterFamo != nodeRota.end()) 
					&& !(famoNode == coupReps && iterCono != coupNode.end())){
					if(iterOffs != nodeRota.end() && iterFamo != nodeRota.end()
						/*&& (iterOffs->second - iterFamo->second).norm() < 1.0E-10*/){
						//only one pattern of nodeRota
					}
					else{
						if(iterOffs != nodeRota.end()){
							tempRota *= (iterOffs->second).transpose();
						}
						if(iterFamo != nodeRota.end()){
							tempRota *= iterFamo->second;
						}
					}
				}
				for(long tj = 0; tj < 3; tj ++){
					for(long tk = 0; tk < 3; tk ++){
						tempList.emplace_back(
							3 * ti + tj, 3 * iterScpr.col() + tk, tempRota(tj,tk)
						);
					}
				}
			}
		}
		prolOper[tl].resize(3 * scalProl[tl].rows(), 3 * scalProl[tl].cols());
		prolOper[tl].setFromTriplets(tempList.begin(), tempList.end());
	}
	for(long ti = mgpi.maxiLeve; ti >= ((geomMult == 1) ? 0 : mgpi.maxiLeve); ti --){
		origStif[ti] = prolOper[ti].transpose() * origStif[ti + 1] * prolOper[ti];
	}
	//
	consFlag = Eigen::VectorXi::Ones(3 * nodeCoor.size());
	Eigen::VectorXd dispForc_0 = Eigen::VectorXd::Zero(3 * nodeCoor.size());
	for(const auto& iterCodo : consDofv){
		long tempNode = iterCodo.first / 3;
		long tempCdof = iterCodo.first % 3;
		//earlTran.transpose() * consFlag : double * long
		consFlag(3 * nodeLepo[tempNode][1] + tempCdof) = 0;
		dispForc_0(iterCodo.first) = iterCodo.second;
	}
	dispForc_0 = earlTran.transpose() * dispForc_0;//NO: nodeRota * 
	long tempNumb_0 = 0;
	dispForc = Eigen::VectorXd::Zero(consDofv.size());
	for(long ti = 0; ti < origStif[mgpi.maxiLeve].rows(); ti ++){
		if(consFlag(ti) == 0){
			dispForc(tempNumb_0) = dispForc_0(ti);
			tempNumb_0 ++;
		}
	}
	dispForc = dispForc.block(0,0,tempNumb_0,1);
	//
	consForc = Eigen::VectorXd::Zero(3 * nodeCoor.size());
	for(const auto& iterExfo : exteForc){
		consForc(iterExfo.first) += iterExfo.second;
	}
	//NO NEED: nodeRota * 
	consForc = prolOper[mgpi.maxiLeve].transpose() * earlTran.transpose() * consForc;
	//
	consOper.resize(mgpi.maxiLeve + 1);
	mgpi.consStif.resize(mgpi.maxiLeve + 1);
	for(long ti = mgpi.maxiLeve; ti >= ((geomMult == 1) ? 0 : mgpi.maxiLeve); ti --){
		std::vector<Eigen::Triplet<double>> consList;
		long tempNumb = 0;
		for(long tj = 0; tj < origStif[ti].rows(); tj ++){
			if(consFlag(tj) == 1){
				consList.emplace_back(tempNumb, tj, 1.0);
				tempNumb ++;
			}
		}
		consOper[ti].resize(tempNumb, origStif[ti].rows());
		consOper[ti].setFromTriplets(consList.begin(), consList.end());
		consList.clear();
		mgpi.consStif[ti] = consOper[ti] * origStif[ti] * consOper[ti].transpose();
		if(ti == mgpi.maxiLeve){
			long tempNumb = 0;
			for(long tj = 0; tj < origStif[ti].rows(); tj ++){
				if(consFlag(tj) == 0){
					consList.emplace_back(tj, tempNumb, 1.0);
					tempNumb ++;
				}
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> consMatr_1;
			consMatr_1.resize(origStif[ti].cols(), tempNumb);
			consMatr_1.setFromTriplets(consList.begin(), consList.end());
			consList.clear();
			consForc = consOper[ti] * consForc 
				- consOper[ti] * (origStif[ti] * (consMatr_1 * dispForc));
		}
	}
	//
	if(geomMult == 1){
		mgpi.realProl.resize(mgpi.maxiLeve);
		for(long ti = 0; ti < mgpi.maxiLeve; ti ++){
			mgpi.realProl[ti] = consOper[ti + 1] * prolOper[ti] * consOper[ti].transpose();
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		mgpi.ESTABLISH();
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
	return 1;
}

long MULTIGRID::ADDITIONAL_FORCE(Eigen::VectorXd &addiForc){
	addiForc = consOper[mgpi.maxiLeve] * prolOper[mgpi.maxiLeve].transpose() 
		* earlTran.transpose() * addiForc;
	return 1;
}

long MULTIGRID::OUTP_SUB1(const Eigen::VectorXd &resuDisp, Eigen::VectorXd &outpDisp){
	outpDisp = Eigen::VectorXd::Zero(origStif[mgpi.maxiLeve].rows());
	long tempNumb = 0;
	for(long ti = 0; ti < origStif[mgpi.maxiLeve].rows(); ti ++){
		if(consFlag(ti) == 1){
			outpDisp(ti) = resuDisp[tempNumb];
			tempNumb ++;
		}
	}
	for(const auto& iterCodo : consDofv){
		long tempNode = iterCodo.first / 3;
		long tempCdof = iterCodo.first % 3;
		if(nodeLepo[tempNode][0] <= mgpi.maxiLeve){
			outpDisp(3 * nodeLepo[tempNode][1] + tempCdof) = iterCodo.second;
		}
	}
	outpDisp = earlTran * prolOper[mgpi.maxiLeve] * outpDisp;
	return 1;
}

long MULTIGRID::OUTP_SUB2(const Eigen::VectorXd &outpDisp, long fileIden){
	OUTP_SUB2(outpDisp, std::to_string(fileIden));
	return 1;
}

long MULTIGRID::OUTP_SUB2(const Eigen::VectorXd &outpDisp, std::string fileIden){
	std::stringstream tempStre;
	tempStre << "resuDisp_" << fileIden << ".txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long ti = 0; ti < nodeCoor.size(); ti ++){
		Eigen::Vector3d disp_ti;
		disp_ti << outpDisp(3 * ti + 0), outpDisp(3 * ti + 1), outpDisp(3 * ti + 2);
		auto iterNoro = nodeRota.find(ti);
		if(iterNoro != nodeRota.end()){
			disp_ti = iterNoro->second * disp_ti;
		}
		tempOfst << std::setw(30) << disp_ti(0) << std::setw(30) << disp_ti(1) 
			<< std::setw(30) << disp_ti(2) << std::endl;
	}
	tempOfst.close();
	return 1;
}

Eigen::VectorXd MULTIGRID::OUTP_DISP(Eigen::VectorXd resuDisp, long fileIden){
	Eigen::VectorXd outpDisp;
	OUTP_SUB1(resuDisp, outpDisp);
	OUTP_SUB2(outpDisp, fileIden);
	return outpDisp;
}

long MULTIGRID::STRESS_RECOVERY(Eigen::VectorXd outpDisp, long fileIden){
	OUTPUT_TIME("MULTIGRID::STRESS_RECOVERY");
	for(long ti = 0; ti < nodeCoor.size(); ti ++){
		auto iterNoro = nodeRota.find(ti);
		if(iterNoro != nodeRota.end()){
			outpDisp.block(3 * ti + 0,0,3,1) = 
				iterNoro->second * outpDisp.block(3 * ti + 0,0,3,1).eval();
		}
	}
	//
	VECTOR3D nodeStre;
	VECT_RESI(nodeStre, nodeCoor.size(), 6);
	long tempCoun = 0;
	for(const auto &iterElem : elemVect){
		if(tempCoun % 10000 == 0){
			std::cout << "Element stress #" << tempCoun << std::endl;
		}
		tempCoun ++;
		if(iterElem.children.size() > 0){
			continue;
		}
		//
		Eigen::MatrixXd exyz = Eigen::MatrixXd::Zero(8, 3);
		Eigen::Vector<double,24> elemDisp;
		for(long tj = 0; tj < 8; tj ++){
			auto iterNoco = nodeCoor.find(iterElem.cornNode[tj]);
			for(long tk = 0; tk < 3; tk ++){
				elemDisp(tj * 3 + tk) = outpDisp(3 * iterElem.cornNode[tj] + tk);
				exyz(tj,tk) = (iterNoco->second)[tk];
			}
		}
		//
		Eigen::Matrix<double,8,6> inpoStre_0 = Eigen::MatrixXd::Zero(8,6);
		Eigen::Matrix<double,8,8> streMatr_0 = Eigen::MatrixXd::Zero(8,8);
		for(long tj = 0; tj < trilQuad.numbNgip; tj ++){
			Eigen::MatrixXd ngipJact = trilQuad.pnpxNgip[tj] * exyz;
			Eigen::MatrixXd B_L = Eigen::MatrixXd::Zero(6, 24);
			for(long tk = 0; tk < 8; tk ++){
				Eigen::Vector3d pnpX = ngipJact.fullPivLu().solve(
					trilQuad.pnpxNgip[tj].block(0,tk,3,1)
				);
				B_L.block(0, 3 * tk, 6, 3) << 
					pnpX(0), 0.0, 0.0, 
					0.0, pnpX(1), 0.0, 
					0.0, 0.0, pnpX(2), 
					pnpX(1), pnpX(0), 0.0, 
					0.0, pnpX(2), pnpX(1), 
					pnpX(2), 0.0, pnpX(0);
			}
			Eigen::MatrixXd stre_tj = (elatMatr * B_L * elemDisp).transpose();
			inpoStre_0 += trilQuad.niwfNgip[tj] * ngipJact.determinant() 
				* trilQuad.shfuNgip.block(tj,0,1,8).transpose() * stre_tj;
			streMatr_0 += trilQuad.niwfNgip[tj] * ngipJact.determinant() 
				* trilQuad.shfuNgip.block(tj,0,1,8).transpose() 
				* trilQuad.shfuNgip.block(tj,0,1,8);
		}
		Eigen::Matrix<double,8,6> elnoStre_0 = streMatr_0.fullPivLu().solve(inpoStre_0);
		for(long tj = 0; tj < 8; tj ++){
			for(long tk = 0; tk < 6; tk ++){
				nodeStre[iterElem.cornNode[tj]][tk].emplace_back(elnoStre_0(tj,tk));
			}
		}
		//
		for(long tj = 0; tj < 8; tj ++){
			auto iterFico = finoCono.find(iterElem.cornNode[tj]);
			if(iterFico != finoCono.end()){
				for(const auto &iterCono : iterFico->second){
					for(long tk = 0; tk < 6; tk ++){
						nodeStre[iterCono][tk].emplace_back(elnoStre_0(tj,tk));
					}
				}
			}
		}
		//
		for(long th = 0; th < 2; th ++){
			const VECTOR2L &hexaLifa = (th == 0) ? hexaLine : hexaFace;
			for(long tj = 0; tj < hexaLifa.size(); tj ++){
				std::vector<long> inpuNode(hexaLifa[tj].size());
				Eigen::Matrix<double,1,6> tempStre = Eigen::MatrixXd::Zero(1,6);
				for(long tk = 0; tk < hexaLifa[tj].size(); tk ++){
					inpuNode[tk] = iterElem.cornNode[hexaLifa[tj][tk]];
					tempStre += elnoStre_0.block(hexaLifa[tj][tk],0,1,6);
				}
				sort(inpuNode.begin(), inpuNode.end());
				auto iterCofi = conoFino.find(inpuNode);
				if(iterCofi != conoFino.end()){
					tempStre /= inpuNode.size();
					for(long tk = 0; tk < 6; tk ++){
						nodeStre[iterCofi->second][tk].emplace_back(tempStre(0,tk));
					}
				}
			}
		}
	}
	//
	std::stringstream tempStre;
	tempStre << "resuStre_" << fileIden << ".txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long ti = 0; ti < nodeCoor.size(); ti ++){
		std::vector<double> stre_ti(6, 0.0);
		for(long tj = 0; tj < 6; tj ++){
			stre_ti[tj] = 
				std::accumulate(nodeStre[ti][tj].begin(), nodeStre[ti][tj].end(), 0.0)
				/ nodeStre[ti][tj].size();
			tempOfst << std::setw(30) << stre_ti[tj];
		}
		double equiStre = sqrt((pow(stre_ti[0] - stre_ti[1], 2.0) 
			+ pow(stre_ti[1] - stre_ti[2], 2.0) + pow(stre_ti[0] - stre_ti[2], 2.0) 
			+ 6.0 * (pow(stre_ti[3], 2.0) + pow(stre_ti[4], 2.0) + pow(stre_ti[5], 2.0))) 
			/ 2.0);
		tempOfst << std::setw(30) << equiStre << std::endl;
	}
	tempOfst.close();
	return 1;
}

#endif
