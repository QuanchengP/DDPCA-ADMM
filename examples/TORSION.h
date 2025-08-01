
#include "../MCONTACT.h"

class TORSION: public MCONTACT{
public:
	//
	double axiaLeng;//axial length
	double inneRadi;//inner radius
	double outeRadi;//outer radius
	std::vector<long> diviNumb;//number of elements along different directions
	long globInho;//global inhomogeneous
	long globHomo;//global homogeneous
	std::vector<long> domaNumb;//number of decomposed domains
	//
	double torqLoad;//torque load
	std::vector<std::vector<CURVEDS>> blocSurf;
	long coloSett;//constraint load settings
	TORSION(long tempColo);
	//cylindrical coordinate average
	long COOR_AVER(const std::vector<COOR> &inpuCoor, COOR &outpCoor);
	long SUBR_PLSU(long tg, long ei, std::map<std::vector<long>,COOR> &planSurf);
	long SUBR_COLO(long tg);//displacement constraint and load
	long ESTA_SURF();//establish surface
	long MESH_DD();//domain decomposition
	long MESH_NODD();//no domain decomposition
	long SOLVE(long appsCont = 1, long noddPrec = 1);//solve for displacement
	long SOLVE_DD(long appsCont);
	long SOLVE_NODD(long noddPrec);
};

TORSION::TORSION(long tempColo){
	//
	mkdir("Torsion", 0777);
	outpDire = "Torsion/";
	//
	muscSett = (1 << 1);
	domaNumb = {1, 32, 8};//each diviNumb is divisible by each domaNumb
	doleMcsc.assign(domaNumb[0] * domaNumb[1] * domaNumb[2], 1);
	//
	axiaLeng = 0.1;
	inneRadi = 0.015;
	outeRadi = 0.025;//I_p = PI / 32.0 * (D^4 - d^4) = 0.005340707511103
	diviNumb = {1, 32, 8};//8, 32, 64
	globInho = 1;//0
	globHomo = 4;//0
	torqLoad = 20.0;//N*m, u = T*l/G/I_p * R = 1.159111630361142e-06
	coloSett = tempColo;
}

long TORSION::COOR_AVER(const std::vector<COOR> &inpuCoor, COOR &outpCoor){
	COOR cyliCoor(0.0, 0.0, 0.0);
	long tempFlag_0 = 0;
	long tempFlag_1 = 0;
	for(long tk = 0; tk < inpuCoor.size(); tk ++){
		cyliCoor[0] += sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + 
			inpuCoor[tk][1] * inpuCoor[tk][1]
		);
		double tempAngl = atan2(inpuCoor[tk][1], inpuCoor[tk][0]);
		cyliCoor[1] += tempAngl;
		cyliCoor[2] += inpuCoor[tk][2];
		if(tempAngl > PI / 2.0){
			tempFlag_0 ++;
		}
		if(tempAngl < - PI / 2.0){
			tempFlag_1 ++;
		}
	}
	if(tempFlag_0 > 0 && tempFlag_1 > 0){
		cyliCoor[1] += tempFlag_1 * (PI * 2.0);
	}
	outpCoor.resize(3);
	outpCoor[0] = cyliCoor[0] / inpuCoor.size() * cos(cyliCoor[1] / inpuCoor.size());
	outpCoor[1] = cyliCoor[0] / inpuCoor.size() * sin(cyliCoor[1] / inpuCoor.size());
	outpCoor[2] = cyliCoor[2] / inpuCoor.size();
	return 1;
}

long TORSION::SUBR_PLSU(long tg, long ei, std::map<std::vector<long>,COOR> &planSurf){
	std::vector<long> inpuNode(8);
	std::vector<COOR> inpuCoor(8);
	for(long tk = 0; tk < 8; tk ++){
		inpuNode[tk] = multGrid[tg].elemVect[ei].cornNode[tk];
		auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
		inpuCoor[tk] = iterNoco->second;
	}
	COOR outpCoor;
	COOR_AVER(inpuCoor, outpCoor);
	sort(inpuNode.begin(), inpuNode.end());
	planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
		inpuNode, outpCoor
	));
	//
	for(long tj = 0; tj < hexaLine.size(); tj ++){
		std::vector<long> inpuNode(hexaLine[tj].size());
		std::vector<COOR> inpuCoor(hexaLine[tj].size());
		for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
			inpuNode[tk] = multGrid[tg].elemVect[ei].cornNode[hexaLine[tj][tk]];
			auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
			inpuCoor[tk] = iterNoco->second;
		}
		COOR outpCoor;
		COOR_AVER(inpuCoor, outpCoor);
		sort(inpuNode.begin(), inpuNode.end());
		planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
			inpuNode, outpCoor
		));
	}
	//
	for(long tj = 0; tj < hexaFace.size(); tj ++){
		std::vector<long> inpuNode(hexaFace[tj].size());
		std::vector<COOR> inpuCoor(hexaFace[tj].size());
		for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
			inpuNode[tk] = multGrid[tg].elemVect[ei].cornNode[hexaFace[tj][tk]];
			auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
			inpuCoor[tk] = iterNoco->second;
		}
		COOR outpCoor;
		COOR_AVER(inpuCoor, outpCoor);
		sort(inpuNode.begin(), inpuNode.end());
		planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
			inpuNode, outpCoor
		));
	}
	return 1;
}

long TORSION::SUBR_COLO(long tg){
	//displacement constraint must at first
	for(const auto& iterNoco : multGrid[tg].nodeCoor){
		if((iterNoco.second)[2] <= 1.0E-10){
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
		}
	}
	//load must at next
	double torqScal = 2.0 * torqLoad / (pow(outeRadi, 4.0) - pow(inneRadi, 4.0)) / PI;
	for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
		if(multGrid[tg].elemVect[ti].children.size() > 0){
			continue;
		}
		for(long tj = 0; tj < hexaFace.size(); tj ++){
			std::vector<long> inpuNode(hexaFace[tj].size());
			std::vector<Eigen::Vector3d> elemCoor(4);
			Eigen::Matrix<double,12,1> elemCoor_1;
			bool tempFlag = true;
			for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
				inpuNode[tk] = multGrid[tg].elemVect[ti].cornNode[hexaFace[tj][tk]];
				auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
				elemCoor[tk] << (iterNoco->second)[0], 
					(iterNoco->second)[1], (iterNoco->second)[2];
				elemCoor_1.block(3*tk+0,0,3,1) << (iterNoco->second)[0], 
					(iterNoco->second)[1], (iterNoco->second)[2];
				if((iterNoco->second)[2] < axiaLeng - 1.0E-10){
					tempFlag = false;
					break;
				}
			}
			if(tempFlag == false){
				continue;
			}
			Eigen::Matrix<double,12,1> tempForc = Eigen::MatrixXd::Zero(12,1);
			for(long tk = 0; tk < biliQuad.numbNgip; tk ++){
				double weigFact;
				BIQU_TRJA(biliQuad.nacoNgip[tk], elemCoor, weigFact);
				Eigen::VectorXd tempCoor = biliQuad.shfuNgip[tk] * elemCoor_1;
				double tempAngl = atan2(tempCoor(1), tempCoor(0)) + PI / 2.0;
				double tempAmpl = torqScal * sqrt(pow(tempCoor(0), 2.0) + pow(tempCoor(1), 2.0));
				Eigen::Vector3d tempPres;
				tempPres << tempAmpl * cos(tempAngl), tempAmpl * sin(tempAngl), 0.0;
				tempForc += biliQuad.niwfNgip[tk] 
					* biliQuad.shfuNgip[tk].transpose() * tempPres * weigFact;
			}
			for(long tk = 0; tk < 4; tk ++){
				multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 0, tempForc(3 * tk + 0));
				multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 1, tempForc(3 * tk + 1));
				multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 2, tempForc(3 * tk + 2));
			}
		}
	}
	//
	multGrid[tg].coupReps = -1;
	return 1;
}

long TORSION::ESTA_SURF(){
	std::vector<long> diviNumb_es = {
		diviNumb[0] * (1 << (globHomo)), 
		diviNumb[1] * (1 << (globHomo)), 
		diviNumb[2] * (1 << (globInho + globHomo))
	};
	std::vector<long> diviReal = {
		diviNumb_es[0] / domaNumb[0], 
		diviNumb_es[1] / domaNumb[1], 
		diviNumb_es[2] / domaNumb[2]
	};
	VECT_RESI(blocSurf, domaNumb[0] * domaNumb[1] * domaNumb[2], 6);
	#pragma omp parallel for
	for(long tg = 0; tg < blocSurf.size(); tg ++){
		VECT_RESI(blocSurf[tg][0].indiPoin, diviReal[1] + 1, diviReal[2] + 1);
		VECT_RESI(blocSurf[tg][1].indiPoin, diviReal[1] + 1, diviReal[2] + 1);
		VECT_RESI(blocSurf[tg][2].indiPoin, diviReal[2] + 1, diviReal[0] + 1);
		VECT_RESI(blocSurf[tg][3].indiPoin, diviReal[2] + 1, diviReal[0] + 1);
		VECT_RESI(blocSurf[tg][4].indiPoin, diviReal[0] + 1, diviReal[1] + 1);
		VECT_RESI(blocSurf[tg][5].indiPoin, diviReal[0] + 1, diviReal[1] + 1);
		long tg_0 = tg / (domaNumb[1] * domaNumb[2]);
		long tg_1 = (tg % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
		long tg_2 = (tg % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
		//nodes
		for(long ti = 0; ti <= diviReal[0]; ti ++){
			long ti_real = tg_0 * diviReal[0] + ti;
			double radi_ti = inneRadi + (outeRadi - inneRadi) 
				/ (double)(diviNumb_es[0]) * ti_real;
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				long tj_real = tg_1 * diviReal[1] + tj;
				double angl_tj = 0.0 + (2.0 * PI - 0.0) / (double)(diviNumb_es[1]) * tj_real;
				for(long tk =  0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					COOR tempCoor(radi_ti * cos(angl_tj), radi_ti * sin(angl_tj), 
						axiaLeng / diviNumb_es[2] * (double)tk_real
					);
					if(ti == 0){
						blocSurf[tg][0].INSERT(tj, tk, tempCoor);
					}
					if(ti == diviReal[0]){
						blocSurf[tg][1].INSERT(tj, tk, tempCoor);
					}
					if(tj == 0){
						blocSurf[tg][2].INSERT(tk, ti, tempCoor);
					}
					if(tj == diviReal[1]){
						blocSurf[tg][3].INSERT(tk, ti, tempCoor);
					}
					if(tk == 0){
						blocSurf[tg][4].INSERT(ti, tj, tempCoor);
					}
					if(tk == diviReal[2]){
						blocSurf[tg][5].INSERT(ti, tj, tempCoor);
					}
				}
			}
		}
	}
	return 1;
}

long TORSION::MESH_NODD(){
	multGrid.resize(1);
	long tg = 0;
	//nodes
	VECTOR3L blocNode;
	blocNode.resize(diviNumb[0] + 1);
	for(long ti = 0; ti <= diviNumb[0]; ti ++){
		blocNode[ti].resize(diviNumb[1] + 1);
		double radi_ti = inneRadi + (outeRadi - inneRadi) / (double)(diviNumb[0]) * ti;
		for(long tj = 0; tj <= diviNumb[1]; tj ++){
			blocNode[ti][tj].resize(diviNumb[2] + 1);
			double angl_tj = 0.0 + (2.0 * PI - 0.0) / (double)(diviNumb[1]) * tj;
			for(long tk =  0; tk <= diviNumb[2]; tk ++){
				COOR tempCoor(radi_ti * cos(angl_tj), radi_ti * sin(angl_tj), 
					axiaLeng / diviNumb[2] * (double)tk
				);
				blocNode[ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
			}
		}
	}
	//elements
	for(long ti = 0; ti < blocNode.size() - 1; ti ++){
		for(long tj = 0; tj < blocNode[ti].size() - 1; tj ++){
			for(long tk = 0; tk < blocNode[ti][tj].size() - 1; tk ++){
				TREE_ELEM tempElem;
				tempElem.parent = -1;
				tempElem.cornNode[0] = blocNode[ti][tj][tk];
				tempElem.cornNode[1] = blocNode[ti + 1][tj][tk];
				tempElem.cornNode[2] = blocNode[ti + 1][tj + 1][tk];
				tempElem.cornNode[3] = blocNode[ti][tj + 1][tk];
				tempElem.cornNode[4] = blocNode[ti][tj][tk + 1];
				tempElem.cornNode[5] = blocNode[ti + 1][tj][tk + 1];
				tempElem.cornNode[6] = blocNode[ti + 1][tj + 1][tk + 1];
				tempElem.cornNode[7] = blocNode[ti][tj + 1][tk + 1];
				tempElem.level = 0;
				tempElem.refiPatt = 7;
				tempElem.children.resize(0);
				long elemNumb = multGrid[tg].ADD_ELEMENT(tempElem);
			}
		}
	}
	//global refinement
	std::set<long> spliElem;
	std::map<long,std::set<long>> spliFlag;
	std::map<std::vector<long>,COOR> planSurf;
	//global refinement
	for(long tr = 0; tr < globInho + globHomo; tr ++){
		spliElem.clear();
		for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
			if(multGrid[tg].elemVect[ti].children.size() > 0){
				continue;
			}
			spliElem.insert(ti);
			if(tr == 0){
				multGrid[tg].elemVect[ti].refiPatt = 6;
			}
			else if(tr < globInho){
				multGrid[tg].elemVect[ti].refiPatt = 2;
			}
			else{
				multGrid[tg].elemVect[ti].refiPatt = 0;
			}
		}
		planSurf.clear();
		for(const auto &iterSpel : spliElem){
			SUBR_PLSU(tg, iterSpel, planSurf);
		}
		multGrid[tg].REFINE(spliElem, spliFlag, planSurf);
	}
	//
	multGrid[tg].OUTPUT_ELEMENT(tg);
	//
	SUBR_COLO(tg);
	return 1;
}

long TORSION::MESH_DD(){
	multGrid.resize(domaNumb[0] * domaNumb[1] * domaNumb[2]);
	std::vector<long> diviReal = {diviNumb[0] / domaNumb[0], 
		diviNumb[1] / domaNumb[1], diviNumb[2] / domaNumb[2]
	};
	#pragma omp parallel for
	for(long tg = 0; tg < multGrid.size(); tg ++){
		long tg_0 = tg / (domaNumb[1] * domaNumb[2]);
		long tg_1 = (tg % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
		long tg_2 = (tg % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
		//nodes
		VECTOR3L blocNode;
		blocNode.resize(diviReal[0] + 1);
		for(long ti = 0; ti <= diviReal[0]; ti ++){
			blocNode[ti].resize(diviReal[1] + 1);
			long ti_real = tg_0 * diviReal[0] + ti;
			double radi_ti = inneRadi + (outeRadi - inneRadi) / (double)(diviNumb[0]) * ti_real;
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				blocNode[ti][tj].resize(diviReal[2] + 1);
				long tj_real = tg_1 * diviReal[1] + tj;
				double angl_tj = 0.0 + (2.0 * PI - 0.0) / (double)(diviNumb[1]) * tj_real;
				for(long tk =  0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					COOR tempCoor(radi_ti * cos(angl_tj), radi_ti * sin(angl_tj), 
						axiaLeng / diviNumb[2] * (double)tk_real
					);
					blocNode[ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		//elements
		for(long ti = 0; ti < blocNode.size() - 1; ti ++){
			for(long tj = 0; tj < blocNode[ti].size() - 1; tj ++){
				for(long tk = 0; tk < blocNode[ti][tj].size() - 1; tk ++){
					TREE_ELEM tempElem;
					tempElem.parent = -1;
					tempElem.cornNode[0] = blocNode[ti][tj][tk];
					tempElem.cornNode[1] = blocNode[ti + 1][tj][tk];
					tempElem.cornNode[2] = blocNode[ti + 1][tj + 1][tk];
					tempElem.cornNode[3] = blocNode[ti][tj + 1][tk];
					tempElem.cornNode[4] = blocNode[ti][tj][tk + 1];
					tempElem.cornNode[5] = blocNode[ti + 1][tj][tk + 1];
					tempElem.cornNode[6] = blocNode[ti + 1][tj + 1][tk + 1];
					tempElem.cornNode[7] = blocNode[ti][tj + 1][tk + 1];
					tempElem.level = 0;
					tempElem.refiPatt = 7;
					tempElem.children.resize(0);
					long elemNumb = multGrid[tg].ADD_ELEMENT(tempElem);
				}
			}
		}
		//global refinement
		std::set<long> spliElem;
		std::map<long,std::set<long>> spliFlag;
		std::map<std::vector<long>,COOR> planSurf;
		for(long tr = 0; tr < globInho + globHomo; tr ++){
			spliElem.clear();
			for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
				if(multGrid[tg].elemVect[ti].children.size() > 0){
					continue;
				}
				spliElem.insert(ti);
				if(tr < globInho){
					multGrid[tg].elemVect[ti].refiPatt = 6;
				}
				else{
					multGrid[tg].elemVect[ti].refiPatt = 0;
				}
			}
			planSurf.clear();
			for(const auto &iterSpel : spliElem){
				SUBR_PLSU(tg, iterSpel, planSurf);
			}
			multGrid[tg].REFINE(spliElem, spliFlag, planSurf);
		}
		//
		multGrid[tg].OUTPUT_ELEMENT(tg);
		//
		SUBR_COLO(tg);
	}
	return 1;
}

long TORSION::SOLVE(long appsCont, long noddPrec){
	if(coloSett == 0){
		MESH_NODD();
		SOLVE_NODD(noddPrec);
	}
	else{
		ESTA_SURF();
		MESH_DD();
		SOLVE_DD(appsCont);
	}
	return 1;
}

long TORSION::SOLVE_NODD(long noddPrec){
	long tg = 0;
	multGrid[tg].TRANSFER();
	multGrid[tg].STIF_MATR();
	Eigen::VectorXd resuSolu;
	multGrid[tg].CONSTRAINT(1);
	multGrid[tg].mgpi.CG_SOLV(noddPrec, multGrid[tg].consForc, resuSolu);
	Eigen::VectorXd outpDisp;
	multGrid[tg].OUTP_SUB1(resuSolu, outpDisp);
	multGrid[tg].OUTP_SUB2(outpDisp, tg);
	multGrid[tg].STRESS_RECOVERY(outpDisp, tg);
	return 1;
}

long TORSION::SOLVE_DD(long appsCont){
	//
	std::vector<long> xyzN = {
		(domaNumb[0] - 1) * domaNumb[1] * domaNumb[2], 
		(domaNumb[1] - 0) * domaNumb[0] * domaNumb[2], 
		(domaNumb[2] - 1) * domaNumb[0] * domaNumb[1]
	};
	//
	searCont.resize(xyzN[0] + xyzN[1] + xyzN[2]);
	penaFact_n.resize(xyzN[0] + xyzN[1] + xyzN[2]);
	penaFact_f.resize(xyzN[0] + xyzN[1] + xyzN[2]);
	fricCoef.resize(xyzN[0] + xyzN[1] + xyzN[2]);
	contBody.resize(xyzN[0] + xyzN[1] + xyzN[2]);
	for(long tg_0 = 0; tg_0 < domaNumb[0]; tg_0 ++){
		for(long tg_1 = 0; tg_1 < domaNumb[1]; tg_1 ++){
			for(long tg_2 = 0; tg_2 < domaNumb[2]; tg_2 ++){
				long tg_m = tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2;
				if(tg_0 <= domaNumb[0] - 2){
					long tg_s = tg_m + domaNumb[1] * domaNumb[2];
					long ts = tg_m;
					contBody[ts] = {tg_m, tg_s};
				}
				if(tg_1 <= domaNumb[1] - 1){
					long tg_s = (tg_m + domaNumb[2]) % multGrid.size();
					long ts = xyzN[0] 
						+ tg_1 * domaNumb[0] * domaNumb[2] + tg_0 * domaNumb[2] + tg_2;
					contBody[ts] = {tg_m, tg_s};
				}
				if(tg_2 <= domaNumb[2] - 2){
					long tg_s = tg_m + 1;
					long ts = xyzN[0] + xyzN[1] 
						+ tg_2 * domaNumb[0] * domaNumb[1] + tg_0 * domaNumb[1] + tg_1;
					contBody[ts] = {tg_m, tg_s};
				}
			}
		}
	}
	for(long ts = 0; ts < searCont.size(); ts ++){
		searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
		searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
		penaFact_n[ts] = 210.0E9 * 1000.0;
		penaFact_f[ts] = 210.0E9 * 1000.0;
		fricCoef[ts] = -1.0;
	}
	//
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		const CURVEDS *mastSurf, *slavSurf;
		if(ts < xyzN[0]){
			mastSurf = &(blocSurf[contBody[ts][0]][1]);
			slavSurf = &(blocSurf[contBody[ts][1]][0]);
		}
		else if(ts < xyzN[0] + xyzN[1]){
			mastSurf = &(blocSurf[contBody[ts][0]][3]);
			slavSurf = &(blocSurf[contBody[ts][1]][2]);
		}
		else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
			mastSurf = &(blocSurf[contBody[ts][0]][5]);
			slavSurf = &(blocSurf[contBody[ts][1]][4]);
		}
		//
		EFACE_SURFACE iterEfsu_0(searCont[ts].mastGrid, mastSurf);
		while(iterEfsu_0.INCREMENT() == 1){
			searCont[ts].mastSegm.emplace_back(iterEfsu_0.currNode);
		}
		EFACE_SURFACE iterEfsu_1(searCont[ts].slavGrid, slavSurf);
		while(iterEfsu_1.INCREMENT() == 1){
			searCont[ts].slavSegm.emplace_back(iterEfsu_1.currNode);
		}
		searCont[ts].OUTPUT_COSE(ts);
		//
		VECTOR2D mastCoor(2);
		for(long ti = 0; ti < searCont[ts].mastSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].mastSegm[ti][tj];
				COOR tempCoor = multGrid[contBody[ts][0]].nodeCoor[node_tj];
				double radi_tj = sqrt(pow(tempCoor[0], 2.0) + pow(tempCoor[1], 2.0));
				double angl_tj = atan2(tempCoor[1], tempCoor[0]);
				if(angl_tj < 0.0){
					angl_tj += 2.0 * PI;
				}
				if(tempCoor[0] > 0.0 && abs(tempCoor[1]) < 1.0E-10){
					angl_tj = 0.0;
				}
				if(ts < xyzN[0]){
					tempXico += angl_tj;
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += radi_tj;
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
					tempXico += radi_tj;
					tempEtac += angl_tj;
				}
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		if(ts < xyzN[0]){
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[1] / domaNumb[1] * (1 << (globHomo)), 
				diviNumb[2] / domaNumb[2] * (1 << (globInho + globHomo))
			});
		}
		else if(ts < xyzN[0] + xyzN[1]){
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[0] / domaNumb[0] * (1 << (globHomo)), 
				diviNumb[2] / domaNumb[2] * (1 << (globInho + globHomo))
			});
		}
		else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[0] / domaNumb[0] * (1 << (globHomo)), 
				diviNumb[1] / domaNumb[1] * (1 << (globHomo))
			});
		}
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				COOR tempCoor = multGrid[contBody[ts][1]].nodeCoor[node_tj];
				double radi_tj = sqrt(pow(tempCoor[0], 2.0) + pow(tempCoor[1], 2.0));
				double angl_tj = atan2(tempCoor[1], tempCoor[0]);
				if(angl_tj < 0.0){
					angl_tj += 2.0 * PI;
				}
				if(tempCoor[0] > 0.0 && abs(tempCoor[1]) < 1.0E-10){
					angl_tj = 0.0;
				}
				if(ts < xyzN[0]){
					tempXico += angl_tj;
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += radi_tj;
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
					tempXico += radi_tj;
					tempEtac += angl_tj;
				}
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
	}
	//
	ESTABLISH();
	if(appsCont == 0){
		APPS();
	}
	else{
		//initialization
		// for(long tv = 0; tv < multGrid.size(); tv ++){//for MONITOR
			// double tempDisp = 0.0 + (1.159E-6 / domaNumb[2]) * (tv % domaNumb[2]);
			// for(const auto &iterMgnc : multGrid[tv].nodeCoor){
				// COOR tempCoor = iterMgnc.second;
				// double tempAngl = atan2(tempCoor[1], tempCoor[0]) + PI / 2.0;
				// (resuDisp[tv])(3 * iterMgnc.first + 0) = tempDisp * cos(tempAngl);
				// (resuDisp[tv])(3 * iterMgnc.first + 1) = tempDisp * sin(tempAngl);
			// }
		// }
		// #pragma omp parallel for
		// for(long tv = 0; tv < multGrid.size(); tv ++){
			// multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
		// }
		// long test; std::cin >> test;
		// for(long ts = 0; ts < searCont.size(); ts ++){//necessary
			// for(long tv = 0; tv < 2; tv ++){
				// long tempBody = contBody[ts][tv];
				// double tempDisp = 0.0 + (1.159E-6 / domaNumb[2]) * (tempBody % domaNumb[2]);
				// for(const auto iterNoco : nodeCont[ts][tv]){
					// COOR tempCoor = multGrid[tempBody].nodeCoor[iterNoco.first];
					// double tempAngl = atan2(tempCoor[1], tempCoor[0]) + PI / 2.0;
					// (inteAuxi[ts][tv])(3 * iterNoco.second + 0) = tempDisp * cos(tempAngl);
					// (inteAuxi[ts][tv])(3 * iterNoco.second + 1) = tempDisp * sin(tempAngl);
				// }
			// }
		// }
		CONTACT_ANALYSIS();
		//
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
			multGrid[tv].STRESS_RECOVERY(resuDisp[tv], tv);
		}
	}
	return 1;
}
