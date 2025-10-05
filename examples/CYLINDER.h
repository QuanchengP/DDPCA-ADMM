
#include "../MCONTACT.h"

constexpr double struScalFact = 1.0;//structure scale factor

class CYLINDER: public MCONTACT{
public:
	long copyNumb;//axial number of domains for each cylinder
	//0 - lower cylinder, 1/2 - middle cylinder, 3 - upper cylinder
	std::vector<double> radi;//radius of cylinder
	std::vector<double> leng;//length of cylinder
	std::vector<double> diviAngl;//relative to - PI / 2.0
	std::vector<std::vector<Eigen::Vector2d>> diviPoin;//auxiliary point for meshing
	VECTOR2L diviNumb;//number of elements along different directions
	long globInho;//level of inhomogeneous global refinement
	long globHomo;//level of homogeneous global refinement
	long locaLeve;//level of local refinement
	VECTOR2L buckNumb;//number for bucket search
	double bandWidt;//predicted contact band width > real value
	double loadInte;//load intensity
	//cylindrical surface: 0 - lower cylinder, 1/2 - middle cylinder, 3 - upper cylinder
	std::vector<CURVEDS> cyliSurf;
	std::vector<CURVEDS> cyliSurf_1;//left
	std::vector<CURVEDS> cyliSurf_2;//right
	std::vector<CURVEDS> inteSurf;//interface between two middle cylinders
	//
	CYLINDER();
	long COOR_TRAN_0(COOR &tempCoor);//coordinate transfer for lower cylinder
	long COOR_TRAN_2(COOR &tempCoor);//coordinate transfer for middle cylinder
	long ESTA_SURF();//establish surface
	long MESH();//generate mesh
	long SOLVE(long appsCont = 1);//solve for displacement and contact pressure
	double charFact = 25.0;
};

CYLINDER::CYLINDER(){
	//
	mkdir("Cylinder", 0777);
	outpDire = "Cylinder/";
	//
	copyNumb = 16;
	muscSett = (1 << 0);
	doleMcsc.assign(4 * 2 * copyNumb, 0);
	//
	radi = {0.02 * struScalFact, 0.022 * struScalFact, 0.022 * struScalFact, 0.02 * struScalFact};
	diviAngl = {- 5.0 / 8.0 * PI, - 3.0 / 8.0 * PI};
	diviPoin.resize(4);
	for(long tg = 0; tg < 4; tg ++){
		diviPoin[tg].resize(3);
		diviPoin[tg][0] << - radi[tg] / 3.0, 0.0;
		diviPoin[tg][1] << - radi[tg] / 5.0, - radi[tg] / 2.0;
		diviPoin[tg][2] << radi[tg] / 5.0, - radi[tg] / 2.0;
	}
	//globInho + globHomo must be larger than 1, 
	//so that diviNumb[][1] becomes an even number for load
	globInho = 3;
	globHomo = 0;
	locaLeve = 7;
	//
	bandWidt = 100.0E-6 * struScalFact;//half contact band width
	loadInte = -50.0E3 * struScalFact;
}

long CYLINDER::COOR_TRAN_0(COOR &tempCoor){
	tempCoor[0] = - tempCoor[0];
	tempCoor[1] = - tempCoor[1] - radi[0] - radi[1] - radi[2] - radi[3];
	return 1;
}

long CYLINDER::COOR_TRAN_2(COOR &tempCoor){
	tempCoor[0] = - tempCoor[0];
	tempCoor[1] = - tempCoor[1] - radi[2] - radi[3];
	return 1;
}

long CYLINDER::ESTA_SURF(){
	cyliSurf.resize(4);
	inteSurf.resize(3);
	cyliSurf_1.resize(4);
	cyliSurf_2.resize(4);
	//
	for(long tg = 0; tg < 4; tg ++){
		long totaDivi_0 = diviNumb[tg][1] * (1 << (globInho + globHomo + locaLeve));
		long totaDivi_1 = diviNumb[tg][3] * (1 << (globHomo + locaLeve));
		cyliSurf[tg].indiPoin.resize(totaDivi_0 + 1);
		for(long ti = 0; ti <= totaDivi_0; ti ++){
			cyliSurf[tg].indiPoin[ti].resize(totaDivi_1 + 1);
			double angl_ti = diviAngl[0] + (diviAngl[1] - diviAngl[0]) / (double)totaDivi_0 * ti;
			for(long tj = 0; tj <= totaDivi_1; tj ++){
				COOR tempCoor(radi[tg] * cos(angl_ti), radi[tg] * sin(angl_ti), 
					leng[tg] / (double)totaDivi_1 * tj
				);
				if(tg == 0){
					COOR_TRAN_0(tempCoor);
				}
				else if(tg == 1){
					tempCoor[1] = tempCoor[1] - radi[0] - radi[1];
				}
				else if(tg == 2){
					COOR_TRAN_2(tempCoor);
				}
				cyliSurf[tg].INSERT(ti, tj, tempCoor);
			}
		}
	}
	//
	for(long tg = 0; tg < 4; tg ++){
		long totaDivi_0 = diviNumb[tg][0] * (1 << (globInho + globHomo));
		long totaDivi_1 = diviNumb[tg][3] * (1 << (globHomo));
		cyliSurf_1[tg].indiPoin.resize(totaDivi_0 + 1);
		cyliSurf_2[tg].indiPoin.resize(totaDivi_0 + 1);
		for(long ti = 0; ti <= totaDivi_0; ti ++){
			cyliSurf_1[tg].indiPoin[ti].resize(totaDivi_1 + 1);
			cyliSurf_2[tg].indiPoin[totaDivi_0 - ti].resize(totaDivi_1 + 1);
			double angl_ti = - PI + (diviAngl[0] + PI) / (double)totaDivi_0 * ti;
			for(long tj = 0; tj <= totaDivi_1; tj ++){
				COOR tempCoor(radi[tg] * cos(angl_ti), radi[tg] * sin(angl_ti), 
					leng[tg] / (double)totaDivi_1 * tj
				);
				if(tg == 0){
					COOR_TRAN_0(tempCoor);
				}
				else if(tg == 1){
					tempCoor[1] = tempCoor[1] - radi[0] - radi[1];
				}
				else if(tg == 2){
					COOR_TRAN_2(tempCoor);
				}
				cyliSurf_1[tg].INSERT(ti, tj, tempCoor);
				tempCoor[0] = - tempCoor[0];
				cyliSurf_2[tg].INSERT(totaDivi_0 - ti, tj, tempCoor);
			}
		}
	}
	//
	for(long ta = 0; ta < 3; ta ++){
		long totaDivi_0;
		long totaDivi_1 = diviNumb[1][3] * (1 << globHomo);
		if(ta == 1){
			totaDivi_0 = diviNumb[1][1] * (1 << (globInho + globHomo));
		}
		else{
			totaDivi_0 = diviNumb[1][2] * (1 << (globInho + globHomo));
		}
		inteSurf[ta].indiPoin.resize(totaDivi_0 + 1);
		for(long ti = 0; ti <= totaDivi_0; ti ++){
			inteSurf[ta].indiPoin[ti].resize(totaDivi_1 + 1);
			double xcoo_ti;
			if(ta == 0){
				xcoo_ti = - radi[1] + (diviPoin[1][0](0) + radi[1]) / (double)totaDivi_0 * ti;
			}
			else if(ta == 1){
				xcoo_ti = diviPoin[1][0](0) 
					+ (- diviPoin[1][0](0) - diviPoin[1][0](0)) / (double)totaDivi_0 * ti;
			}
			else if(ta == 2){
				xcoo_ti = - diviPoin[1][0](0) 
					+ (radi[1] + diviPoin[1][0](0)) / (double)totaDivi_0 * ti;
			}
			for(long tj = 0; tj <= totaDivi_1; tj ++){
				COOR tempCoor(xcoo_ti, - radi[3] - radi[2], 
					leng[1] / (double)totaDivi_1 * tj
				);
				inteSurf[ta].INSERT(ti, tj, tempCoor);
			}
		}
	}
	return 1;
}

long CYLINDER::MESH(){
	doleMcsc.assign(4 * 2 * copyNumb, 2);
	diviNumb = {{2, 2, 1, 16 / copyNumb}, {2, 2, 1, 16 / copyNumb}, 
		{2, 2, 1, 16 / copyNumb}, {2, 2, 1, 16 / copyNumb}
	};
	//leng must be equal
	leng = {0.02 * struScalFact / copyNumb, 0.02 * struScalFact / copyNumb, 
		0.02 * struScalFact / copyNumb, 0.02 * struScalFact / copyNumb
	};
	buckNumb.resize(11);
	//buckNumb can not be too large, the smaller is the safer (although lower efficiency)
	buckNumb[0] = {8, 
		std::min(diviNumb[0][3], diviNumb[1][3]) * (1 << (globHomo + locaLeve - 1))
	};
	buckNumb[1] = {8, 
		std::min(diviNumb[2][3], diviNumb[3][3]) * (1 << (globHomo + locaLeve - 1))
	};
	buckNumb[2] = {(diviNumb[1][2] + diviNumb[1][1] / 2) * (1 << (globInho + globHomo - 1)), 
		diviNumb[1][3] * (1 << std::max((long)0, globHomo - 1))
	};
	for(long ts = 0; ts < 4; ts ++){
		buckNumb[3 + ts] = {
			(diviNumb[ts][2] + diviNumb[ts][1] / 2) * (1 << (globInho + globHomo - 1)), 
			(diviNumb[ts][0] + diviNumb[ts][1]) * (1 << (globInho + globHomo - 1))
		};
	}
	for(long ts = 0; ts < 4; ts ++){
		buckNumb[7 + ts] = {
			(diviNumb[ts][0] + diviNumb[ts][1]) * (1 << (globInho + globHomo - 1)), 
			1 << (globHomo + locaLeve - 1)//matching meshes
		};
	}
	//
	ESTA_SURF();
	//
	multGrid.resize(4 * 2 * copyNumb);
	#pragma omp parallel for
	for(long tg = 0; tg < 4; tg ++){
		//boundary nodes of block
		std::vector<Eigen::Vector2d> uppeLine_0, uppeLine_1, uppeLine_2, downLine_0, downLine_1;
		uppeLine_0.resize(diviNumb[tg][0] + 1);
		downLine_0.resize(diviNumb[tg][0] + 1);
		uppeLine_1.resize(diviNumb[tg][1] / 2 + 1);
		downLine_1.resize(diviNumb[tg][1] / 2 + 1);
		uppeLine_2.resize(diviNumb[tg][1] / 2 + 1);
		for(long ti = 0; ti <= diviNumb[tg][0]; ti ++){
			uppeLine_0[ti] = (1.0 - (double)ti / diviNumb[tg][0]) * diviPoin[tg][0] 
				+ (double)ti / diviNumb[tg][0] * diviPoin[tg][1];
			double tempAngl = - PI + (diviAngl[0] + PI) / diviNumb[tg][0] * (double)ti;
			downLine_0[ti] << radi[tg] * cos(tempAngl), radi[tg] * sin(tempAngl);
		}
		for(long ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
			uppeLine_1[ti] = (1.0 - (double)ti / diviNumb[tg][1]) * diviPoin[tg][1] 
				+ (double)ti / diviNumb[tg][1] * diviPoin[tg][2];
			double tempAngl = diviAngl[0] 
				+ (diviAngl[1] - diviAngl[0]) / diviNumb[tg][1] * (double)ti;
			downLine_1[ti] << radi[tg] * cos(tempAngl), radi[tg] * sin(tempAngl);
			uppeLine_2[ti](0) = (1.0 - (double)ti / diviNumb[tg][1]) * (-radi[tg] / 3.0)
				+ (double)ti / diviNumb[tg][1] * (radi[tg] / 3.0);
			uppeLine_2[ti](1) = 0.0;
		}
		//node of block 0~2
		VECTOR4L blocNode(3);
		blocNode[0].resize(diviNumb[tg][0] + 1);
		for(long ti = 0; ti <= diviNumb[tg][0]; ti ++){
			blocNode[0][ti].resize(diviNumb[tg][2] + 1);
			for(long tj = 0; tj <= diviNumb[tg][2]; tj ++){
				blocNode[0][ti][tj].resize(diviNumb[tg][3] + 1);
				Eigen::Vector2d poin_ij = 
					(1.0 - (double)tj / diviNumb[tg][2]) * downLine_0[ti] 
					+ (double)tj / diviNumb[tg][2] * uppeLine_0[ti];
				for(long tk =  0; tk <= diviNumb[tg][3]; tk ++){
					COOR tempCoor(poin_ij(0), poin_ij(1), 
						leng[tg] / diviNumb[tg][3] * (double)tk
					);
					if(tg == 0){
						COOR_TRAN_0(tempCoor);
					}
					else if(tg == 1){
						tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
					}
					else if(tg == 2){
						COOR_TRAN_2(tempCoor);
					}
					blocNode[0][ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		blocNode[1].resize(diviNumb[tg][1] / 2 + 1);
		for(long ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
			blocNode[1][ti].resize(diviNumb[tg][2] + 1);
			for(long tj = 0; tj <= diviNumb[tg][2]; tj ++){
				blocNode[1][ti][tj].resize(diviNumb[tg][3] + 1);
				Eigen::Vector2d poin_ij = 
					(1.0 - (double)tj / diviNumb[tg][2]) * downLine_1[ti] 
					+ (double)tj / diviNumb[tg][2] * uppeLine_1[ti];
				for(long tk = 0; tk <= diviNumb[tg][3]; tk ++){
					COOR tempCoor(poin_ij(0), poin_ij(1), 
						leng[tg] / diviNumb[tg][3] * (double)tk
					);
					if(tg == 0){
						COOR_TRAN_0(tempCoor);
					}
					else if(tg == 1){
						tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
					}
					else if(tg == 2){
						COOR_TRAN_2(tempCoor);
					}
					blocNode[1][ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		blocNode[2].resize(diviNumb[tg][1] / 2 + 1);
		for(long ti = 0; ti <= diviNumb[tg][1] / 2; ti ++){
			blocNode[2][ti].resize(diviNumb[tg][0] + 1);
			for(long tj = 0; tj <= diviNumb[tg][0]; tj ++){
				blocNode[2][ti][tj].resize(diviNumb[tg][3] + 1);
				Eigen::Vector2d poin_ij = 
					(1.0 - (double)tj / diviNumb[tg][0]) * uppeLine_1[ti] 
					+ (double)tj / diviNumb[tg][0] * uppeLine_2[ti];
				for(long tk = 0; tk <= diviNumb[tg][3]; tk ++){
					COOR tempCoor(poin_ij(0), poin_ij(1), 
						leng[tg] / diviNumb[tg][3] * (double)tk
					);
					if(tg == 0){
						COOR_TRAN_0(tempCoor);
					}
					else if(tg == 1){
						tempCoor[1] = tempCoor[1] - radi[3] - radi[2];
					}
					else if(tg == 2){
						COOR_TRAN_2(tempCoor);
					}
					blocNode[2][ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		//elements
		for(long ti = 0; ti < blocNode.size(); ti ++){
			for(long tj = 0; tj < blocNode[ti].size() - 1; tj ++){
				for(long tk = 0; tk < blocNode[ti][tj].size() - 1; tk ++){
					for(long tm = 0; tm < blocNode[ti][tj][tk].size() - 1; tm ++){
						TREE_ELEM tempElem;
						tempElem.parent = -1;
						tempElem.cornNode[0] = blocNode[ti][tj][tk][tm];
						tempElem.cornNode[1] = blocNode[ti][tj + 1][tk][tm];
						tempElem.cornNode[2] = blocNode[ti][tj + 1][tk + 1][tm];
						tempElem.cornNode[3] = blocNode[ti][tj][tk + 1][tm];
						tempElem.cornNode[4] = blocNode[ti][tj][tk][tm + 1];
						tempElem.cornNode[5] = blocNode[ti][tj + 1][tk][tm + 1];
						tempElem.cornNode[6] = blocNode[ti][tj + 1][tk + 1][tm + 1];
						tempElem.cornNode[7] = blocNode[ti][tj][tk + 1][tm + 1];
						tempElem.level = 0;
						tempElem.refiPatt = 7;
						tempElem.children.resize(0);
						long elemNumb = multGrid[tg].ADD_ELEMENT(tempElem);
					}
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
				multGrid[tg].elemVect[ti].refiPatt = (tr < globInho) ? 1 : 0;
			}
			planSurf.clear();
			cyliSurf[tg].REFINE(spliElem, multGrid[tg].elemVect, 
				multGrid[tg].nodeCoor, planSurf
			);
			cyliSurf_1[tg].REFINE(spliElem, multGrid[tg].elemVect, 
				multGrid[tg].nodeCoor, planSurf
			);
			cyliSurf_2[tg].REFINE(spliElem, multGrid[tg].elemVect, 
				multGrid[tg].nodeCoor, planSurf
			);
			if(tg == 1 || tg == 2){
				for(long ti = 0; ti < inteSurf.size(); ti ++){
					inteSurf[ti].REFINE(spliElem, multGrid[tg].elemVect, 
						multGrid[tg].nodeCoor, planSurf
					);
				}
			}
			multGrid[tg].REFINE(spliElem, spliFlag, planSurf);
		}
		//local refinement
		std::vector<double> refeRadi = {- radi[3] - radi[2] - radi[1], 
			- radi[3] - radi[2] - radi[1], - radi[3], - radi[3]
		};
		for(long tr = 0; tr < locaLeve; tr ++){
			spliElem.clear();
			//***********************************************************************************
			for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
				if(multGrid[tg].elemVect[ti].children.size() > 0){
					continue;
				}
				const auto &tempElem = multGrid[tg].elemVect[ti];
				bool tempFlag = false;
				for(long tj = 0; tj < 8; tj ++){
					auto iterNoco = multGrid[tg].nodeCoor.find(tempElem.cornNode[tj]);
					if(tg == 0 && abs((iterNoco->second)[0]) <= bandWidt
						&& (iterNoco->second)[1] >= refeRadi[tg] - 2.0 * bandWidt){
						tempFlag = true;
						break;
					}
					else if(tg == 1 && abs((iterNoco->second)[0]) <= bandWidt
						&& (iterNoco->second)[1] <= refeRadi[tg] + 2.0 * bandWidt){
						tempFlag = true;
						break;
					}
					else if(tg == 2 && abs((iterNoco->second)[0]) <= bandWidt
						&& (iterNoco->second)[1] >= refeRadi[tg] - 2.0 * bandWidt){
						tempFlag = true;
						break;
					}
					else if(tg == 3 && abs((iterNoco->second)[0]) <= bandWidt
						&& (iterNoco->second)[1] <= refeRadi[tg] + 2.0 * bandWidt){
						tempFlag = true;
						break;
					}
				}
				if(tempFlag == true){
					spliElem.insert(ti);
					multGrid[tg].elemVect[ti].refiPatt = 0;
				}
			}
			//***********************************************************************************
			// EFACE_SURFACE iterEfsu(&(multGrid[tg]), &(cyliSurf[tg]));
			// while(iterEfsu.INCREMENT() == 1){
				// std::vector<long> tempNode = iterEfsu.currNode;
				// bool tempFlag = false;
				// for(long ti = 0; ti < 4; ti ++){
					// auto iterNoco = multGrid[tg].nodeCoor.find(tempNode[ti]);
					// if(abs((iterNoco->second)[0]) <= bandWidt){
						// tempFlag = true;
						// break;
					// }
				// }
				// //no need: || multGrid[tg].elemVect[iterEfsu.ti].children.size() > 0
				// if(tempFlag == true){
					// spliElem.insert(iterEfsu.ti);
					// multGrid[tg].elemVect[iterEfsu.ti].refiPatt = 0;
				// }
			// }
			//***********************************************************************************
			planSurf.clear();
			cyliSurf[tg].REFINE(spliElem, multGrid[tg].elemVect, 
				multGrid[tg].nodeCoor, planSurf
			);
			multGrid[tg].REFINE(spliElem, spliFlag, planSurf);
		}
		//
		multGrid[tg].OUTPUT_ELEMENT(tg);
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			if(tg == 0 && 
				(iterNoco.second)[1] <= - radi[0] - radi[1] - radi[2] - radi[3] + 1.0E-10){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
			}
			else if((tg == 1 || tg == 2) 
				&& abs((iterNoco.second)[1] + radi[2] + radi[3]) <= 1.0E-10){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
			}
			else if(tg == 3 && (iterNoco.second)[1] >= -1.0E-10){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
			}
		}
		//load must at next
		if(tg == 3){
			double loadIncr = loadInte * leng[tg] / (diviNumb[tg][3] * (1 << (globHomo)));
			for(const auto& iterNoco : multGrid[tg].nodeCoor){
				if((iterNoco.second)[1] >= -1.0E-10 && abs((iterNoco.second)[0]) <= 1.0E-10){
					double tempFact = 0.5;
					if((iterNoco.second)[2] <= 1.0E-10 
						|| (iterNoco.second)[2] >= leng[tg] - 1.0E-10){
						tempFact = 0.25;
					}
					long free_tm = 3 * iterNoco.first + 1;
					multGrid[tg].LOAD_ACCU(free_tm, tempFact * loadIncr);
				}
			}
		}
		//
		multGrid[tg].coupReps = -1;
	}
	//
	#pragma omp parallel for
	for(long tc = 0; tc < 4; tc ++){
		long tg = 4 + tc;
		multGrid[tg].COPY(multGrid[tc]);
		Eigen::Matrix3d tempRota;
		tempRota << -1.0, 0.0, 0.0, 
			0.0, 1.0, 0.0, 
			0.0, 0.0, -1.0;
		Eigen::Vector3d tempTran;
		tempTran << 0.0, 0.0, leng[tc];
		multGrid[tg].RIGI_ROTR(tempRota, tempTran);
		multGrid[tg].OUTPUT_ELEMENT(tg);
	}
	//
	#pragma omp parallel for
	for(long tb = 1; tb < copyNumb; tb ++){
		for(long tc = 0; tc < 8; tc ++){
			long tg = tb * 8 + tc;
			multGrid[tg].COPY(multGrid[tc]);
			Eigen::Matrix3d tempRota;
			tempRota << 1.0, 0.0, 0.0, 
				0.0, 1.0, 0.0, 
				0.0, 0.0, 1.0;
			Eigen::Vector3d tempTran;
			tempTran << 0.0, 0.0, tb * leng[tc % 4];
			multGrid[tg].RIGI_ROTR(tempRota, tempTran);
			multGrid[tg].OUTPUT_ELEMENT(tg);
		}
	}
	return 1;
}

long CYLINDER::SOLVE(long appsCont){
	//
	MESH();
	double charLeng = GET_CHAR_LENG();
	//
	long totaSear = 6 * copyNumb + 8 * (copyNumb - 1) + 4 * copyNumb;
	searCont.resize(totaSear);
	penaFact_n.resize(totaSear);
	penaFact_f.resize(totaSear);
	fricCoef.resize(totaSear);
	contBody.resize(totaSear);
	for(long ta = 0; ta < copyNumb; ta ++){
		contBody[6 * ta + 0] = {ta * 8 + 0, ta * 8 + 5};
		contBody[6 * ta + 1] = {ta * 8 + 4, ta * 8 + 1};
		contBody[6 * ta + 2] = {ta * 8 + 2, ta * 8 + 7};
		contBody[6 * ta + 3] = {ta * 8 + 6, ta * 8 + 3};
		contBody[6 * ta + 4] = {ta * 8 + 5, ta * 8 + 2};
		contBody[6 * ta + 5] = {ta * 8 + 1, ta * 8 + 6};
		//
		for(long ti = 0; ti < 6; ti ++){
			searCont[6 * ta + ti].mastGrid = &(multGrid[contBody[6 * ta + ti][0]]);
			searCont[6 * ta + ti].slavGrid = &(multGrid[contBody[6 * ta + ti][1]]);
			penaFact_n[6 * ta + ti] = 210.0E9 * charFact / charLeng;
			penaFact_f[6 * ta + ti] = 210.0E9 * charFact / charLeng;
			fricCoef[6 * ta + ti] = 0.0;
		}
	}
	for(long ta = 0; ta < copyNumb - 1; ta ++){
		for(long tb = 0; tb < 8; tb ++){
			long ts = 6 * copyNumb + 8 * ta + tb;
			contBody[ts] = {ta * 8 + tb, (ta + 1) * 8 + tb};
			//
			searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
			searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
			penaFact_n[ts] = 210.0E9 * charFact / charLeng;
			penaFact_f[ts] = 210.0E9 * charFact / charLeng;
			fricCoef[ts] = -1.0;
		}
	}
	for(long ta = 0; ta < copyNumb; ta ++){
		for(long tb = 0; tb < 4; tb ++){
			long ts = 6 * copyNumb + 8 * (copyNumb - 1) + 4 * ta + tb;
			contBody[ts] = {ta * 8 + tb, ta * 8 + tb + 4};
			//
			searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
			searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
			penaFact_n[ts] = 210.0E9 * charFact / charLeng;
			penaFact_f[ts] = 210.0E9 * charFact / charLeng;
			fricCoef[ts] = -1.0;
		}
	}
	//
	#pragma omp parallel for
	for(long ts = 0; ts < 6; ts ++){
		long tg_mast = contBody[ts][0];
		long tg_slav = contBody[ts][1];
		//
		if(ts <= 3){
			EFACE_SURFACE iterEfsu_0(&(multGrid[tg_mast]), &(cyliSurf[tg_mast % 4]));
			while(iterEfsu_0.INCREMENT() == 1){
				std::vector<long> tempNode = iterEfsu_0.currNode;
				bool tempFlag = true;
				for(long ti = 0; ti < 4; ti ++){
					auto iterNoco = multGrid[tg_mast].nodeCoor.find(tempNode[ti]);
					if(abs((iterNoco->second)[0]) > bandWidt){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == true){
					searCont[ts].mastSegm.emplace_back(tempNode);
				}
			}
			EFACE_SURFACE iterEfsu_1(&(multGrid[tg_slav]), &(cyliSurf[tg_slav % 4]));
			while(iterEfsu_1.INCREMENT() == 1){
				std::vector<long> tempNode = iterEfsu_1.currNode;
				bool tempFlag = true;
				for(long ti = 0; ti < 4; ti ++){
					auto iterNoco = multGrid[tg_slav].nodeCoor.find(tempNode[ti]);
					if(abs((iterNoco->second)[0]) > bandWidt){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == true){
					searCont[ts].slavSegm.emplace_back(tempNode);
				}
			}
		}
		else{
			for(long ta = 0; ta < inteSurf.size(); ta ++){
				EFACE_SURFACE iterEfsu_0(&(multGrid[tg_mast]), &(inteSurf[ta]));
				while(iterEfsu_0.INCREMENT() == 1){
					std::vector<long> tempNode = iterEfsu_0.currNode;
					searCont[ts].mastSegm.emplace_back(tempNode);
				}
				EFACE_SURFACE iterEfsu_1(&(multGrid[tg_slav]), &(inteSurf[ta]));
				while(iterEfsu_1.INCREMENT() == 1){
					std::vector<long> tempNode = iterEfsu_1.currNode;
					searCont[ts].slavSegm.emplace_back(tempNode);
				}
			}
		}
		searCont[ts].OUTPUT_COSE(ts);
		//
		VECTOR2D mastCoor(2);
		for(long ti = 0; ti < searCont[ts].mastSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].mastSegm[ti][tj];
				tempXico += (multGrid[tg_mast].nodeCoor[node_tj])[0];
				tempEtac += (multGrid[tg_mast].nodeCoor[node_tj])[2];
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, {buckNumb[ts / 2][0], buckNumb[ts / 2][1]});
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				tempXico += (multGrid[tg_slav].nodeCoor[node_tj])[0];
				tempEtac += (multGrid[tg_slav].nodeCoor[node_tj])[2];
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
	}
	#pragma omp parallel for
	for(long ta = 1; ta < copyNumb; ta ++){
		for(long tb = 0; tb < 6; tb ++){
			long ts = 6 * ta + tb;
			searCont[ts].COPY(searCont[tb]);
			double zcoo_ab = ta * leng[0];//!!!!!!!!!!
			for(auto &iterInpo : searCont[ts].intePoin){
				(iterInpo.contPoin)[0](2) += zcoo_ab;
				(iterInpo.contPoin)[1](2) += zcoo_ab;
			}
			searCont[ts].OUTPUT_COSE(ts);
			searCont[ts].OUTPUT_INPO(ts);
		}
	}
	//
	if(copyNumb > 1){
		#pragma omp parallel for
		for(long tb = 0; tb < 8; tb ++){
			long ts = 6 * copyNumb + tb;
			//
			long tg_mast = tb;
			long tg_slav = tb + 8;
			for(long ti = 0; ti < multGrid[tg_mast].elemVect.size(); ti ++){
				if(multGrid[tg_mast].elemVect[ti].children.size() > 0){
					continue;
				}
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					bool tempFlag = true;
					std::vector<long> tempNode(4);
					for(long tk = 0; tk < 4; tk ++){
						tempNode[tk] = multGrid[tg_mast].elemVect[ti].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tg_mast].nodeCoor.find(tempNode[tk]);
						if(abs((iterNoco->second)[2] - leng[tb % 4]) > 1.0E-10){
							tempFlag = false;
							break;
						}
					}
					if(tempFlag == true){
						searCont[ts].mastSegm.emplace_back(tempNode);
					}
				}
			}
			for(long ti = 0; ti < multGrid[tg_slav].elemVect.size(); ti ++){
				if(multGrid[tg_slav].elemVect[ti].children.size() > 0){
					continue;
				}
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					bool tempFlag = true;
					std::vector<long> tempNode(4);
					for(long tk = 0; tk < 4; tk ++){
						tempNode[tk] = multGrid[tg_slav].elemVect[ti].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tg_slav].nodeCoor.find(tempNode[tk]);
						if(abs((iterNoco->second)[2] - leng[tb % 4]) > 1.0E-10){
							tempFlag = false;
							break;
						}
					}
					if(tempFlag == true){
						searCont[ts].slavSegm.emplace_back(tempNode);
					}
				}
			}
			searCont[ts].OUTPUT_COSE(ts);
			//
			VECTOR2D mastCoor(2);
			for(long ti = 0; ti < searCont[ts].mastSegm.size(); ti ++){
				double tempXico = 0.0;
				double tempEtac = 0.0;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = searCont[ts].mastSegm[ti][tj];
					tempXico += (multGrid[tg_mast].nodeCoor[node_tj])[0];
					tempEtac += (multGrid[tg_mast].nodeCoor[node_tj])[1];
				}
				mastCoor[0].emplace_back(tempXico / 4.0);
				mastCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].BUCKET_SORT(mastCoor, {buckNumb[3 + tb%4][0], buckNumb[3 + tb%4][1]});
			VECTOR2D slavCoor(2);
			for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
				double tempXico = 0.0;
				double tempEtac = 0.0;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = searCont[ts].slavSegm[ti][tj];
					tempXico += (multGrid[tg_slav].nodeCoor[node_tj])[0];
					tempEtac += (multGrid[tg_slav].nodeCoor[node_tj])[1];
				}
				slavCoor[0].emplace_back(tempXico / 4.0);
				slavCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].CONTACT_SEARCH(slavCoor);
			searCont[ts].OUTPUT_INPO(ts);
		}
		#pragma omp parallel for
		for(long ta = 1; ta < copyNumb - 1; ta ++){
			for(long tb = 0; tb < 8; tb ++){
				long ts = 6 * copyNumb + 8 * ta + tb;
				searCont[ts].COPY(searCont[6 * copyNumb + tb]);
				double zcoo_ab = ta * leng[tb % 4];
				for(auto &iterInpo : searCont[ts].intePoin){
					(iterInpo.contPoin)[0](2) += zcoo_ab;
					(iterInpo.contPoin)[1](2) += zcoo_ab;
				}
				searCont[ts].OUTPUT_COSE(ts);
				searCont[ts].OUTPUT_INPO(ts);
			}
		}
	}
	#pragma omp parallel for
	for(long tb = 0; tb < 4; tb ++){
		long ts = 6 * copyNumb + 8 * (copyNumb - 1) + tb;
		//
		long tg_mast = tb;
		long tg_slav = tb + 4;
		for(long ti = 0; ti < multGrid[tg_mast].elemVect.size(); ti ++){
			if(multGrid[tg_mast].elemVect[ti].children.size() > 0){
				continue;
			}
			for(long tj = 0; tj < hexaFace.size(); tj ++){
				bool tempFlag = true;
				std::vector<long> tempNode(4);
				for(long tk = 0; tk < 4; tk ++){
					tempNode[tk] = multGrid[tg_mast].elemVect[ti].cornNode[hexaFace[tj][tk]];
					auto iterNoco = multGrid[tg_mast].nodeCoor.find(tempNode[tk]);
					if(abs((iterNoco->second)[0]) > 1.0E-10){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == true){
					searCont[ts].mastSegm.emplace_back(tempNode);
				}
			}
		}
		for(long ti = 0; ti < multGrid[tg_slav].elemVect.size(); ti ++){
			if(multGrid[tg_slav].elemVect[ti].children.size() > 0){
				continue;
			}
			for(long tj = 0; tj < hexaFace.size(); tj ++){
				bool tempFlag = true;
				std::vector<long> tempNode(4);
				for(long tk = 0; tk < 4; tk ++){
					tempNode[tk] = multGrid[tg_slav].elemVect[ti].cornNode[hexaFace[tj][tk]];
					auto iterNoco = multGrid[tg_slav].nodeCoor.find(tempNode[tk]);
					if(abs((iterNoco->second)[0]) > 1.0E-10){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == true){
					searCont[ts].slavSegm.emplace_back(tempNode);
				}
			}
		}
		searCont[ts].OUTPUT_COSE(ts);
		//
		VECTOR2D mastCoor(2);
		for(long ti = 0; ti < searCont[ts].mastSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].mastSegm[ti][tj];
				tempXico += (multGrid[tg_mast].nodeCoor[node_tj])[1];
				tempEtac += (multGrid[tg_mast].nodeCoor[node_tj])[2];
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, {buckNumb[7 + tb][0], buckNumb[7 + tb][1]});
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				tempXico += (multGrid[tg_slav].nodeCoor[node_tj])[1];
				tempEtac += (multGrid[tg_slav].nodeCoor[node_tj])[2];
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
	}
	#pragma omp parallel for
	for(long ta = 1; ta < copyNumb; ta ++){
		for(long tb = 0; tb < 4; tb ++){
			long ts = 6 * copyNumb + 8 * (copyNumb - 1) + ta * 4 + tb;
			searCont[ts].COPY(searCont[6 * copyNumb + 8 * (copyNumb - 1) + tb]);
			double zcoo_ab = ta * leng[tb];
			for(auto &iterInpo : searCont[ts].intePoin){
				(iterInpo.contPoin)[0](2) += zcoo_ab;
				(iterInpo.contPoin)[1](2) += zcoo_ab;
			}
			searCont[ts].OUTPUT_COSE(ts);
			searCont[ts].OUTPUT_INPO(ts);
		}
	}
	//
	if(appsCont <= 1){
		ESTABLISH();
	}
	if(appsCont == 0){
		APPS();
	}
	else{
		if(appsCont == 1){
			CONTACT_ANALYSIS();
		}
		else{
			LAGRANGE(appsCont - 1);
		}
		//
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
			multGrid[tv].STRESS_RECOVERY(resuDisp[tv], tv);
		}
	}
	return 1;
}
