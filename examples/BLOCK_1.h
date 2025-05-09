
#include "../MCONTACT.h"

//non-matching interface with cross corner problem
class BLOCK_1: public MCONTACT{
public:
	std::vector<double> leng;//edge length of each block
	std::vector<long> diviNumb;//number of elements along the edge of each block
	long globLeve;;//multigrid level, total levels of global mesh refinement
	std::vector<long> domaNumb;//control the number of decomposed domains
	//0 - upper surface of bottom block, 1 - bottom surface of middle block
	//2 - upper surface of middle block, 3 - bottom surface of upper block
	//4 - upper surface of upper block
	std::vector<CURVEDS> bmupSurf;
	std::vector<std::vector<CURVEDS>> blocSurf;// each domain has six surfaces
	//
	BLOCK_1();
	long ESTA_SURF();//establish surface
	long MESH();//generate mesh
	Eigen::Vector3d loadPres;//load pressure
	long LOAD_SUB(long tb, long tg);//load subroutine for upper surface of middle/lower block
	long LOAD_SUB2(long tg);//load subroutine for upper surface of upper block
	long SOLVE(long appsCont = 1);//solve for displacement and contact pressure
};

BLOCK_1::BLOCK_1(){
	//
	mkdir("Block", 0777);
	outpDire = "Block/";
	//
	muscSett = (1 << 1);
	//every diviNumb is divisible by all domaNumb
	domaNumb = {3, 3, 3};//only support 1~3
	doleMcsc.assign(3 * domaNumb[0] * domaNumb[1] * domaNumb[2], 1);
	//
	leng = {0.03, 0.025, 0.02};
	diviNumb = {6, 6, 6};
	globLeve = 3;
	loadPres << 0.0E3, 0.0, -1.0E7;
}

long BLOCK_1::ESTA_SURF(){
	bmupSurf.resize(5);
	//
	long totaDivi = diviNumb[0] * (1 << (globLeve));
	bmupSurf[0].indiPoin.resize(totaDivi + 1);
	for(long ti = 0; ti <= totaDivi; ti ++){
		bmupSurf[0].indiPoin[ti].resize(totaDivi + 1);
		for(long tj = 0; tj <= totaDivi; tj ++){
			COOR tempCoor(- leng[0] / 2.0 + leng[0] / totaDivi * (double)ti,
				- leng[0] / 2.0 + leng[0] / totaDivi * (double)tj, leng[0]);
			bmupSurf[0].INSERT(ti, tj, tempCoor);
		}
	}
	//
	totaDivi = diviNumb[1] * (1 << (globLeve));
	bmupSurf[1].indiPoin.resize(totaDivi + 1);
	bmupSurf[2].indiPoin.resize(totaDivi + 1);
	for(long ti = 0; ti <= totaDivi; ti ++){
		bmupSurf[1].indiPoin[ti].resize(totaDivi + 1);
		bmupSurf[2].indiPoin[ti].resize(totaDivi + 1);
		for(long tj = 0; tj <= totaDivi; tj ++){
			//
			COOR tempCoor(- leng[1] / 2.0 + leng[1] / totaDivi * (double)ti,
				- leng[1] / 2.0 + leng[1] / totaDivi * (double)tj, leng[0]);
			bmupSurf[1].INSERT(ti, tj, tempCoor);
			//
			tempCoor[2] += leng[1];
			bmupSurf[2].INSERT(ti, tj, tempCoor);
		}
	}
	//
	totaDivi = diviNumb[2] * (1 << (globLeve));
	bmupSurf[3].indiPoin.resize(totaDivi + 1);
	bmupSurf[4].indiPoin.resize(totaDivi + 1);
	for(long ti = 0; ti <= totaDivi; ti ++){
		bmupSurf[3].indiPoin[ti].resize(totaDivi + 1);
		bmupSurf[4].indiPoin[ti].resize(totaDivi + 1);
		for(long tj = 0; tj <= totaDivi; tj ++){
			//
			COOR tempCoor(- leng[2] / 2.0 + leng[2] / totaDivi * (double)ti,
				- leng[2] / 2.0 + leng[2] / totaDivi * (double)tj, leng[0] + leng[1]);
			bmupSurf[3].INSERT(ti, tj, tempCoor);
			//
			tempCoor[2] += leng[2];
			bmupSurf[4].INSERT(ti, tj, tempCoor);
		}
	}
	std::vector<long> diviNumb_es = {
		diviNumb[0] * (1 << (globLeve)), 
		diviNumb[1] * (1 << (globLeve)), 
		diviNumb[2] * (1 << (globLeve))
	};
	VECT_RESI(blocSurf, 3 * domaNumb[0] * domaNumb[1] * domaNumb[2], 6);
	#pragma omp parallel for
	for(long tg = 0; tg < blocSurf.size(); tg ++){
		//
		long tb = tg / (domaNumb[0] * domaNumb[1] * domaNumb[2]);
		long tg_b = tg % (domaNumb[0] * domaNumb[1] * domaNumb[2]);
		long tg_0 = tg_b / (domaNumb[1] * domaNumb[2]);
		long tg_1 = (tg_b % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
		long tg_2 = (tg_b % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
		std::vector<long> diviReal = {diviNumb_es[tb] / domaNumb[0], 
			diviNumb_es[tb] / domaNumb[1], diviNumb_es[tb] / domaNumb[2]
		};
		VECT_RESI(blocSurf[tg][0].indiPoin, diviReal[1] + 1, diviReal[2] + 1);
		VECT_RESI(blocSurf[tg][1].indiPoin, diviReal[1] + 1, diviReal[2] + 1);
		VECT_RESI(blocSurf[tg][2].indiPoin, diviReal[2] + 1, diviReal[0] + 1);
		VECT_RESI(blocSurf[tg][3].indiPoin, diviReal[2] + 1, diviReal[0] + 1);
		VECT_RESI(blocSurf[tg][4].indiPoin, diviReal[0] + 1, diviReal[1] + 1);
		VECT_RESI(blocSurf[tg][5].indiPoin, diviReal[0] + 1, diviReal[1] + 1);
		//nodes
		for(long ti = 0; ti <= diviReal[0]; ti ++){
			long ti_real = tg_0 * diviReal[0] + ti;
			double xCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb_es[tb] * (double)ti_real;
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				long tj_real = tg_1 * diviReal[1] + tj;
				double yCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb_es[tb] * (double)tj_real;
				for(long tk = 0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					double zCoo = leng[tb] / diviNumb_es[tb] * (double)tk_real;
					if(tb == 1){
						zCoo += leng[0];
					}
					else if(tb == 2){
						zCoo += leng[0] + leng[1];
					}
					COOR tempCoor(xCoo, yCoo, zCoo);
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

long BLOCK_1::MESH(){
	ESTA_SURF();
	multGrid.resize(3 * domaNumb[0] * domaNumb[1] * domaNumb[2]);
	#pragma omp parallel for
	for(long tg = 0; tg < multGrid.size(); tg ++){
		long tb = tg / (domaNumb[0] * domaNumb[1] * domaNumb[2]);
		long tg_b = tg % (domaNumb[0] * domaNumb[1] * domaNumb[2]);
		long tg_0 = tg_b / (domaNumb[1] * domaNumb[2]);
		long tg_1 = (tg_b % (domaNumb[1] * domaNumb[2])) / domaNumb[2];
		long tg_2 = (tg_b % (domaNumb[1] * domaNumb[2])) % domaNumb[2];
		std::vector<long> diviReal = {diviNumb[tb] / domaNumb[0], 
			diviNumb[tb] / domaNumb[1], diviNumb[tb] / domaNumb[2]
		};
		//nodes
		long tempNode[diviReal[0] + 1][diviReal[1] + 1][diviReal[2] + 1];
		for(long ti = 0; ti <= diviReal[0]; ti ++){
			long ti_real = tg_0 * diviReal[0] + ti;
			double xCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (double)ti_real;
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				long tj_real = tg_1 * diviReal[1] + tj;
				double yCoo = - leng[tb] / 2.0 + leng[tb] / diviNumb[tb] * (double)tj_real;
				for(long tk = 0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					double zCoo = leng[tb] / diviNumb[tb] * (double)tk_real;
					if(tb == 1){
						zCoo += leng[0];
					}
					else if(tb == 2){
						zCoo += leng[0] + leng[1];
					}
					COOR tempCoor(xCoo, yCoo, zCoo);
					tempNode[ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		//elements
		for(long ti = 0; ti < diviReal[0]; ti ++){
			for(long tj = 0; tj < diviReal[1]; tj ++){
				for(long tk = 0; tk < diviReal[2]; tk ++){
					TREE_ELEM tempElem;
					tempElem.parent = -1;
					tempElem.cornNode[0] = tempNode[ti][tj][tk];
					tempElem.cornNode[1] = tempNode[ti + 1][tj][tk];
					tempElem.cornNode[2] = tempNode[ti + 1][tj + 1][tk];
					tempElem.cornNode[3] = tempNode[ti][tj + 1][tk];
					tempElem.cornNode[4] = tempNode[ti][tj][tk + 1];
					tempElem.cornNode[5] = tempNode[ti + 1][tj][tk + 1];
					tempElem.cornNode[6] = tempNode[ti + 1][tj + 1][tk + 1];
					tempElem.cornNode[7] = tempNode[ti][tj + 1][tk + 1];
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
		for(long tr = 0; tr < globLeve; tr ++){
			//global refinement level tr
			spliElem.clear();
			for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
				if(multGrid[tg].elemVect[ti].children.size() > 0){
					continue;
				}
				spliElem.insert(ti);
				multGrid[tg].elemVect[ti].refiPatt = 0;
			}
			planSurf.clear();
			multGrid[tg].REFINE(spliElem, spliFlag, planSurf);
		}
		//
		multGrid[tg].OUTPUT_ELEMENT(tg);
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			if((iterNoco.second)[2] <= 1.0E-10){//displacement constraint
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
				// multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
				// multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
			}
			if((iterNoco.second)[0] <= - leng[tb] / 2.0 + 1.0E-12){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
			}
			if((iterNoco.second)[1] <= - leng[tb] / 2.0 + 1.0E-12){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
			}
		}
		//load must at next
		if(tg_2 == domaNumb[2] - 1){
			if(tb == 2){
				LOAD_SUB2(tg);
			}
			else{
				LOAD_SUB(tb, tg);
			}
		}
		//
		multGrid[tg].coupReps = -1;
	}
	return 1;
}

long BLOCK_1::LOAD_SUB(long tb, long tg){
	//
	double zCoo = (tb == 0) ? leng[0] : leng[0] + leng[1];
	std::vector<std::vector<Eigen::Vector3d>> pseuElem;
	VECT_RESI(pseuElem, 4, 4);
	pseuElem[0][0] << - leng[0 + tb] / 2.0, - leng[0 + tb] / 2.0, zCoo;
	pseuElem[0][1] << - leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[0][2] << leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[0][3] << leng[0 + tb] / 2.0, - leng[0 + tb] / 2.0, zCoo;
	pseuElem[1][0] << - leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[1][1] << - leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	pseuElem[1][2] << - leng[1 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	pseuElem[1][3] << - leng[1 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[2][0] << leng[1 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[2][1] << leng[1 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	pseuElem[2][2] << leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	pseuElem[2][3] << leng[0 + tb] / 2.0, - leng[1 + tb] / 2.0, zCoo;
	pseuElem[3][0] << - leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	pseuElem[3][1] << - leng[0 + tb] / 2.0, leng[0 + tb] / 2.0, zCoo;
	pseuElem[3][2] << leng[0 + tb] / 2.0, leng[0 + tb] / 2.0, zCoo;
	pseuElem[3][3] << leng[0 + tb] / 2.0, leng[1 + tb] / 2.0, zCoo;
	//
	EFACE_SURFACE iterEfsu(&(multGrid[tg]), &(bmupSurf[(tb == 0) ? 0 : 2]));
	while(iterEfsu.INCREMENT() == 1){
		std::vector<long> mastNode = iterEfsu.currNode;
		std::vector<Eigen::Vector3d> mastCorn(4);
		for(long ti = 0; ti < 4; ti ++){
			auto iterNoco = multGrid[tg].nodeCoor.find(mastNode[ti]);
			mastCorn[ti] << (iterNoco->second)[0], (iterNoco->second)[1], (iterNoco->second)[2];
		}
		for(long ti = 0; ti < pseuElem.size(); ti ++){
			std::vector<Eigen::Vector3d> slavCorn = pseuElem[ti];
			std::vector<Eigen::Vector2d> listXiet; 
			std::vector<double> listWeig;
			CSEARCH tempCsea;
			tempCsea.SI_SUB(mastCorn, slavCorn, listXiet, listWeig);
			for(long tj = 0; tj < listXiet.size(); tj ++){
				Eigen::Vector<double,4> N_e;
				for(long tk = 0; tk < 4; tk ++){
					N_e(tk) = (1.0 + biliQuad.nacoCorn[tk][0] * listXiet[tj](0)) 
						* (1.0 + biliQuad.nacoCorn[tk][1] * listXiet[tj](1)) / 4.0;
				}
				double weigFact;
				BIQU_TRJA(listXiet[tj], mastCorn, weigFact);
				Eigen::Matrix<double,3,12> N_e_0;
				N_e_0 << 
					N_e(0), 0.0, 0.0, N_e(1), 0.0, 0.0, N_e(2), 0.0, 0.0, N_e(3), 0.0, 0.0,
					0.0, N_e(0), 0.0, 0.0, N_e(1), 0.0, 0.0, N_e(2), 0.0, 0.0, N_e(3), 0.0, 
					0.0, 0.0, N_e(0), 0.0, 0.0, N_e(1), 0.0, 0.0, N_e(2), 0.0, 0.0, N_e(3);
				Eigen::MatrixXd F_e = listWeig[tj] * weigFact * N_e_0.transpose() * loadPres;
				for(long tk = 0; tk < 4; tk ++){
					for(long tm = 0; tm < 3; tm ++){
						long free_tm = 3 * mastNode[tk] + tm;
						multGrid[tg].LOAD_ACCU(free_tm, F_e(3 * tk + tm));
					}
				}
			}
		}
	}
	return 1;
}

long BLOCK_1::LOAD_SUB2(long tg){
	EFACE_SURFACE iterEfsu(&(multGrid[tg]), &(bmupSurf[4]));
	while(iterEfsu.INCREMENT() == 1){
		std::vector<long> tempNode = iterEfsu.currNode;
		std::vector<Eigen::Vector3d> elemCoor(4);
		for(long tk = 0; tk < 4; tk ++){
			auto iterNoco = multGrid[tg].nodeCoor.find(tempNode[tk]);
			elemCoor[tk] << (iterNoco->second)[0], 
				(iterNoco->second)[1], (iterNoco->second)[2];
		}
		Eigen::Matrix<double,12,1> tempForc = Eigen::MatrixXd::Zero(12,1);
		for(long tk = 0; tk < biliQuad.numbNgip; tk ++){
			double weigFact;
			BIQU_TRJA(biliQuad.nacoNgip[tk], elemCoor, weigFact);
			tempForc += biliQuad.niwfNgip[tk] * biliQuad.shfuNgip[tk].transpose() 
				* loadPres * weigFact;
		}
		for(long tk = 0; tk < 4; tk ++){
			for(long tm = 0; tm < 3; tm ++){
				long free_tk = 3 * tempNode[tk] + tm;
				multGrid[tg].LOAD_ACCU(free_tk, tempForc(3 * tk + tm));
			}
		}
	}
	return 1;
}

long BLOCK_1::SOLVE(long appsCont){
	//
	MESH();
	//
	std::vector<std::vector<std::vector<long>>> copaPrec(2);
	for(long tb = 0; tb < 2; tb ++){
		for(long tg_0 = 0; tg_0 < domaNumb[0]; tg_0 ++){
			double x_min = - leng[tb] / 2.0 + leng[tb] / domaNumb[0] * tg_0;
			double x_max = - leng[tb] / 2.0 + leng[tb] / domaNumb[0] * (tg_0 + 1);
			for(long tg_1 = 0; tg_1 < domaNumb[1]; tg_1 ++){
				double y_min = - leng[tb] / 2.0 + leng[tb] / domaNumb[1] * tg_1;
				double y_max = - leng[tb] / 2.0 + leng[tb] / domaNumb[1] * (tg_1 + 1);
				for(long tg_0_ = 0; tg_0_ < domaNumb[0]; tg_0_ ++){
					double x_min_ = 
						- leng[1 + tb] / 2.0 + leng[1 + tb] / domaNumb[0] * tg_0_;
					double x_max_ = 
						- leng[1 + tb] / 2.0 + leng[1 + tb] / domaNumb[0] * (tg_0_ + 1);
					if(x_min_ > x_max - 1.0E-10 || x_max_ < x_min + 1.0E-10){
						continue;
					}
					for(long tg_1_ = 0; tg_1_ < domaNumb[1]; tg_1_ ++){
						double y_min_ = 
							- leng[1 + tb] / 2.0 + leng[1 + tb] / domaNumb[1] * tg_1_;
						double y_max_ = 
							- leng[1 + tb] / 2.0 + leng[1 + tb] / domaNumb[1] * (tg_1_ + 1);
						if(y_min_ > y_max - 1.0E-10 || y_max_ < y_min + 1.0E-10){
							continue;
						}
						std::vector<long> tempIndi = {
							tb * domaNumb[0] * domaNumb[1] * domaNumb[2] 
							+ tg_0 * domaNumb[1] * domaNumb[2] 
							+ tg_1 * domaNumb[2] + domaNumb[2] - 1, 
							(tb + 1) * domaNumb[0] * domaNumb[1] * domaNumb[2] 
							+ tg_0_ * domaNumb[1] * domaNumb[2] 
							+ tg_1_ * domaNumb[2] + 0
						};
						copaPrec[tb].emplace_back(tempIndi);
					}
				}
			}
		}
	}
	//
	long summCopa = (copaPrec[0]).size() + (copaPrec[1]).size();
	std::vector<long> xyzN = {
		3 * (domaNumb[0] - 1) * domaNumb[1] * domaNumb[2], 
		3 * (domaNumb[1] - 1) * domaNumb[0] * domaNumb[2], 
		3 * (domaNumb[2] - 1) * domaNumb[0] * domaNumb[1], 
		summCopa
	};
	//
	searCont.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]);
	penaFact_n.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]);
	penaFact_f.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]);
	fricCoef.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]);
	contBody.resize(xyzN[0] + xyzN[1] + xyzN[2] + xyzN[3]);
	for(long tb = 0; tb < 3; tb ++){
		for(long tg_0 = 0; tg_0 < domaNumb[0]; tg_0 ++){
			for(long tg_1 = 0; tg_1 < domaNumb[1]; tg_1 ++){
				for(long tg_2 = 0; tg_2 < domaNumb[2]; tg_2 ++){
					long tg_m = tb * domaNumb[0] * domaNumb[1] * domaNumb[2]
						+ tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2;
					if(tg_0 <= domaNumb[0] - 2){
						long tg_s = tg_m + domaNumb[1] * domaNumb[2];
						long ts = tg_0 * domaNumb[1] * domaNumb[2] + tg_1 * domaNumb[2] + tg_2 
							+ tb * (domaNumb[0] - 1) * domaNumb[1] * domaNumb[2];
						contBody[ts] = {tg_m, tg_s};
					}
					if(tg_1 <= domaNumb[1] - 2){
						long tg_s = tg_m + domaNumb[2];
						long ts = xyzN[0] 
							+ tg_1 * domaNumb[0] * domaNumb[2] + tg_0 * domaNumb[2] + tg_2 
							+ tb * (domaNumb[1] - 1) * domaNumb[0] * domaNumb[2];
						contBody[ts] = {tg_m, tg_s};
					}
					if(tg_2 <= domaNumb[2] - 2){
						long tg_s = tg_m + 1;
						long ts = xyzN[0] + xyzN[1] 
							+ tg_2 * domaNumb[0] * domaNumb[1] + tg_0 * domaNumb[1] + tg_1 
							+ tb * (domaNumb[2] - 1) * domaNumb[0] * domaNumb[1];
						contBody[ts] = {tg_m, tg_s};
					}
				}
			}
		}
	}
	for(long ti = 0; ti < copaPrec[0].size(); ti ++){
		contBody[xyzN[0] + xyzN[1] + xyzN[2] + ti] = {copaPrec[0][ti][0], copaPrec[0][ti][1]};
	}
	for(long ti = 0; ti < copaPrec[1].size(); ti ++){
		contBody[xyzN[0] + xyzN[1] + xyzN[2] + copaPrec[0].size() + ti] = 
			{copaPrec[1][ti][0], copaPrec[1][ti][1]};
	}
	for(long ts = 0; ts < searCont.size(); ts ++){
		searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
		searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
		penaFact_n[ts] = 210.0E9 * 1000.0;
		penaFact_f[ts] = 210.0E9 * 1000.0;
		if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
			fricCoef[ts] = -1.0;
		}
		else{
			fricCoef[ts] = 0.0;
		}
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
		else if(ts < xyzN[0] + xyzN[1] + xyzN[2] + copaPrec[0].size()){
			mastSurf = &(bmupSurf[0]);
			slavSurf = &(bmupSurf[1]);
		}
		else{
			mastSurf = &(bmupSurf[2]);
			slavSurf = &(bmupSurf[3]);
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
				if(ts < xyzN[0]){
					tempXico += tempCoor[1];
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += tempCoor[0];
					tempEtac += tempCoor[2];
				}
				else/* if(ts < xyzN[0] + xyzN[1] + xyzN[2])*/{
					tempXico += tempCoor[0];
					tempEtac += tempCoor[1];
				}
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		if(ts < xyzN[0]){
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[1] / domaNumb[1] * (1 << (globLeve)), 
				diviNumb[2] / domaNumb[2] * (1 << (globLeve))
			});
		}
		else if(ts < xyzN[0] + xyzN[1]){
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[0] / domaNumb[0] * (1 << (globLeve)), 
				diviNumb[2] / domaNumb[2] * (1 << (globLeve))
			});
		}
		else/* if(ts < xyzN[0] + xyzN[1] + xyzN[2])*/{
			searCont[ts].BUCKET_SORT(mastCoor, {
				diviNumb[0] / domaNumb[0] * (1 << (globLeve)), 
				diviNumb[1] / domaNumb[1] * (1 << (globLeve))
			});
		}
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				COOR tempCoor = multGrid[contBody[ts][1]].nodeCoor[node_tj];
				if(ts < xyzN[0]){
					tempXico += tempCoor[1];
					tempEtac += tempCoor[2];
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += tempCoor[0];
					tempEtac += tempCoor[2];
				}
				else/* if(ts < xyzN[0] + xyzN[1] + xyzN[2])*/{
					tempXico += tempCoor[0];
					tempEtac += tempCoor[1];
				}
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
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
