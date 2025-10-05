/*
 * The GPL-3.0 license
 *
 * Copyright (c) [2025] [Quancheng Peng/Shenyang Ligong University]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * \file   BEAM.h
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  flexible beam example.
 */

#include "../MCONTACT.h"

constexpr double struScalFact = 1.0;//structure scale factor

class BEAM: public MCONTACT{
public:
	//
	std::vector<double> leng;//geometric length
	double lengFact;//length reduction factor from fixed end to free end
	double angl;//geometric angle
	std::vector<long> diviNumb;//number of elements along different directions
	long globLeve;//multigrid level, total levels of global mesh refinement
	std::vector<long> domaNumb;//number of decomposed domains
	//
	double loadInte;//centerline load intensity
	std::vector<std::vector<CURVEDS>> blocSurf;//interfaces between domains
	long coloSett;//constraint load settings
	BEAM(long tempColo);
	long SUBR_ROTA(Eigen::Vector3d &tempCoor, long dire);//rotate about X-axis
	long COOR_ADJU(long tg);//coordinate adjustment
	long SUBR_COLO(long tg, long loadType);//subroutine for displacement constraint and load
	long ESTA_SURF();//establish surface
	long MESH_DD(long loadType);//domain decomposition
	long MESH_NODD(long loadType);//no domain decomposition
	long SOLVE(long appsCont = 1, long noddPrec = 1, long loadType = 0);//solve for displacement
	long SOLVE_DD(long appsCont);//domain decomposition
	long SOLVE_NODD(long noddPrec);//no domain decomposition
	double charFact = 25.0;
};

BEAM::BEAM(long tempColo){
	//
	mkdir("Beam", 0777);
	outpDire = "Beam/";
	//
	muscSett = (1 << 1);
	domaNumb = {32, 4, 2};//each diviNumb is divisible by each domaNumb
	doleMcsc.assign(domaNumb[0] * domaNumb[1] * domaNumb[2], 1);
	//
	leng = {1.0 * struScalFact, 0.12 * struScalFact, 0.06 * struScalFact};
	lengFact = 1.0 / 3.0;
	angl = 45.0 * PI / 180.0;
	diviNumb = {64, 4, 2};//must be even
	globLeve = 4;
	loadInte = - 8000.0;
	coloSett = tempColo;
}

long BEAM::SUBR_ROTA(Eigen::Vector3d &tempCoor, long dire){
	double angl_temp = dire * (angl * tempCoor(0) / leng[0]);
	Eigen::Matrix3d rotaMatr;
	rotaMatr << 1.0, 0.0, 0.0, 
		0.0, cos(angl_temp), - sin(angl_temp), 
		0.0, sin(angl_temp), cos(angl_temp);
	tempCoor = rotaMatr * tempCoor;
	return 1;
}

long BEAM::COOR_ADJU(long tg){
	multGrid[tg].coorNode.clear();
	for(auto& iterNoco : multGrid[tg].nodeCoor){
		Eigen::Vector3d tempCoor;
		tempCoor << (iterNoco.second)[0], (iterNoco.second)[1], (iterNoco.second)[2];
		SUBR_ROTA(tempCoor, 1);
		iterNoco.second = {tempCoor(0), tempCoor(1), tempCoor(2)};
		multGrid[tg].coorNode.emplace(iterNoco.second, iterNoco.first);
	}
	return 1;
}

long BEAM::SUBR_COLO(long tg, long loadType){
	//displacement constraint must at first
	for(const auto& iterNoco : multGrid[tg].nodeCoor){
		if((iterNoco.second)[0] <= 1.0E-10){//displacement constraint
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
			multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
		}
	}
	//load must at next
	for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
		if(multGrid[tg].elemVect[ti].children.size() > 0){
			continue;
		}
		if(loadType == 0){
			for(long tj = 0; tj < hexaLine.size(); tj ++){
				std::vector<long> inpuNode(hexaLine[tj].size());
				std::vector<COOR> inpuCoor(hexaLine[tj].size());
				bool tempFlag = true;
				for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
					inpuNode[tk] = multGrid[tg].elemVect[ti].cornNode[hexaLine[tj][tk]];
					auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
					inpuCoor[tk] = iterNoco->second;
					if(abs((iterNoco->second)[1]) > 1.0E-10 
						|| abs((iterNoco->second)[2]) > 1.0E-10){
						tempFlag = false;
						break;
					}
				}
				if(tempFlag == false){
					continue;
				}
				// / 4.0: must (four elements share one edge)
				double forcValu = loadInte * abs(inpuCoor[0][0] - inpuCoor[1][0]) / 2.0 / 4.0;
				for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
					if(inpuCoor[tk][0] > 1.0E-10){
						multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 2, forcValu);
					}
				}
			}
		}
		else{
			double faceInte = leng[0] * loadInte / (leng[1] * leng[2] * pow(1.0 - lengFact, 2.0));
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
					if((iterNoco->second)[0] <= leng[0] - 1.0E-12){
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
					Eigen::Vector3d tempPres;
					tempPres << 0.0, 0.0, faceInte;
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
	}
	//
	multGrid[tg].coupReps = -1;
	return 1;
}

long BEAM::ESTA_SURF(){
	std::vector<long> diviNumb_es = {
		diviNumb[0] * (1 << (globLeve)), 
		diviNumb[1] * (1 << (globLeve)), 
		diviNumb[2] * (1 << (globLeve))
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
			double xCoo = leng[0] / diviNumb_es[0] * ti_real;
			double heig_ti = leng[1] * (1.0- (double)ti_real / diviNumb_es[0] * lengFact);
			double widt_ti = leng[2] * (1.0- (double)ti_real / diviNumb_es[0] * lengFact);
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				long tj_real = tg_1 * diviReal[1] + tj;
				double yCoo = - heig_ti / 2.0 + heig_ti / diviNumb_es[1] * (double)tj_real;
				for(long tk =  0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					double zCoo = - widt_ti / 2.0 + widt_ti / diviNumb_es[2] * (double)tk_real;
					Eigen::Vector3d origCoor;
					origCoor << xCoo, yCoo, zCoo;
					SUBR_ROTA(origCoor, 1);
					COOR tempCoor = {origCoor(0), origCoor(1), origCoor(2)};
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

long BEAM::MESH_NODD(long loadType){
	multGrid.resize(1);
	long tg = 0;
	//nodes
	long tempNode[diviNumb[0] + 1][diviNumb[1] + 1][diviNumb[2] + 1];
	for(long ti = 0; ti <= diviNumb[0]; ti ++){
		double heig_ti = leng[1] * (1.0- (double)ti / diviNumb[0] * lengFact);
		double widt_ti = leng[2] * (1.0- (double)ti / diviNumb[0] * lengFact);
		for(long tj = 0; tj <= diviNumb[1]; tj ++){
			double yCoo = - heig_ti / 2.0 + heig_ti / diviNumb[1] * tj;
			for(long tk = 0; tk <= diviNumb[2]; tk ++){
				double zCoo = - widt_ti / 2.0 + widt_ti / diviNumb[2] * tk;
				COOR tempCoor(leng[0] / diviNumb[0] * ti, yCoo, zCoo);
				tempNode[ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
			}
		}
	}
	//elements
	for(long ti = 0; ti < diviNumb[0]; ti ++){
		for(long tj = 0; tj < diviNumb[1]; tj ++){
			for(long tk = 0; tk < diviNumb[2]; tk ++){
				TREE_ELEM tempElem;
				tempElem.parent = -1;
				tempElem.cornNode[0] = tempNode[ti][tj][tk];
				tempElem.cornNode[1] = tempNode[ti][tj + 1][tk];
				tempElem.cornNode[2] = tempNode[ti][tj + 1][tk + 1];
				tempElem.cornNode[3] = tempNode[ti][tj][tk + 1];
				tempElem.cornNode[4] = tempNode[ti + 1][tj][tk];
				tempElem.cornNode[5] = tempNode[ti + 1][tj + 1][tk];
				tempElem.cornNode[6] = tempNode[ti + 1][tj + 1][tk + 1];
				tempElem.cornNode[7] = tempNode[ti + 1][tj][tk + 1];
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
	COOR_ADJU(tg);
	multGrid[tg].OUTPUT_ELEMENT(tg);
	//
	SUBR_COLO(tg, loadType);
	return 1;
}

long BEAM::MESH_DD(long loadType){
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
			double xCoo = leng[0] / diviNumb[0] * ti_real;
			double heig_ti = leng[1] * (1.0- (double)ti_real / diviNumb[0] * lengFact);
			double widt_ti = leng[2] * (1.0- (double)ti_real / diviNumb[0] * lengFact);
			for(long tj = 0; tj <= diviReal[1]; tj ++){
				blocNode[ti][tj].resize(diviReal[2] + 1);
				long tj_real = tg_1 * diviReal[1] + tj;
				double yCoo = - heig_ti / 2.0 + heig_ti / diviNumb[1] * (double)tj_real;
				for(long tk =  0; tk <= diviReal[2]; tk ++){
					long tk_real = tg_2 * diviReal[2] + tk;
					double zCoo = - widt_ti / 2.0 + widt_ti / diviNumb[2] * (double)tk_real;
					COOR tempCoor(xCoo, yCoo, zCoo);
					blocNode[ti][tj][tk] = multGrid[tg].TRY_ADD_NODE(tempCoor);
				}
			}
		}
		//elements
		for(long ti = 0; ti < diviReal[0]; ti ++){
			for(long tj = 0; tj < diviReal[1]; tj ++){
				for(long tk = 0; tk < diviReal[2]; tk ++){
					TREE_ELEM tempElem;
					tempElem.parent = -1;
					tempElem.cornNode[0] = blocNode[ti][tj][tk];
					tempElem.cornNode[1] = blocNode[ti][tj + 1][tk];
					tempElem.cornNode[2] = blocNode[ti][tj + 1][tk + 1];
					tempElem.cornNode[3] = blocNode[ti][tj][tk + 1];
					tempElem.cornNode[4] = blocNode[ti + 1][tj][tk];
					tempElem.cornNode[5] = blocNode[ti + 1][tj + 1][tk];
					tempElem.cornNode[6] = blocNode[ti + 1][tj + 1][tk + 1];
					tempElem.cornNode[7] = blocNode[ti + 1][tj][tk + 1];
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
		COOR_ADJU(tg);
		multGrid[tg].OUTPUT_ELEMENT(tg);
		//
		SUBR_COLO(tg, loadType);
	}
	return 1;
}

long BEAM::SOLVE(long appsCont, long noddPrec, long loadType){
	if(coloSett == 0){
		MESH_NODD(loadType);
		SOLVE_NODD(noddPrec);
	}
	else{
		ESTA_SURF();
		MESH_DD(loadType);
		SOLVE_DD(appsCont);
	}
	return 1;
}

long BEAM::SOLVE_NODD(long noddPrec){
	long tg = 0;
	//
	multGrid[tg].TRANSFER();
	multGrid[tg].STIF_MATR();
	Eigen::VectorXd resuSolu;
	multGrid[tg].CONSTRAINT(1);
	multGrid[tg].mgpi.CG_SOLV(noddPrec, multGrid[tg].consForc, resuSolu);
	//not advisable
	// DIRE_SOLV direSolv;
	// direSolv.compute(multGrid[tg].mgpi.consStif[multGrid[tg].mgpi.maxiLeve]);
	// resuSolu = direSolv.solve(multGrid[tg].consForc);
	// multGrid[tg].mgpi.MULT_SOLV(multGrid[tg].consForc, resuSolu);
	// multGrid[tg].mgpi.GMRES_SOLV(multGrid[tg].consForc, resuSolu);
	Eigen::VectorXd outpDisp;
	multGrid[tg].OUTP_SUB1(resuSolu, outpDisp);
	multGrid[tg].OUTP_SUB2(outpDisp, tg);
	multGrid[tg].STRESS_RECOVERY(outpDisp, tg);
	return 1;
}

long BEAM::SOLVE_DD(long appsCont){
	//
	double charLeng = GET_CHAR_LENG();
	//
	std::vector<long> xyzN = {
		(domaNumb[0] - 1) * domaNumb[1] * domaNumb[2], 
		(domaNumb[1] - 1) * domaNumb[0] * domaNumb[2], 
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
				if(tg_1 <= domaNumb[1] - 2){
					long tg_s = tg_m + domaNumb[2];
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
		penaFact_n[ts] = 210.0E9 * charFact / charLeng;
		penaFact_f[ts] = 210.0E9 * charFact / charLeng;
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
				Eigen::Vector3d origCoor;
				origCoor << tempCoor[0], tempCoor[1], tempCoor[2];
				SUBR_ROTA(origCoor, -1);
				if(ts < xyzN[0]){
					tempXico += origCoor(1);
					tempEtac += origCoor(2);
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += origCoor(0);
					tempEtac += origCoor(2);
				}
				else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
					tempXico += origCoor(0);
					tempEtac += origCoor(1);
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
		else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
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
				Eigen::Vector3d origCoor;
				origCoor << tempCoor[0], tempCoor[1], tempCoor[2];
				SUBR_ROTA(origCoor, -1);
				if(ts < xyzN[0]){
					tempXico += origCoor(1);
					tempEtac += origCoor(2);
				}
				else if(ts < xyzN[0] + xyzN[1]){
					tempXico += origCoor(0);
					tempEtac += origCoor(2);
				}
				else if(ts < xyzN[0] + xyzN[1] + xyzN[2]){
					tempXico += origCoor(0);
					tempEtac += origCoor(1);
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
		// initialization
		// for(long tv = 0; tv < multGrid.size(); tv ++){//for MONITOR
			// double tempDisp_y = 0.0 + (4.0E-4 / domaNumb[0]) * (tv / domaNumb[1]);
			// double tempDisp_z = 0.0 + (-0.002888 / domaNumb[0]) * (tv / domaNumb[1]);
			// for(const auto &iterMgnc : multGrid[tv].nodeCoor){
				// (resuDisp[tv])(3 * iterMgnc.first + 1) = tempDisp_y;
				// (resuDisp[tv])(3 * iterMgnc.first + 2) = tempDisp_z;
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
				// double tempDisp_y = 0.0 + (4.0E-4 / domaNumb[0]) * (tempBody / domaNumb[1]);
				// double tempDisp_z = 0.0 + (-0.002888 / domaNumb[0]) * (tempBody / domaNumb[1]);
				// for(const auto iterNoco : nodeCont[ts][tv]){
					// (inteAuxi[ts][tv])(3 * iterNoco.second + 1) = tempDisp_y;
					// (inteAuxi[ts][tv])(3 * iterNoco.second + 2) = tempDisp_z;
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
