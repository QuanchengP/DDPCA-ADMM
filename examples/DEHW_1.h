
#include "DEHWSURF.h"
#include "../MCONTACT.h"

class DEHW_1: public MCONTACT{
public:
	//****************************************************************************************
	long coloSett;//constraint load settings: is/not self locking
	DEHW_1(long tempColo);
	//****************************************************************************************
	DEHWSURF dehwSurf;//discretization of worm/wheel tooth surface
	double centErro;//center distance error
	Eigen::Vector2d analAngl;//rotating angle of worm and worm wheel
	std::vector<double> distCrit;//critical gap used for adaptive mesh refinement
	std::vector<double> loadIncr;//load intensity of the worm inner hub
	//****************************************************************************************
	//when mesh refinement, average inpuCoor to get outpCoor for worm
	long COOR_AVER(const std::vector<COOR> &inpuCoor, COOR &outpCoor);
	//when mesh refinement, average inpuCoor to get outpCoor for wheel
	long COOR_AVER_1(const std::vector<COOR> &inpuCoor, COOR &outpCoor);
	long WORM_MESH(long isnoDode);//generate grid of worm
	long WHEE_MESH_DD();//generate grid of worm wheel: domain decomposition
	long WHEE_MESH_NODD();//generate grid of worm wheel: no demain  decomposition
	//update wodeAucu
	long UPDA_WODE(const std::vector<long> &inpuNode, 
		const COOR &outpCoor, long face_tw, long face_sub, long worm_tv
	);
	//update whdeAucu
	long UPDA_WHDE(const std::vector<long> &inpuNode, 
		const COOR &outpCoor, long face_tw, long whee_tv
	);
	long SUBR_COLO_WORM(long tg);//subroutin of displacement constraint and load for worm
	long SUBR_COLO_WHEE(long tg);//subroutin of displacement constraint and load for wheel
	//****************************************************************************************
	std::vector<double> wodeAuan;//worm decomposition auxiliary angle
	std::vector<std::vector<CURVEDS>> wodeAucu;//worm decomposition auxiliary curved surface
	std::vector<CURVEDS> whdeAucu;//wheel decomposition auxiliary curved surface
	long CONT_INTE_DD();//establish contact interfaces: domain decomposition
	long CONT_INTE_NODD();//establish contact interfaces: no demain  decomposition
	//solve for displacement and contact pressure
	//appsCont(automatic penalty parameter selection): 0 - eigen analysis, 1 - contact analysis
	//isnoDode: 0 - no DD, 1 - DD(domain decomposition)
	long SOLVE(long appsCont = 1, long isnoDode = 1);
};

DEHW_1::DEHW_1(long tempColo){
	//
	mkdir("Dehw", 0777);
	outpDire = "Dehw/";
	coloSett = tempColo;
}

long DEHW_1::COOR_AVER(const std::vector<COOR> &inpuCoor, COOR &outpCoor){
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

long DEHW_1::COOR_AVER_1(const std::vector<COOR> &inpuCoor, COOR &outpCoor){
	COOR cyliCoor(0.0, 0.0, 0.0);
	long tempFlag_0 = 0;
	long tempFlag_1 = 0;
	bool inneHubf = true;
	for(long tk = 0; tk < inpuCoor.size(); tk ++){
		double tempRadi = sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + 
			inpuCoor[tk][1] * inpuCoor[tk][1]
		);
		if(abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
			inneHubf = false;
		}
		double tempAngl = atan2(inpuCoor[tk][1], inpuCoor[tk][0]);
		cyliCoor[1] += tempAngl;
		if(tempAngl > PI / 2.0){
			tempFlag_0 ++;
		}
		if(tempAngl < - PI / 2.0){
			tempFlag_1 ++;
		}
	}
	if(inneHubf == true){
		COOR_AVER(inpuCoor, outpCoor);
		return 1;
	}
	if(tempFlag_0 > 0 && tempFlag_1 > 0){
		cyliCoor[1] += tempFlag_1 * (PI * 2.0);
	}
	//
	double toruRadi = 0.0;
	double toruAngl = 0.0;
	for(long tk = 0; tk < inpuCoor.size(); tk ++){
		double tempRadi = dehwSurf.a_h2 + centErro - 
			sqrt(inpuCoor[tk][0] * inpuCoor[tk][0] + inpuCoor[tk][1] * inpuCoor[tk][1]);
		toruRadi += sqrt(pow(tempRadi, 2.0) + pow(inpuCoor[tk][2], 2.0));
		toruAngl += atan2(inpuCoor[tk][2], tempRadi);
	}
	toruRadi /= inpuCoor.size();
	toruAngl /= inpuCoor.size();
	//
	cyliCoor[0] = dehwSurf.a_h2 + centErro - toruRadi * cos(toruAngl);
	cyliCoor[1] /= inpuCoor.size();
	cyliCoor[2] = toruRadi * sin(toruAngl);
	outpCoor.resize(3);
	outpCoor[0] = cyliCoor[0] * cos(cyliCoor[1]);
	outpCoor[1] = cyliCoor[0] * sin(cyliCoor[1]);
	outpCoor[2] = cyliCoor[2];
	return 1;
}

long DEHW_1::SUBR_COLO_WORM(long tg){
	if(tg == -1){
		long tg_endi = (multGrid.size() == 2) ? 1 : dehwSurf.gridNumb[0][6];
		//caculate the total area
		double totaArea = 0.0;
		for(long tg = 0; tg < tg_endi; tg ++){
			for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
				if(multGrid[tg].elemVect[ti].children.size() > 0){
					continue;
				}
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					bool tempFlag = true;
					std::vector<Eigen::Vector3d> elemCoor(4);
					for(long tk = 0; tk < 4; tk ++){
						long tempNode = multGrid[tg].elemVect[ti].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tg].nodeCoor.find(tempNode);
						elemCoor[tk] << (iterNoco->second)[0], 
							(iterNoco->second)[1], (iterNoco->second)[2];
						double x_loca = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
						double y_loca = - (iterNoco->second)[2];
						double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
						if(abs(tempRadi - dehwSurf.inneRadi[0]) > 1.0E-10){
							tempFlag = false;
							break;
						}
					}
					if(tempFlag == false){
						continue;
					}
					double tempArea = 0.0;
					for(long tk = 0; tk < biliQuad.numbNgip; tk ++){
						double weigFact;
						BIQU_TRJA(biliQuad.nacoNgip[tk], elemCoor, weigFact);
						tempArea += biliQuad.niwfNgip[tk] * weigFact;
					}
					totaArea += tempArea;
				}
			}
		}
		std::cout << "The total area of worm inner hub: " << totaArea;
		OUTPUT_TIME("");
		loadIncr[0] = dehwSurf.inpuTorq / dehwSurf.inneRadi[0] / totaArea;
	}
	else if(coloSett == 1){
		//rotate the nodal coordinate system
		//NO: node coupling!
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			double x_loca = (iterNoco.second)[0] + (dehwSurf.a_h2 + centErro);
			double y_loca = - (iterNoco.second)[2];
			double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
			if(abs(tempRadi - dehwSurf.inneRadi[0]) <= 1.0E-10){
				double tempAngl = atan2(y_loca, x_loca);
				Eigen::Matrix3d tempRota;
				tempRota << cos(tempAngl), -sin(tempAngl), 0.0, 
					0.0, 0.0, 1.0, 
					-sin(tempAngl), -cos(tempAngl), 0.0;
				multGrid[tg].nodeRota.emplace(iterNoco.first, tempRota);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
				//
				// multGrid[tg].coupNode.emplace(iterNoco.first);
			}
		}
		// multGrid[tg].coupReps = -1;
		// double miniCoor = 1.0E20;
		// for(const auto iterCono : multGrid[tg].coupNode){
			// auto iterNoco = multGrid[tg].nodeCoor.find(iterCono);
			// if((iterNoco->second)[1] < miniCoor){
				// miniCoor = (iterNoco->second)[1];
				// multGrid[tg].coupReps = iterNoco->first;
			// }
		// }
		// multGrid[tg].coupNode.erase(multGrid[tg].coupReps);
		//load must at next
		double totaForc = 0.0;
		for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
			if(multGrid[tg].elemVect[ti].children.size() > 0){
				continue;
			}
			for(long tj = 0; tj < hexaFace.size(); tj ++){
				bool tempFlag = true;
				std::vector<long> inpuNode(4);
				std::vector<Eigen::Vector3d> elemCoor(4);
				for(long tk = 0; tk < 4; tk ++){
					inpuNode[tk] = multGrid[tg].elemVect[ti].cornNode[hexaFace[tj][tk]];
					auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
					elemCoor[tk] << (iterNoco->second)[0], 
						(iterNoco->second)[1], (iterNoco->second)[2];
					double x_loca = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
					double y_loca = - (iterNoco->second)[2];
					double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
					if(abs(tempRadi - dehwSurf.inneRadi[0]) > 1.0E-10){
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
					Eigen::Vector3d tempLoad;
					tempLoad << 0.0, loadIncr[0], 0.0;
					tempForc += biliQuad.niwfNgip[tk] 
						* biliQuad.shfuNgip[tk].transpose() * tempLoad * weigFact;
				}
				for(long tk = 0; tk < 4; tk ++){
					multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 1, tempForc(3 * tk + 1));
					totaForc += tempForc(3 * tk + 1);
				}
			}
		}
		std::cout << "The total force of worm inner hub: " << totaForc;
		OUTPUT_TIME("");
	}
	else{
		//rotate the nodal coordinate system
		//NO: node coupling!
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			double x_loca = (iterNoco.second)[0] + (dehwSurf.a_h2 + centErro);
			double y_loca = - (iterNoco.second)[2];
			double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
			if(abs(tempRadi - dehwSurf.inneRadi[0]) <= 1.0E-10){
				double tempAngl = atan2(y_loca, x_loca);
				Eigen::Matrix3d tempRota;
				tempRota << cos(tempAngl), -sin(tempAngl), 0.0, 
					0.0, 0.0, 1.0, 
					-sin(tempAngl), -cos(tempAngl), 0.0;
				multGrid[tg].nodeRota.emplace(iterNoco.first, tempRota);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
				multGrid[tg].LOAD_ACCU(3 * iterNoco.first + 1, 1.0E-10);
			}
		}
	}
	return 1;
}

long DEHW_1::SUBR_COLO_WHEE(long tg){
	if(tg == -1){
		long tg_star = (multGrid.size() == 2) ? 1 : dehwSurf.gridNumb[0][6];
		//caculate the total area
		double totaArea = 0.0;
		for(long tg = tg_star; tg < multGrid.size(); tg ++){
			for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
				if(multGrid[tg].elemVect[ti].children.size() > 0){
					continue;
				}
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					bool tempFlag = true;
					std::vector<Eigen::Vector3d> elemCoor(4);
					for(long tk = 0; tk < 4; tk ++){
						long tempNode = multGrid[tg].elemVect[ti].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tg].nodeCoor.find(tempNode);
						elemCoor[tk] << (iterNoco->second)[0], 
							(iterNoco->second)[1], (iterNoco->second)[2];
						double x_loca = (iterNoco->second)[0];
						double y_loca = (iterNoco->second)[1];
						double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
						if(abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
							tempFlag = false;
							break;
						}
					}
					if(tempFlag == false){
						continue;
					}
					double tempArea = 0.0;
					for(long tk = 0; tk < biliQuad.numbNgip; tk ++){
						double weigFact;
						BIQU_TRJA(biliQuad.nacoNgip[tk], elemCoor, weigFact);
						tempArea += biliQuad.niwfNgip[tk] * weigFact;
					}
					totaArea += tempArea;
				}
			}
		}
		std::cout << "The total area of worm wheel inner hub: " << totaArea;
		OUTPUT_TIME("");
		loadIncr[1] = - dehwSurf.inpuTorq * dehwSurf.i_h2 / dehwSurf.inneRadi[1] / totaArea;
	}
	else if(coloSett == 1){
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			double tempRadi = sqrt(pow((iterNoco.second)[0], 2.0) 
				+ pow((iterNoco.second)[1], 2.0)
			);
			if(abs(tempRadi - dehwSurf.inneRadi[1]) <= 1.0E-10){
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 1, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
			}
		}
	}
	else{
		//rotate the nodal coordinate system
		//NO: node coupling!
		//displacement constraint must at first
		for(const auto& iterNoco : multGrid[tg].nodeCoor){
			double x_loca = (iterNoco.second)[0];
			double y_loca = (iterNoco.second)[1];
			double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
			if(abs(tempRadi - dehwSurf.inneRadi[1]) <= 1.0E-10){
				double tempAngl = atan2(y_loca, x_loca);
				Eigen::Matrix3d tempRota;
				tempRota << cos(tempAngl), -sin(tempAngl), 0.0, 
					sin(tempAngl), cos(tempAngl), 0.0, 
					0.0, 0.0, 1.0;
				multGrid[tg].nodeRota.emplace(iterNoco.first, tempRota);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 0, 0.0);
				multGrid[tg].consDofv.emplace(3 * iterNoco.first + 2, 0.0);
			}
		}
		//load must at next
		double totaForc = 0.0;
		for(long ti = 0; ti < multGrid[tg].elemVect.size(); ti ++){
			if(multGrid[tg].elemVect[ti].children.size() > 0){
				continue;
			}
			for(long tj = 0; tj < hexaFace.size(); tj ++){
				bool tempFlag = true;
				std::vector<long> inpuNode(4);
				std::vector<Eigen::Vector3d> elemCoor(4);
				for(long tk = 0; tk < 4; tk ++){
					inpuNode[tk] = multGrid[tg].elemVect[ti].cornNode[hexaFace[tj][tk]];
					auto iterNoco = multGrid[tg].nodeCoor.find(inpuNode[tk]);
					elemCoor[tk] << (iterNoco->second)[0], 
						(iterNoco->second)[1], (iterNoco->second)[2];
					double x_loca = (iterNoco->second)[0];
					double y_loca = (iterNoco->second)[1];
					double tempRadi = sqrt(pow(x_loca, 2.0) + pow(y_loca, 2.0));
					if(abs(tempRadi - dehwSurf.inneRadi[1]) > 1.0E-10){
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
					Eigen::Vector3d tempLoad;
					tempLoad << 0.0, loadIncr[1], 0.0;
					tempForc += biliQuad.niwfNgip[tk] 
						* biliQuad.shfuNgip[tk].transpose() * tempLoad * weigFact;
				}
				for(long tk = 0; tk < 4; tk ++){
					multGrid[tg].LOAD_ACCU(3 * inpuNode[tk] + 1, tempForc(3 * tk + 1));
					totaForc += tempForc(3 * tk + 1);
				}
			}
		}
		std::cout << "The total force of wheel inner hub: " << totaForc;
		OUTPUT_TIME("");
	}
	return 1;
}

long DEHW_1::WORM_MESH(long isnoDode){
	OUTPUT_TIME("DEHW_1::WORM_MESH");
	//
	Eigen::Matrix3d wormRota_0;
	wormRota_0 << cos(analAngl[0]), -sin(analAngl[0]), 0.0, 
		sin(analAngl[0]), cos(analAngl[0]), 0.0, 
		0.0, 0.0, 1.0;
	Eigen::Matrix3d wormRota_1;
	wormRota_1 << 1.0, 0.0, 0.0, 
		0.0, 0.0, 1.0, 
		0.0, -1.0, 0.0;
	wormRota_1 = wormRota_1 * wormRota_0;
	Eigen::Vector3d wormTran;
	wormTran << - (dehwSurf.a_h2 + centErro), 0.0, 0.0;
	if(isnoDode == 1){
		wodeAuan.resize(dehwSurf.gridNumb[0][6] - 1);
		VECT_RESI(wodeAucu, dehwSurf.gridNumb[0][6], 2);
	}
	long wodeFact_0 = (1 << (dehwSurf.globHomo));
	long wodeFact_1 = (1 << (dehwSurf.globInho + dehwSurf.globHomo));
	//points
	long tw_endi = (isnoDode == 0) ? 1 : dehwSurf.gridNumb[0][6];
	#pragma omp parallel for
	for(long tw = 0; tw < tw_endi; tw ++){
		long numb_tw, numbStar;
		if(isnoDode == 0){
			numbStar = 0;
			numb_tw = dehwSurf.gridNumb[0][4] * (dehwSurf.gridNumb[0][6] - 2) 
				+ dehwSurf.gridNumb[0][5] * 2;
		}
		else{
			if(tw == 0 || tw == dehwSurf.gridNumb[0][6] - 1){
				numb_tw = dehwSurf.gridNumb[0][5];
				if(tw == 0){
					numbStar = 0;
				}
				else{
					numbStar = dehwSurf.gridNumb[0][5] + (tw - 1) * dehwSurf.gridNumb[0][4];
				}
			}
			else{
				numb_tw = dehwSurf.gridNumb[0][4];
				numbStar = dehwSurf.gridNumb[0][5] + (tw - 1) * dehwSurf.gridNumb[0][4];
			}
			for(long ti = 0; ti < 2; ti ++){
				VECT_RESI(wodeAucu[tw][ti].indiPoin, 
					(dehwSurf.gridNumb[0][1] + dehwSurf.gridNumb[0][2]) * wodeFact_0 + 1, 
					numb_tw * wodeFact_1 + 1
				);
			}
			if(tw >= 1){
				long ti_real = numbStar 
					* (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
				double tempAuan = dehwSurf.curvCoor[0][ti_real][0](0);
				while(tempAuan > PI){
					tempAuan -= 2.0 * PI;
				}
				wodeAuan[tw - 1] = tempAuan;
			}
		}
		VECTOR4L blocPoin;
		VECT_RESI(blocPoin, 4, numb_tw + 1);
		for(long ti = 0; ti <= numb_tw; ti ++){
			VECT_RESI(blocPoin[0][ti], dehwSurf.gridNumb[0][1] + 1, 
				dehwSurf.gridNumb[0][0] + 1
			);
			VECT_RESI(blocPoin[1][ti], dehwSurf.gridNumb[0][2] + 1, 
				dehwSurf.gridNumb[0][0] / 2 + 1
			);
			VECT_RESI(blocPoin[2][ti], dehwSurf.gridNumb[0][3] + 1, 
				2 * dehwSurf.gridNumb[0][2] + 1
			);
			VECT_RESI(blocPoin[3][ti], dehwSurf.gridNumb[0][2] + 1, 
				dehwSurf.gridNumb[0][0] / 2 + 1
			);
			long ti_real = (numbStar + ti) 
				* (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
			Eigen::Vector<Eigen::Matrix<double,3,2>,Eigen::Dynamic> wormProf;
			wormProf.resize(dehwSurf.gridNumb[0][3] + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][3]; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
				wormProf(tj).block(0,0,3,1) << dehwSurf.wormTosu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wormTosu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wormTosu.indiPoin[ti_real][tj_real][2];
				wormProf(tj).block(0,1,3,1) << dehwSurf.wormToba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wormToba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wormToba.indiPoin[ti_real][tj_real][2];
			}
			ti_real = (numbStar + ti) 
				* (1 << (dehwSurf.globInho + dehwSurf.globHomo));
			Eigen::MatrixXd rootProf_1, rootProf_2;
			rootProf_1.resize(3, dehwSurf.gridNumb[0][0] / 2 + 1);
			rootProf_2.resize(3, dehwSurf.gridNumb[0][0] / 2 + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][0] / 2; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo));
				rootProf_1.block(0,tj,3,1) << dehwSurf.wormRtsu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wormRtsu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wormRtsu.indiPoin[ti_real][tj_real][2];
				rootProf_2.block(0,tj,3,1) << dehwSurf.wormRtba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wormRtba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wormRtba.indiPoin[ti_real][tj_real][2];
			}
			Eigen::Vector2d profRadi;
			profRadi << 
				sqrt(rootProf_1(0,0) * rootProf_1(0,0) + rootProf_1(1,0) * rootProf_1(1,0)), 
				sqrt(rootProf_2(0,0) * rootProf_2(0,0) + rootProf_2(1,0) * rootProf_2(1,0));
			Eigen::Vector2d tranRadi;
			tranRadi(0) = profRadi(0) - PI / 4.0 * dehwSurf.m_t;
			tranRadi(1) = profRadi(1) - PI / 4.0 * dehwSurf.m_t;
			Eigen::Matrix<double,3,4> blocCoor;
			blocCoor << 
				rootProf_1(0,0) * dehwSurf.inneRadi[0] / profRadi(0), 
				rootProf_2(0,0) * dehwSurf.inneRadi[0] / profRadi(1), 
				rootProf_2(0,0) * tranRadi(1) / profRadi(1), 
				rootProf_1(0,0) * tranRadi(0) / profRadi(0), 
				rootProf_1(1,0) * dehwSurf.inneRadi[0] / profRadi(0), 
				rootProf_2(1,0) * dehwSurf.inneRadi[0] / profRadi(1), 
				rootProf_2(1,0) * tranRadi(1) / profRadi(1), 
				rootProf_1(1,0) * tranRadi(0) / profRadi(0), 
				rootProf_1(2,0), rootProf_2(2,0), rootProf_2(2,0), rootProf_1(2,0);
			//0
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][1]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[0][0]; tk ++){
					Eigen::Matrix<double,3,2> downUp;
					downUp.block(0,0,3,1) = 
						(1.0 - (double)tk / dehwSurf.gridNumb[0][0]) 
						* blocCoor.block(0,0,3,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[0][0] * blocCoor.block(0,1,3,1).eval();
					downUp.block(0,1,3,1) = 
						(1.0 - (double)tk / dehwSurf.gridNumb[0][0]) 
						* blocCoor.block(0,3,3,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[0][0] * blocCoor.block(0,2,3,1).eval();
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tj / dehwSurf.gridNumb[0][1]) 
						* downUp.block(0,0,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[0][1] * downUp.block(0,1,3,1).eval();
					COOR inpuCoor = {tempCoor(0), tempCoor(1), tempCoor(2)};
					blocPoin[0][ti][tj][tk] = multGrid[tw].TRY_ADD_NODE(inpuCoor);
					if(isnoDode == 1 && tk == 0){
						wodeAucu[tw][0].INSERT(tj * wodeFact_0, ti * wodeFact_1, inpuCoor);
					}
					else if(isnoDode == 1 && tk == dehwSurf.gridNumb[0][0]){
						wodeAucu[tw][1].INSERT(tj * wodeFact_0, ti * wodeFact_1, inpuCoor);
					}
				}
			}
			//1
			Eigen::MatrixXd lineCoor;
			lineCoor.resize(3,2);
			lineCoor.block(0,0,3,1) = blocCoor.block(0,3,3,1).eval();
			lineCoor.block(0,1,3,1) = 0.5 * blocCoor.block(0,3,3,1).eval() 
				+ 0.5 * blocCoor.block(0,2,3,1).eval();
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[0][0] / 2; tk ++){
					Eigen::Matrix<double,3,2> downUp;
					downUp.block(0,0,3,1) = 
						(1.0 - (double)tk / (dehwSurf.gridNumb[0][0] / 2.0)) 
						* lineCoor.block(0,0,3,1).eval() 
						+ (double)tk / (dehwSurf.gridNumb[0][0]/2.0) 
						* lineCoor.block(0,1,3,1).eval();
					downUp.block(0,1,3,1) = rootProf_1.block(0,tk,3,1);
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tj / dehwSurf.gridNumb[0][2]) 
						* downUp.block(0,0,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[0][2] * downUp.block(0,1,3,1).eval();
					COOR inpuCoor = {tempCoor(0), tempCoor(1), tempCoor(2)};
					blocPoin[1][ti][tj][tk] = multGrid[tw].TRY_ADD_NODE(inpuCoor);
					if(isnoDode == 1 && tk == 0){
						wodeAucu[tw][0].INSERT(
							(dehwSurf.gridNumb[0][1] + tj) * wodeFact_0, 
							ti * wodeFact_1, inpuCoor
						);
					}
				}
			}
			//2
			lineCoor.block(0,0,3,1) = lineCoor.block(0,1,3,1).eval();
			lineCoor.block(0,1,3,1) = 0.5 * wormProf(dehwSurf.gridNumb[0][3]).block(0,0,3,1) 
				+ 0.5 * wormProf(dehwSurf.gridNumb[0][3]).block(0,1,3,1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][3]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[0][2]; tk ++){
					Eigen::Matrix<double,3,2> downUp;
					downUp.block(0,0,3,1) = wormProf(tj).block(0,0,3,1);
					downUp.block(0,1,3,1) = 
						(1.0 - (double)tj / dehwSurf.gridNumb[0][3]) 
						* lineCoor.block(0,0,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[0][3] 
						* lineCoor.block(0,1,3,1).eval();
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tk / dehwSurf.gridNumb[0][2]) 
						* downUp.block(0,0,3,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[0][2] * downUp.block(0,1,3,1).eval();
					blocPoin[2][ti][tj][tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
				for(long tk = 1; tk <= dehwSurf.gridNumb[0][2]; tk ++){
					Eigen::Matrix<double,3,2> downUp;
					downUp.block(0,0,3,1) = 
						(1.0 - (double)tj / dehwSurf.gridNumb[0][3]) 
						* lineCoor.block(0,0,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[0][3] * lineCoor.block(0,1,3,1).eval();
					downUp.block(0,1,3,1) = wormProf(tj).block(0,1,3,1);
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tk / dehwSurf.gridNumb[0][2]) 
						* downUp.block(0,0,3,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[0][2] * downUp.block(0,1,3,1).eval();
					blocPoin[2][ti][tj][dehwSurf.gridNumb[0][2] + tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//3
			lineCoor.block(0,1,3,1) = blocCoor.block(0,2,3,1).eval();
			for(long tj = 0; tj <= dehwSurf.gridNumb[0][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[0][0] / 2; tk ++){
					Eigen::Matrix<double,3,2> downUp;
					downUp.block(0,0,3,1) = 
						(1.0 - (double)tk / (dehwSurf.gridNumb[0][0] / 2.0)) 
						* lineCoor.block(0,0,3,1).eval() 
						+ (double)tk / (dehwSurf.gridNumb[0][0]/2.0) 
						* lineCoor.block(0,1,3,1).eval();
					downUp.block(0,1,3,1) = 
						rootProf_2.block(0,dehwSurf.gridNumb[0][0] / 2 - tk,3,1);
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tj / dehwSurf.gridNumb[0][2]) 
						* downUp.block(0,0,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[0][2] * downUp.block(0,1,3,1).eval();
					COOR inpuCoor = {tempCoor(0), tempCoor(1), tempCoor(2)};
					blocPoin[3][ti][tj][tk] = multGrid[tw].TRY_ADD_NODE(inpuCoor);
					if(isnoDode == 1 && tk == dehwSurf.gridNumb[0][0] / 2){
						wodeAucu[tw][1].INSERT(
							(dehwSurf.gridNumb[0][1] + tj) * wodeFact_0, 
							ti * wodeFact_1, inpuCoor
						);
					}
				}
			}
		}
		//volumes
		for(long tj = 0; tj < blocPoin.size(); tj ++){
			for(long tk = 0; tk < blocPoin[tj].size() - 1; tk ++){
				for(long tm = 0; tm < blocPoin[tj][tk].size() - 1; tm ++){
					for(long tn = 0; tn < blocPoin[tj][tk][tm].size() - 1; tn ++){
						TREE_ELEM tempElem;
						tempElem.parent = -1;
						tempElem.cornNode[0] = blocPoin[tj][tk][tm][tn];
						tempElem.cornNode[1] = blocPoin[tj][tk][tm + 1][tn];
						tempElem.cornNode[2] = blocPoin[tj][tk][tm + 1][tn + 1];
						tempElem.cornNode[3] = blocPoin[tj][tk][tm][tn + 1];
						tempElem.cornNode[4] = blocPoin[tj][tk + 1][tm][tn];
						tempElem.cornNode[5] = blocPoin[tj][tk + 1][tm + 1][tn];
						tempElem.cornNode[6] = blocPoin[tj][tk + 1][tm + 1][tn + 1];
						tempElem.cornNode[7] = blocPoin[tj][tk + 1][tm][tn + 1];
						tempElem.level = 0;
						tempElem.refiPatt = 7;
						tempElem.children.resize(0);
						long elemNumb = multGrid[tw].ADD_ELEMENT(tempElem);
					}
				}
			}
		}
		//global refinement
		std::set<long> spliElem;
		std::map<long,std::set<long>> spliFlag;
		std::map<std::vector<long>,COOR> planSurf;
		for(long tr = 0; tr < dehwSurf.globInho + dehwSurf.globHomo; tr ++){
			spliElem.clear();
			for(long ti = 0; ti < multGrid[tw].elemVect.size(); ti ++){
				if(multGrid[tw].elemVect[ti].children.size() > 0){
					continue;
				}
				spliElem.insert(ti);
				if(tr < dehwSurf.globInho){
					multGrid[tw].elemVect[ti].refiPatt = 6;
				}
				else{
					multGrid[tw].elemVect[ti].refiPatt = 0;
				}
			}
			planSurf.clear();
			std::vector<CURVEDS*> wormSurf = {
				&(dehwSurf.wormTosu), &(dehwSurf.wormToba), 
				&(dehwSurf.wormRtsu), &(dehwSurf.wormRtba)
			};
			for(long ts = 0; ts < wormSurf.size(); ts ++){
				(* wormSurf[ts]).REFINE(spliElem, multGrid[tw].elemVect, 
					multGrid[tw].nodeCoor, planSurf
				);
			}
			for(const auto &iterSpel : spliElem){
				//
				std::vector<long> inpuNode(8);
				std::vector<COOR> inpuCoor(8);
				for(long tk = 0; tk < 8; tk ++){
					inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[tk];
					auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
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
						inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[hexaLine[tj][tk]];
						auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
						inpuCoor[tk] = iterNoco->second;
					}
					COOR outpCoor;
					COOR_AVER(inpuCoor, outpCoor);
					sort(inpuNode.begin(), inpuNode.end());
					//the same key: new key will be abandoned
					planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
						inpuNode, outpCoor
					));
				}
				//
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					std::vector<long> inpuNode(hexaFace[tj].size());
					std::vector<COOR> inpuCoor(hexaFace[tj].size());
					for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
						inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
						inpuCoor[tk] = iterNoco->second;
					}
					COOR outpCoor;
					COOR_AVER(inpuCoor, outpCoor);
					sort(inpuNode.begin(), inpuNode.end());
					//the same key: new key will be abandoned
					planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
						inpuNode, outpCoor
					));
				}
			}
			multGrid[tw].REFINE(spliElem, spliFlag, planSurf);
			if(isnoDode == 1){
				for(const auto &iterPlsu : planSurf){
					UPDA_WODE(iterPlsu.first, iterPlsu.second, tw, 0, tw);
					UPDA_WODE(iterPlsu.first, iterPlsu.second, tw, 1, tw);
				}
			}
		}
		//
		multGrid[tw].RIGI_ROTR(wormRota_1, wormTran);
		if(isnoDode == 1 && isnoDode == 1){
			wodeAucu[tw][0].RIGI_ROTR(wormRota_1, wormTran);
			wodeAucu[tw][1].RIGI_ROTR(wormRota_1, wormTran);
		}
	}
	return 1;
}

long DEHW_1::WHEE_MESH_DD(){
	OUTPUT_TIME("DEHW_1::WHEE_MESH_DD");
	//
	double wheeAngl_1 = analAngl[1] - 2.0 * PI / dehwSurf.z[1] * 2.0;
	Eigen::Matrix3d wheeRota;
	wheeRota << cos(wheeAngl_1), -sin(wheeAngl_1), 0.0, 
		sin(wheeAngl_1), cos(wheeAngl_1), 0.0, 
		0.0, 0.0, 1.0;
	Eigen::Vector3d wheeTran;
	wheeTran << 0.0, 0.0, 0.0;
	long numbFace = dehwSurf.gridNumb[1][4] / dehwSurf.gridNumb[1][6];
	long whdeFact_0 = (1 << (dehwSurf.globHomo));
	long whdeFact_1 = (1 << (dehwSurf.globInho + dehwSurf.globHomo));
	whdeAucu.resize(dehwSurf.gridNumb[1][6]);
	for(long ti = 0; ti < dehwSurf.gridNumb[1][6]; ti ++){
		VECT_RESI(whdeAucu[ti].indiPoin, 
			(dehwSurf.gridNumb[1][1] + dehwSurf.gridNumb[1][2]) * whdeFact_0 + 1, 
			numbFace * whdeFact_1 + 1
		);
	}
	//
	#pragma omp parallel for
	for(long tw = dehwSurf.gridNumb[0][6]; tw < multGrid.size(); tw ++){
		long toot_tw = (tw - dehwSurf.gridNumb[0][6]) / dehwSurf.gridNumb[1][6];
		long face_tw = (tw - dehwSurf.gridNumb[0][6]) % dehwSurf.gridNumb[1][6];
		//point
		VECTOR4L blocPoin;
		VECT_RESI(blocPoin, 4, numbFace + 1);
		for(long ti = 0; ti <= numbFace; ti ++){
			VECT_RESI(blocPoin[0][ti], dehwSurf.gridNumb[1][1] + 1, 
				dehwSurf.gridNumb[1][0] + 1
			);
			VECT_RESI(blocPoin[1][ti], dehwSurf.gridNumb[1][2] + 1, 
				dehwSurf.gridNumb[1][0] / 2 + 1
			);
			VECT_RESI(blocPoin[2][ti], dehwSurf.gridNumb[1][3] + 1, 
				2 * dehwSurf.gridNumb[1][2] + 1
			);
			VECT_RESI(blocPoin[3][ti], dehwSurf.gridNumb[1][2] + 1, 
				dehwSurf.gridNumb[1][0] / 2 + 1
			);
			long ti_real = face_tw * numbFace + ti;
			ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
			Eigen::Vector<Eigen::Matrix<double,3,2>,Eigen::Dynamic> wheeProf;
			wheeProf.resize(dehwSurf.gridNumb[1][3] + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
				wheeProf(tj).block(0,0,3,1) << dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][2];
				wheeProf(tj).block(0,1,3,1) << dehwSurf.wheeToba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeToba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeToba.indiPoin[ti_real][tj_real][2];
			}
			//transition into "unfolded cone surface"
			double tempAlph_3 = - dehwSurf.curvCoor[1][ti_real][0](0);
			Eigen::Matrix<Eigen::Vector2d,Eigen::Dynamic,Eigen::Dynamic> wheeProf_;
			wheeProf_.resize(dehwSurf.gridNumb[1][3] + 1, 2);
			for(long tk = 0; tk <= 1; tk ++){
				for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
					long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
					Eigen::Vector3d tempCoor;
					if(tk == 0){
						tempCoor << dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][0], 
							dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][1], 
							dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][2];
					}
					else{
						tempCoor << dehwSurf.wheeToba.indiPoin[ti_real][tj_real][0], 
							dehwSurf.wheeToba.indiPoin[ti_real][tj_real][1], 
							dehwSurf.wheeToba.indiPoin[ti_real][tj_real][2];
					}
					wheeProf_(tj,tk)  = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
				}
			}
			ti_real = face_tw * numbFace + ti;
			ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo));
			double r_1f = dehwSurf.a_h2 / cos(tempAlph_3) 
				- (dehwSurf.a_h2 - dehwSurf.d_f[1] / 2.0);
			Eigen::MatrixXd rootProf_0, rootProf_1;
			rootProf_0.resize(2, dehwSurf.gridNumb[1][0] / 2 + 1);
			rootProf_1.resize(2, dehwSurf.gridNumb[1][0] / 2 + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][0] / 2; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo));
				Eigen::Vector3d tempCoor;
				tempCoor << dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][2];
				rootProf_0.block(0,tj,2,1) = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
				tempCoor << dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][2];
				rootProf_1.block(0,tj,2,1) = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
			}
			double tranRadi_0 = r_1f - PI / 4.0 * dehwSurf.m_t;
			Eigen::Vector2d tempAngl_0;
			tempAngl_0(0) = atan2(rootProf_0(1,0), rootProf_0(0,0));
			tempAngl_0(1) = atan2(rootProf_1(1,0), rootProf_1(0,0));
			Eigen::MatrixXd tranProf_0, tranProf_1, inneProf;
			tranProf_0.resize(2, dehwSurf.gridNumb[1][0] + 1);
			tranProf_1.resize(3, dehwSurf.gridNumb[1][0] + 1);
			inneProf.resize(3, dehwSurf.gridNumb[1][0] + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][0]; tj ++){
				double tempAngl_j = tempAngl_0(0) 
					+ (tempAngl_0(1) - tempAngl_0(0)) / dehwSurf.gridNumb[1][0] * tj;
				tranProf_0.block(0,tj,2,1) << tranRadi_0 * cos(tempAngl_j), 
					tranRadi_0 * sin(tempAngl_j);
				tranProf_1.block(0,tj,3,1) = 
					dehwSurf.WHEE_CONE(tranProf_0.block(0,tj,2,1), tempAlph_3);
				double tempRadi = sqrt(pow(tranProf_1(0,tj), 2.0) + pow(tranProf_1(1,tj), 2.0));
				inneProf.block(0,tj,3,1) << dehwSurf.inneRadi[1] / tempRadi * tranProf_1(0,tj), 
					dehwSurf.inneRadi[1] / tempRadi * tranProf_1(1,tj), 
					tranProf_1(2,tj);
			}
			//0
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][1]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0]; tk ++){
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tj / dehwSurf.gridNumb[1][1]) 
						* inneProf.block(0,tk,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][1] 
						* tranProf_1.block(0,tk,3,1).eval();
					blocPoin[0][ti][tj][tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
					if(toot_tw == 0 && tk == dehwSurf.gridNumb[1][0]){
						COOR inpuCoor = {tempCoor(0), tempCoor(1), tempCoor(2)};
						whdeAucu[face_tw].INSERT(tj * whdeFact_0, ti * whdeFact_1, inpuCoor);
					}
				}
			}
			//1
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
						* tranProf_0.block(0,tk,2,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][2] 
						* rootProf_0.block(0,tk,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					blocPoin[1][ti][tj][tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//2
			Eigen::MatrixXd lineCoor;
			lineCoor.resize(2, dehwSurf.gridNumb[1][3] + 1);
			lineCoor.block(0,0,2,1) = tranProf_0.block(0,dehwSurf.gridNumb[1][0] / 2,2,1);
			lineCoor.block(0,dehwSurf.gridNumb[1][3],2,1) = 
				0.5 * (wheeProf_(dehwSurf.gridNumb[1][3],0) 
				+ wheeProf_(dehwSurf.gridNumb[1][3],1));
			for(long tj = 1; tj < dehwSurf.gridNumb[1][3]; tj ++){
				lineCoor.block(0,tj,2,1) = lineCoor.block(0,0,2,1) 
					+ (lineCoor.block(0,dehwSurf.gridNumb[1][3],2,1) 
					- lineCoor.block(0,0,2,1)) / dehwSurf.gridNumb[1][3] * tj;
			}
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][2]; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tk / dehwSurf.gridNumb[1][2]) * wheeProf_(tj,0) 
						+ (double)tk / dehwSurf.gridNumb[1][2] * lineCoor.block(0,tj,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					blocPoin[2][ti][tj][tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
				for(long tk = 1; tk <= dehwSurf.gridNumb[1][2]; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tk / dehwSurf.gridNumb[1][2]) 
						* lineCoor.block(0,tj,2,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[1][2] * wheeProf_(tj,1);
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					blocPoin[2][ti][tj][dehwSurf.gridNumb[1][2] + tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//3
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
					Eigen::Vector2d tempPoin = (1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
						* tranProf_0.block(0,dehwSurf.gridNumb[1][0] / 2 + tk,2,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][2] 
						* rootProf_1.block(0,dehwSurf.gridNumb[1][0] / 2 - tk,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					blocPoin[3][ti][tj][tk] = 
						multGrid[tw].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
					if(toot_tw == 0 && tk == dehwSurf.gridNumb[1][0] / 2){
						COOR inpuCoor = {tempCoor(0), tempCoor(1), tempCoor(2)};
						whdeAucu[face_tw].INSERT(
							(dehwSurf.gridNumb[1][1] + tj) * whdeFact_0, 
							ti * whdeFact_1, inpuCoor
						);
					}
				}
			}
		}
		//volume
		for(long tj = 0; tj < blocPoin.size(); tj ++){
			for(long tk = 0; tk < blocPoin[tj].size() - 1; tk ++){
				for(long tm = 0; tm < blocPoin[tj][tk].size() - 1; tm ++){
					for(long tn = 0; tn < blocPoin[tj][tk][tm].size() - 1; tn ++){
						TREE_ELEM tempElem;
						tempElem.parent = -1;
						tempElem.cornNode[0] = blocPoin[tj][tk][tm][tn];
						tempElem.cornNode[1] = blocPoin[tj][tk][tm + 1][tn];
						tempElem.cornNode[2] = blocPoin[tj][tk][tm + 1][tn + 1];
						tempElem.cornNode[3] = blocPoin[tj][tk][tm][tn + 1];
						tempElem.cornNode[4] = blocPoin[tj][tk + 1][tm][tn];
						tempElem.cornNode[5] = blocPoin[tj][tk + 1][tm + 1][tn];
						tempElem.cornNode[6] = blocPoin[tj][tk + 1][tm + 1][tn + 1];
						tempElem.cornNode[7] = blocPoin[tj][tk + 1][tm][tn + 1];
						tempElem.level = 0;
						tempElem.refiPatt = 7;
						tempElem.children.resize(0);
						long elemNumb = multGrid[tw].ADD_ELEMENT(tempElem);
					}
				}
			}
		}
		//global refinement
		std::set<long> spliElem;
		std::map<long,std::set<long>> spliFlag;
		std::map<std::vector<long>,COOR> planSurf;
		for(long tr = 0; tr < dehwSurf.globInho + dehwSurf.globHomo; tr ++){
			spliElem.clear();
			for(long ti = 0; ti < multGrid[tw].elemVect.size(); ti ++){
				if(multGrid[tw].elemVect[ti].children.size() > 0){
					continue;
				}
				spliElem.insert(ti);
				if(tr < dehwSurf.globInho){
					multGrid[tw].elemVect[ti].refiPatt = 6;
				}
				else{
					multGrid[tw].elemVect[ti].refiPatt = 0;
				}
			}
			planSurf.clear();
			std::vector<CURVEDS*> wheeSurf = {
				&(dehwSurf.wheeTosu), &(dehwSurf.wheeToba), 
				&(dehwSurf.wheeRtsu), &(dehwSurf.wheeRtba)
			};
			for(long ts = 0; ts < wheeSurf.size(); ts ++){
				(* wheeSurf[ts]).REFINE(spliElem, multGrid[tw].elemVect, 
					multGrid[tw].nodeCoor, planSurf
				);
			}
			for(const auto &iterSpel : spliElem){
				//
				std::vector<long> inpuNode(8);
				std::vector<COOR> inpuCoor(8);
				for(long tk = 0; tk < 8; tk ++){
					inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[tk];
					auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
					inpuCoor[tk] = iterNoco->second;
				}
				COOR outpCoor;
				COOR_AVER_1(inpuCoor, outpCoor);
				sort(inpuNode.begin(), inpuNode.end());
				//the same key: new key will be abandoned
				planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
					inpuNode, outpCoor
				));
				//
				for(long tj = 0; tj < hexaLine.size(); tj ++){
					std::vector<long> inpuNode(hexaLine[tj].size());
					std::vector<COOR> inpuCoor(hexaLine[tj].size());
					for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
						inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[hexaLine[tj][tk]];
						auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
						inpuCoor[tk] = iterNoco->second;
					}
					COOR outpCoor;
					COOR_AVER_1(inpuCoor, outpCoor);
					sort(inpuNode.begin(), inpuNode.end());
					//the same key: new key will be abandoned
					planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
						inpuNode, outpCoor
					));
				}
				//
				for(long tj = 0; tj < hexaFace.size(); tj ++){
					std::vector<long> inpuNode(hexaFace[tj].size());
					std::vector<COOR> inpuCoor(hexaFace[tj].size());
					for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
						inpuNode[tk] = multGrid[tw].elemVect[iterSpel].cornNode[hexaFace[tj][tk]];
						auto iterNoco = multGrid[tw].nodeCoor.find(inpuNode[tk]);
						inpuCoor[tk] = iterNoco->second;
					}
					COOR outpCoor;
					COOR_AVER_1(inpuCoor, outpCoor);
					sort(inpuNode.begin(), inpuNode.end());
					//the same key: new key will be abandoned
					planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
						inpuNode, outpCoor
					));
				}
			}
			multGrid[tw].REFINE(spliElem, spliFlag, planSurf);
			if(toot_tw  == 0){
				for(const auto &iterPlsu : planSurf){
					UPDA_WHDE(iterPlsu.first, iterPlsu.second, face_tw, tw);
				}
			}
		}
		//
		Eigen::Matrix3d rotaMatr;
		double tempAngl = 2.0 * PI / dehwSurf.z[1] * (double)toot_tw;
		rotaMatr << cos(tempAngl), -sin(tempAngl), 0.0,
			sin(tempAngl), cos(tempAngl), 0.0,
			0.0, 0.0, 1.0;
		Eigen::Vector3d tranVect = Eigen::Vector3d::Zero(3,1);
		multGrid[tw].RIGI_ROTR(rotaMatr, tranVect);
		//
		multGrid[tw].RIGI_ROTR(wheeRota, wheeTran);
		if(toot_tw  == 0){
			whdeAucu[face_tw].RIGI_ROTR(rotaMatr, tranVect);
			whdeAucu[face_tw].RIGI_ROTR(wheeRota, wheeTran);
		}
	}
	return 1;
}

long DEHW_1::WHEE_MESH_NODD(){
	OUTPUT_TIME("DEHW_1::WHEE_MESH_NODD");
	//
	long numbFace = dehwSurf.gridNumb[1][4];
	long tw_real = 1;
	for(long tw = 0; tw < dehwSurf.gridNumb[1][5]; tw ++){
		long toot_tw = tw;
		long face_tw = 0;
		Eigen::Matrix3d rotaMatr;
		double tempAngl = 2.0 * PI / dehwSurf.z[1] * (double)toot_tw;
		rotaMatr << cos(tempAngl), -sin(tempAngl), 0.0,
			sin(tempAngl), cos(tempAngl), 0.0,
			0.0, 0.0, 1.0;
		//point
		VECTOR4L blocPoin;
		VECT_RESI(blocPoin, 4, numbFace + 1);
		for(long ti = 0; ti <= numbFace; ti ++){
			VECT_RESI(blocPoin[0][ti], dehwSurf.gridNumb[1][1] + 1, 
				dehwSurf.gridNumb[1][0] + 1
			);
			VECT_RESI(blocPoin[1][ti], dehwSurf.gridNumb[1][2] + 1, 
				dehwSurf.gridNumb[1][0] / 2 + 1
			);
			VECT_RESI(blocPoin[2][ti], dehwSurf.gridNumb[1][3] + 1, 
				2 * dehwSurf.gridNumb[1][2] + 1
			);
			VECT_RESI(blocPoin[3][ti], dehwSurf.gridNumb[1][2] + 1, 
				dehwSurf.gridNumb[1][0] / 2 + 1
			);
			long ti_real = face_tw * numbFace + ti;
			ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve));
			Eigen::Vector<Eigen::Matrix<double,3,2>,Eigen::Dynamic> wheeProf;
			wheeProf.resize(dehwSurf.gridNumb[1][3] + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
				wheeProf(tj).block(0,0,3,1) << dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][2];
				wheeProf(tj).block(0,1,3,1) << dehwSurf.wheeToba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeToba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeToba.indiPoin[ti_real][tj_real][2];
			}
			//transition into "unfolded cone surface"
			double tempAlph_3 = - dehwSurf.curvCoor[1][ti_real][0](0);
			Eigen::Matrix<Eigen::Vector2d,Eigen::Dynamic,Eigen::Dynamic> wheeProf_;
			wheeProf_.resize(dehwSurf.gridNumb[1][3] + 1, 2);
			for(long tk = 0; tk <= 1; tk ++){
				for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
					long tj_real = tj * (1 << (dehwSurf.globHomo + dehwSurf.locaLeve));
					Eigen::Vector3d tempCoor;
					if(tk == 0){
						tempCoor << dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][0], 
							dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][1], 
							dehwSurf.wheeTosu.indiPoin[ti_real][tj_real][2];
					}
					else{
						tempCoor << dehwSurf.wheeToba.indiPoin[ti_real][tj_real][0], 
							dehwSurf.wheeToba.indiPoin[ti_real][tj_real][1], 
							dehwSurf.wheeToba.indiPoin[ti_real][tj_real][2];
					}
					wheeProf_(tj,tk)  = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
				}
			}
			ti_real = face_tw * numbFace + ti;
			ti_real *= (1 << (dehwSurf.globInho + dehwSurf.globHomo));
			double r_1f = dehwSurf.a_h2 / cos(tempAlph_3) 
				- (dehwSurf.a_h2 - dehwSurf.d_f[1] / 2.0);
			Eigen::MatrixXd rootProf_0, rootProf_1;
			rootProf_0.resize(2, dehwSurf.gridNumb[1][0] / 2 + 1);
			rootProf_1.resize(2, dehwSurf.gridNumb[1][0] / 2 + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][0] / 2; tj ++){
				long tj_real = tj * (1 << (dehwSurf.globHomo));
				Eigen::Vector3d tempCoor;
				tempCoor << dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeRtsu.indiPoin[ti_real][tj_real][2];
				rootProf_0.block(0,tj,2,1) = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
				tempCoor << dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][0], 
					dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][1], 
					dehwSurf.wheeRtba.indiPoin[ti_real][tj_real][2];
				rootProf_1.block(0,tj,2,1) = dehwSurf.WHEE_UNCONE(tempCoor, tempAlph_3);
			}
			double tranRadi_0 = r_1f - PI / 4.0 * dehwSurf.m_t;
			Eigen::Vector2d tempAngl_0;
			tempAngl_0(0) = atan2(rootProf_0(1,0), rootProf_0(0,0));
			tempAngl_0(1) = atan2(rootProf_1(1,0), rootProf_1(0,0));
			Eigen::MatrixXd tranProf_0, tranProf_1, inneProf;
			tranProf_0.resize(2, dehwSurf.gridNumb[1][0] + 1);
			tranProf_1.resize(3, dehwSurf.gridNumb[1][0] + 1);
			inneProf.resize(3, dehwSurf.gridNumb[1][0] + 1);
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][0]; tj ++){
				double tempAngl_j = tempAngl_0(0) 
					+ (tempAngl_0(1) - tempAngl_0(0)) / dehwSurf.gridNumb[1][0] * tj;
				tranProf_0.block(0,tj,2,1) << tranRadi_0 * cos(tempAngl_j), 
					tranRadi_0 * sin(tempAngl_j);
				tranProf_1.block(0,tj,3,1) = 
					dehwSurf.WHEE_CONE(tranProf_0.block(0,tj,2,1), tempAlph_3);
				double tempRadi = sqrt(pow(tranProf_1(0,tj), 2.0) + pow(tranProf_1(1,tj), 2.0));
				inneProf.block(0,tj,3,1) << dehwSurf.inneRadi[1] / tempRadi * tranProf_1(0,tj), 
					dehwSurf.inneRadi[1] / tempRadi * tranProf_1(1,tj), 
					tranProf_1(2,tj);
			}
			//0
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][1]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0]; tk ++){
					Eigen::Vector3d tempCoor = 
						(1.0 - (double)tj / dehwSurf.gridNumb[1][1]) 
						* inneProf.block(0,tk,3,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][1] 
						* tranProf_1.block(0,tk,3,1).eval();
					tempCoor = rotaMatr * tempCoor;
					blocPoin[0][ti][tj][tk] = 
						multGrid[tw_real].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//1
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
						* tranProf_0.block(0,tk,2,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][2] 
						* rootProf_0.block(0,tk,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					tempCoor = rotaMatr * tempCoor;
					blocPoin[1][ti][tj][tk] = 
						multGrid[tw_real].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//2
			Eigen::MatrixXd lineCoor;
			lineCoor.resize(2, dehwSurf.gridNumb[1][3] + 1);
			lineCoor.block(0,0,2,1) = tranProf_0.block(0,dehwSurf.gridNumb[1][0] / 2,2,1);
			lineCoor.block(0,dehwSurf.gridNumb[1][3],2,1) = 
				0.5 * (wheeProf_(dehwSurf.gridNumb[1][3],0) 
				+ wheeProf_(dehwSurf.gridNumb[1][3],1));
			for(long tj = 1; tj < dehwSurf.gridNumb[1][3]; tj ++){
				lineCoor.block(0,tj,2,1) = lineCoor.block(0,0,2,1) 
					+ (lineCoor.block(0,dehwSurf.gridNumb[1][3],2,1) 
					- lineCoor.block(0,0,2,1)) / dehwSurf.gridNumb[1][3] * tj;
			}
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][3]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][2]; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tk / dehwSurf.gridNumb[1][2]) * wheeProf_(tj,0) 
						+ (double)tk / dehwSurf.gridNumb[1][2] * lineCoor.block(0,tj,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					tempCoor = rotaMatr * tempCoor;
					blocPoin[2][ti][tj][tk] = 
						multGrid[tw_real].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
				for(long tk = 1; tk <= dehwSurf.gridNumb[1][2]; tk ++){
					Eigen::Vector2d tempPoin = 
						(1.0 - (double)tk / dehwSurf.gridNumb[1][2]) 
						* lineCoor.block(0,tj,2,1).eval() 
						+ (double)tk / dehwSurf.gridNumb[1][2] * wheeProf_(tj,1);
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					tempCoor = rotaMatr * tempCoor;
					blocPoin[2][ti][tj][dehwSurf.gridNumb[1][2] + tk] = 
						multGrid[tw_real].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
			//3
			for(long tj = 0; tj <= dehwSurf.gridNumb[1][2]; tj ++){
				for(long tk = 0; tk <= dehwSurf.gridNumb[1][0] / 2; tk ++){
					Eigen::Vector2d tempPoin = (1.0 - (double)tj / dehwSurf.gridNumb[1][2]) 
						* tranProf_0.block(0,dehwSurf.gridNumb[1][0] / 2 + tk,2,1).eval() 
						+ (double)tj / dehwSurf.gridNumb[1][2] 
						* rootProf_1.block(0,dehwSurf.gridNumb[1][0] / 2 - tk,2,1).eval();
					Eigen::Vector3d tempCoor = dehwSurf.WHEE_CONE(tempPoin, tempAlph_3);
					tempCoor = rotaMatr * tempCoor;
					blocPoin[3][ti][tj][tk] = 
						multGrid[tw_real].TRY_ADD_NODE({tempCoor(0), tempCoor(1), tempCoor(2)});
				}
			}
		}
		//volume
		for(long tj = 0; tj < blocPoin.size(); tj ++){
			for(long tk = 0; tk < blocPoin[tj].size() - 1; tk ++){
				for(long tm = 0; tm < blocPoin[tj][tk].size() - 1; tm ++){
					for(long tn = 0; tn < blocPoin[tj][tk][tm].size() - 1; tn ++){
						TREE_ELEM tempElem;
						tempElem.parent = -1;
						tempElem.cornNode[0] = blocPoin[tj][tk][tm][tn];
						tempElem.cornNode[1] = blocPoin[tj][tk][tm + 1][tn];
						tempElem.cornNode[2] = blocPoin[tj][tk][tm + 1][tn + 1];
						tempElem.cornNode[3] = blocPoin[tj][tk][tm][tn + 1];
						tempElem.cornNode[4] = blocPoin[tj][tk + 1][tm][tn];
						tempElem.cornNode[5] = blocPoin[tj][tk + 1][tm + 1][tn];
						tempElem.cornNode[6] = blocPoin[tj][tk + 1][tm + 1][tn + 1];
						tempElem.cornNode[7] = blocPoin[tj][tk + 1][tm][tn + 1];
						tempElem.level = 0;
						tempElem.refiPatt = 7;
						tempElem.children.resize(0);
						long elemNumb = multGrid[tw_real].ADD_ELEMENT(tempElem);
					}
				}
			}
		}
	}
	//global refinement
	std::set<long> spliElem;
	std::map<long,std::set<long>> spliFlag;
	std::map<std::vector<long>,COOR> planSurf;
	for(long tr = 0; tr < dehwSurf.globInho + dehwSurf.globHomo; tr ++){
		spliElem.clear();
		for(long ti = 0; ti < multGrid[tw_real].elemVect.size(); ti ++){
			if(multGrid[tw_real].elemVect[ti].children.size() > 0){
				continue;
			}
			spliElem.insert(ti);
			if(tr < dehwSurf.globInho){
				multGrid[tw_real].elemVect[ti].refiPatt = 6;
			}
			else{
				multGrid[tw_real].elemVect[ti].refiPatt = 0;
			}
		}
		planSurf.clear();
		for(long ti = 0; ti < dehwSurf.gridNumb[1][5]; ti ++){
			Eigen::Matrix3d wheeRota;
			double tempAngl = 2.0 * PI / dehwSurf.z[1] * (double)ti;
			wheeRota << cos(tempAngl), -sin(tempAngl), 0.0,
				sin(tempAngl), cos(tempAngl), 0.0,
				0.0, 0.0, 1.0;
			Eigen::Vector3d wheeTran;
			wheeTran << 0.0, 0.0, 0.0;
			CURVEDS tempSurf_0 = dehwSurf.wheeTosu;
			tempSurf_0.RIGI_ROTR(wheeRota, wheeTran);
			CURVEDS tempSurf_1 = dehwSurf.wheeToba;
			tempSurf_1.RIGI_ROTR(wheeRota, wheeTran);
			CURVEDS tempSurf_2 = dehwSurf.wheeRtsu;
			tempSurf_2.RIGI_ROTR(wheeRota, wheeTran);
			CURVEDS tempSurf_3 = dehwSurf.wheeRtba;
			tempSurf_3.RIGI_ROTR(wheeRota, wheeTran);
			std::vector<CURVEDS*> wheeSurf = {
				&tempSurf_0, &tempSurf_1, &tempSurf_2, &tempSurf_3
			};
			for(long ts = 0; ts < wheeSurf.size(); ts ++){
				(* wheeSurf[ts]).REFINE(spliElem, multGrid[tw_real].elemVect, 
					multGrid[tw_real].nodeCoor, planSurf
				);
			}
		}
		for(const auto &iterSpel : spliElem){
			//
			std::vector<long> inpuNode(8);
			std::vector<COOR> inpuCoor(8);
			for(long tk = 0; tk < 8; tk ++){
				inpuNode[tk] = multGrid[tw_real].elemVect[iterSpel].cornNode[tk];
				auto iterNoco = multGrid[tw_real].nodeCoor.find(inpuNode[tk]);
				inpuCoor[tk] = iterNoco->second;
			}
			COOR outpCoor;
			COOR_AVER_1(inpuCoor, outpCoor);
			sort(inpuNode.begin(), inpuNode.end());
			//the same key: new key will be abandoned
			planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
				inpuNode, outpCoor
			));
			//
			for(long tj = 0; tj < hexaLine.size(); tj ++){
				std::vector<long> inpuNode(hexaLine[tj].size());
				std::vector<COOR> inpuCoor(hexaLine[tj].size());
				for(long tk = 0; tk < hexaLine[tj].size(); tk ++){
					inpuNode[tk] = 
						multGrid[tw_real].elemVect[iterSpel].cornNode[hexaLine[tj][tk]];
					auto iterNoco = multGrid[tw_real].nodeCoor.find(inpuNode[tk]);
					inpuCoor[tk] = iterNoco->second;
				}
				COOR outpCoor;
				COOR_AVER_1(inpuCoor, outpCoor);
				sort(inpuNode.begin(), inpuNode.end());
				//the same key: new key will be abandoned
				planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
					inpuNode, outpCoor
				));
			}
			//
			for(long tj = 0; tj < hexaFace.size(); tj ++){
				std::vector<long> inpuNode(hexaFace[tj].size());
				std::vector<COOR> inpuCoor(hexaFace[tj].size());
				for(long tk = 0; tk < hexaFace[tj].size(); tk ++){
					inpuNode[tk] = 
						multGrid[tw_real].elemVect[iterSpel].cornNode[hexaFace[tj][tk]];
					auto iterNoco = multGrid[tw_real].nodeCoor.find(inpuNode[tk]);
					inpuCoor[tk] = iterNoco->second;
				}
				COOR outpCoor;
				COOR_AVER_1(inpuCoor, outpCoor);
				sort(inpuNode.begin(), inpuNode.end());
				//the same key: new key will be abandoned
				planSurf.insert(std::map<std::vector<long>,COOR>::value_type(
					inpuNode, outpCoor
				));
			}
		}
		multGrid[tw_real].REFINE(spliElem, spliFlag, planSurf);
	}
	//
	double wheeAngl_1 = analAngl[1] - 2.0 * PI / dehwSurf.z[1] * 2.0;
	Eigen::Matrix3d wheeRota;
	wheeRota << cos(wheeAngl_1), -sin(wheeAngl_1), 0.0, 
		sin(wheeAngl_1), cos(wheeAngl_1), 0.0, 
		0.0, 0.0, 1.0;
	Eigen::Vector3d wheeTran;
	wheeTran << 0.0, 0.0, 0.0;
	multGrid[tw_real].RIGI_ROTR(wheeRota, wheeTran);
	return 1;
}

long DEHW_1::UPDA_WODE(const std::vector<long> &inpuNode, 
	const COOR &outpCoor, long face_tw, long face_sub, long worm_tv){
	long ti_aver = 0;
	long tj_aver = 0;
	bool tempFlag = true;
	for(long ic = 0; ic < inpuNode.size(); ic ++){
		auto iterNoco = multGrid[worm_tv].nodeCoor.find(inpuNode[ic]);
		auto iterPoin = wodeAucu[face_tw][face_sub].poinIndi.find(iterNoco->second);
		if(iterPoin == wodeAucu[face_tw][face_sub].poinIndi.end()){
			tempFlag = false;
			break;
		}
		ti_aver += (iterPoin->second)[0];
		tj_aver += (iterPoin->second)[1];
	}
	if(tempFlag == false){
		return 1;
	}
	ti_aver /= inpuNode.size();
	tj_aver /= inpuNode.size();
	wodeAucu[face_tw][face_sub].INSERT(ti_aver, tj_aver, outpCoor);
	return 1;
}

long DEHW_1::UPDA_WHDE(const std::vector<long> &inpuNode, 
	const COOR &outpCoor, long face_tw, long whee_tv){
	long ti_aver = 0;
	long tj_aver = 0;
	bool tempFlag = true;
	for(long ic = 0; ic < inpuNode.size(); ic ++){
		auto iterNoco = multGrid[whee_tv].nodeCoor.find(inpuNode[ic]);
		auto iterPoin = whdeAucu[face_tw].poinIndi.find(iterNoco->second);
		if(iterPoin == whdeAucu[face_tw].poinIndi.end()){
			tempFlag = false;
			break;
		}
		ti_aver += (iterPoin->second)[0];
		tj_aver += (iterPoin->second)[1];
	}
	if(tempFlag == false){
		return 1;
	}
	ti_aver /= inpuNode.size();
	tj_aver /= inpuNode.size();
	whdeAucu[face_tw].INSERT(ti_aver, tj_aver, outpCoor);
	return 1;
}

long DEHW_1::CONT_INTE_DD(){
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD local mesh refinement");
	auto CART_CURV = [=](COOR tempCoor, double &tempXico, double &tempEtac){
		tempXico += tempCoor[1];
		tempEtac += sqrt(pow(tempCoor[0] + (dehwSurf.a_h2 + centErro), 2.0)
			+ pow(tempCoor[2], 2.0));
	};
	std::vector<std::vector<bool>> isnoRefi(4, 
		std::vector<bool>((dehwSurf.circNumb == 8) ? 4 : 3)
	);
	std::vector<std::vector<CURVEDS>> cusuTabl_0(4, std::vector<CURVEDS>(2));
	VECTOR3L contBody_0(4, 
		std::vector<std::vector<long>>(isnoRefi[0].size(), std::vector<long>(2))
	);
	#pragma omp parallel for
	for(long tt = 0; tt < 4; tt ++){
		if(isnoRefi[tt].size() == 4){
			contBody_0[tt][0] = {4 + 8 * tt, dehwSurf.gridNumb[0][6] + 7 + 2 * tt};
			contBody_0[tt][1] = {3 + 8 * tt, dehwSurf.gridNumb[0][6] + 7 + 2 * tt};
			contBody_0[tt][2] = {3 + 8 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
			contBody_0[tt][3] = {2 + 8 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
		}
		else{
			contBody_0[tt][0] = {2 + 4 * tt, dehwSurf.gridNumb[0][6] + 7 + 2 * tt};
			contBody_0[tt][1] = {2 + 4 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
			contBody_0[tt][2] = {1 + 4 * tt, dehwSurf.gridNumb[0][6] + 6 + 2 * tt};
		}
		//
		Eigen::Matrix3d wormRota_0;
		wormRota_0 << cos(analAngl[0]), -sin(analAngl[0]), 0.0, 
			sin(analAngl[0]), cos(analAngl[0]), 0.0, 
			0.0, 0.0, 1.0;
		Eigen::Matrix3d wormRota_1;
		wormRota_1 << 1.0, 0.0, 0.0, 
			0.0, 0.0, 1.0, 
			0.0, -1.0, 0.0;
		wormRota_1 = wormRota_1 * wormRota_0;
		Eigen::Vector3d wormTran;
		wormTran << - (dehwSurf.a_h2 + centErro), 0.0, 0.0;
		CURVEDS mastSurf = dehwSurf.wormTosu;
		mastSurf.RIGI_ROTR(wormRota_1, wormTran);
		double wheeAngl_1 = analAngl[1] + 2.0 * PI / dehwSurf.z[1] * (1.0 + tt);
		Eigen::Matrix3d wheeRota;
		wheeRota << cos(wheeAngl_1), -sin(wheeAngl_1), 0.0, 
			sin(wheeAngl_1), cos(wheeAngl_1), 0.0, 
			0.0, 0.0, 1.0;
		Eigen::Vector3d wheeTran;
		wheeTran << 0.0, 0.0, 0.0;
		CURVEDS slavSurf = dehwSurf.wheeTosu;
		slavSurf.RIGI_ROTR(wheeRota, wheeTran);
		cusuTabl_0[tt][0] = mastSurf;
		cusuTabl_0[tt][1] = slavSurf;
		//
		for(long tr = 0; tr < dehwSurf.locaLeve; tr ++){
			long buckFact = (1 << 
				(dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve - 1 
				- (dehwSurf.locaLeve - tr))
			);
			std::vector<long> buckNumb = {dehwSurf.gridNumb[0][4] * buckFact, 
				dehwSurf.gridNumb[0][3] * buckFact
			};
			for(long tc = 0; tc < isnoRefi[tt].size(); tc ++){
				CSEARCH searCont_0;
				bool tempFlag;
				searCont_0.ADAPTIVE_REFINE(
					&(multGrid[contBody_0[tt][tc][0]]), &(multGrid[contBody_0[tt][tc][1]]), 
					tempFlag, &mastSurf, &slavSurf, 
					dehwSurf.globInho + dehwSurf.globHomo + tr, 
					distCrit[tr], buckNumb, CART_CURV
				);
				isnoRefi[tt][tc] = tempFlag;
			}
		}
	}
	loadIncr.resize(2);
	SUBR_COLO_WORM(-1);
	SUBR_COLO_WHEE(-1);
	#pragma omp parallel for
	for(long tw = 0; tw < multGrid.size(); tw ++){
		if(tw < dehwSurf.gridNumb[0][6]){
			SUBR_COLO_WORM(tw);
		}
		else{
			SUBR_COLO_WHEE(tw);
		}
	}
	#pragma omp parallel for
	for(long tw = 0; tw < multGrid.size(); tw ++){
		multGrid[tw].OUTPUT_ELEMENT(tw);
	}
	//****************************************************************************************
	long totaSear = 0;
	for(long tt = 0; tt < 4; tt ++){
		for(long tc = 0; tc < isnoRefi[tt].size(); tc ++){
			if(isnoRefi[tt][tc] == true){
				totaSear ++;
			}
		}
	}
	std::vector<std::vector<CURVEDS>> cusuTabl(totaSear, std::vector<CURVEDS>(2));
	totaSear += dehwSurf.gridNumb[0][6] - 1;
	totaSear += dehwSurf.gridNumb[0][6] - dehwSurf.circNumb;
	totaSear += dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1);
	totaSear += (dehwSurf.gridNumb[1][5] - 1) * dehwSurf.gridNumb[1][6];
	searCont.resize(totaSear);
	penaFact_n.resize(totaSear);
	penaFact_f.resize(totaSear);
	fricCoef.resize(totaSear);
	contBody.resize(totaSear);
	long searCoun = 0;
	for(long tt = 0; tt < 4; tt ++){
		for(long tc = 0; tc < isnoRefi[tt].size(); tc ++){
			if(isnoRefi[tt][tc] == true){
				contBody[searCoun] = contBody_0[tt][tc];
				searCont[searCoun].mastGrid = &(multGrid[contBody[searCoun][0]]);
				searCont[searCoun].slavGrid = &(multGrid[contBody[searCoun][1]]);
				penaFact_n[searCoun] = (multGrid[contBody[searCoun][0]].mateElas 
					+ multGrid[contBody[searCoun][1]].mateElas) * 500.0;
				penaFact_f[searCoun] = (multGrid[contBody[searCoun][0]].mateElas 
					+ multGrid[contBody[searCoun][1]].mateElas) * 
					((coloSett == 1) ? 5.0 : 50.0);
				fricCoef[searCoun] = (coloSett == 1) ? 0.08 : 0.2;
				cusuTabl[searCoun][0] = cusuTabl_0[tt][0];
				cusuTabl[searCoun][1] = cusuTabl_0[tt][1];
				searCoun ++;
			}
		}
	}
	cusuTabl_0.clear();
	for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - 1; tv ++){
		long ts = searCoun + tv;
		contBody[ts] = {tv, tv + 1};
		searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
		searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
		penaFact_n[ts] = (multGrid[contBody[ts][0]].mateElas 
			+ multGrid[contBody[ts][1]].mateElas) * 500.0;
		penaFact_f[ts] = (multGrid[contBody[ts][0]].mateElas 
			+ multGrid[contBody[ts][1]].mateElas) * 500.0;
		fricCoef[ts] = -1.0;
	}
	for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - dehwSurf.circNumb; tv ++){
		long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 + tv;
		contBody[ts] = {tv, tv + dehwSurf.circNumb};
		searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
		searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
		penaFact_n[ts] = (multGrid[contBody[ts][0]].mateElas 
			+ multGrid[contBody[ts][1]].mateElas) * 500.0;
		penaFact_f[ts] = (multGrid[contBody[ts][0]].mateElas 
			+ multGrid[contBody[ts][1]].mateElas) * 500.0;
		fricCoef[ts] = -1.0;
	}
	for(long ti = 0; ti < dehwSurf.gridNumb[1][5]; ti ++){
		for(long tj = 0; tj < dehwSurf.gridNumb[1][6] - 1; tj ++){
			long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
				+ dehwSurf.gridNumb[0][6] - dehwSurf.circNumb
				+ ti * (dehwSurf.gridNumb[1][6] - 1) + tj;
			long tv_0 = dehwSurf.gridNumb[0][6] + ti * dehwSurf.gridNumb[1][6] + tj;
			long tv_1 = tv_0 + 1;
			contBody[ts] = {tv_0, tv_1};
			searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
			searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
			penaFact_n[ts] = (multGrid[contBody[ts][0]].mateElas 
				+ multGrid[contBody[ts][1]].mateElas) * 500.0;
			penaFact_f[ts] = (multGrid[contBody[ts][0]].mateElas 
				+ multGrid[contBody[ts][1]].mateElas) * 500.0;
			fricCoef[ts] = -1.0;
		}
	}
	for(long ti = 0; ti < dehwSurf.gridNumb[1][5] - 1; ti ++){
		for(long tj = 0; tj < dehwSurf.gridNumb[1][6]; tj ++){
			long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
				+ dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
				+ dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1) 
				+ ti * dehwSurf.gridNumb[1][6] + tj;
			long tv_0 = dehwSurf.gridNumb[0][6] + ti * dehwSurf.gridNumb[1][6] + tj;
			long tv_1 = tv_0 + dehwSurf.gridNumb[1][6];
			contBody[ts] = {tv_0, tv_1};
			searCont[ts].mastGrid = &(multGrid[contBody[ts][0]]);
			searCont[ts].slavGrid = &(multGrid[contBody[ts][1]]);
			penaFact_n[ts] = (multGrid[contBody[ts][0]].mateElas 
				+ multGrid[contBody[ts][1]].mateElas) * 500.0;
			penaFact_f[ts] = (multGrid[contBody[ts][0]].mateElas 
				+ multGrid[contBody[ts][1]].mateElas) * 500.0;
			fricCoef[ts] = -1.0;
		}
	}
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD contact surface");
	long totaLeve = dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve;
	std::vector<long> buckNumb = {dehwSurf.gridNumb[0][4] * (1 << (totaLeve - 1)), 
		dehwSurf.gridNumb[0][3] * (1 << (totaLeve - 1))
	};
	#pragma omp parallel for
	for(long ts = 0; ts < searCoun; ts ++){
		EFACE_SURFACE iterEfsu_0(searCont[ts].mastGrid, &(cusuTabl[ts][0]));
		while(iterEfsu_0.INCREMENT() == 1){
			//NO NEED: if(multGrid[contBody[ts][0]].elemVect[iterEfsu_0.eid].level == totaLeve){
			searCont[ts].mastSegm.emplace_back(iterEfsu_0.currNode);
		}
		EFACE_SURFACE iterEfsu_1(searCont[ts].slavGrid, &(cusuTabl[ts][1]));
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
				auto iterNoco = multGrid[contBody[ts][0]].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				auto iterNoco = multGrid[contBody[ts][1]].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor, distCrit[dehwSurf.locaLeve - 1]);
		searCont[ts].OUTPUT_INPO(ts);
	}
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD surfaces between domains of worm tooth");
	Eigen::Matrix3d wormRota_0;
	wormRota_0 << cos(analAngl[0]), sin(analAngl[0]), 0.0, 
		-sin(analAngl[0]), cos(analAngl[0]), 0.0, 
		0.0, 0.0, 1.0;
	Eigen::Matrix3d wormRota_1;
	wormRota_1 << 1.0, 0.0, 0.0, 
		0.0, 0.0, -1.0, 
		0.0, 1.0, 0.0;
	wormRota_0 = wormRota_0 * wormRota_1;
	Eigen::Vector3d wormTran;
	wormTran << + (dehwSurf.a_h2 + centErro), 0.0, 0.0;
	auto CART_CURV_1 = [=](COOR tempCoor, double tempAngl) -> bool {
		Eigen::Vector3d coorVect;
		coorVect << tempCoor[0], tempCoor[1], tempCoor[2];
		coorVect = wormRota_0 * (coorVect + wormTran);
		double coorAngl = atan2(coorVect(1), coorVect(0));
		//!!!plus, but not minus!!!
		if(abs(coorAngl + tempAngl) < 1.0E-10){
			return true;
		}
		else{
			return false;
		}
	};
	buckNumb = {dehwSurf.gridNumb[0][0] * (1 << (dehwSurf.globHomo - 1)), 
		dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo))
	};
	#pragma omp parallel for
	for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - 1; tv ++){
		long ts = searCoun + tv;
		long tg_mast = contBody[ts][0];
		long tg_slav = contBody[ts][1];
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
					if(CART_CURV_1(iterNoco->second, wodeAuan[tv]) == false){
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
					if(CART_CURV_1(iterNoco->second, wodeAuan[tv]) == false){
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
				auto iterNoco = multGrid[tg_mast].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				auto iterNoco = multGrid[tg_slav].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
	}
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD surfaces between worm tooth");
	#pragma omp parallel for
	for(long tv = 0; tv < dehwSurf.gridNumb[0][6] - dehwSurf.circNumb; tv ++){
		if(tv == 0){
			buckNumb = {dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo)), 
				dehwSurf.gridNumb[0][5] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))
			};
		}
		else{
			buckNumb = {dehwSurf.gridNumb[0][1] * (1 << (dehwSurf.globHomo)), 
				dehwSurf.gridNumb[0][4] * (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))
			};
		}
		//
		long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 + tv;
		long tg_mast = contBody[ts][0];
		long tg_slav = contBody[ts][1];
		EFACE_SURFACE iterEfsu_0(&(multGrid[tg_mast]), &(wodeAucu[tg_mast][0]));
		while(iterEfsu_0.INCREMENT() == 1){
			searCont[ts].mastSegm.emplace_back(iterEfsu_0.currNode);
		}
		EFACE_SURFACE iterEfsu_1(&(multGrid[tg_slav]), &(wodeAucu[tg_slav][1]));
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
				auto iterNoco = multGrid[tg_mast].nodeCoor.find(node_tj);
				double tempX = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
				double tempY = - (iterNoco->second)[2];
				tempXico += sqrt(pow(tempX, 2.0) + pow(tempY, 2.0));
				tempEtac += (iterNoco->second)[1];
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				auto iterNoco = multGrid[tg_slav].nodeCoor.find(node_tj);
				double tempX = (iterNoco->second)[0] + (dehwSurf.a_h2 + centErro);
				double tempY = - (iterNoco->second)[2];
				tempXico += sqrt(pow(tempX, 2.0) + pow(tempY, 2.0));
				tempEtac += (iterNoco->second)[1];
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor);
		searCont[ts].OUTPUT_INPO(ts);
	}
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD surfaces between domains of wheel tooth");
	buckNumb = {dehwSurf.gridNumb[1][0] * (1 << (dehwSurf.globHomo - 1)), 
		dehwSurf.gridNumb[1][1] * (1 << (dehwSurf.globHomo))
	};
	#pragma omp parallel for
	for(long tv_0 = 0; tv_0 < dehwSurf.gridNumb[1][5]; tv_0 ++){
		for(long tv_1 = 0; tv_1 < dehwSurf.gridNumb[1][6] - 1; tv_1 ++){
			long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
				+ dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
				+ tv_0 * (dehwSurf.gridNumb[1][6] - 1) + tv_1;
			long tg_mast = contBody[ts][0];
			long tg_slav = contBody[ts][1];
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
						if(abs((iterNoco->second)[2]) > 1.0E-10){
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
						if(abs((iterNoco->second)[2]) > 1.0E-10){
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
					auto iterNoco = multGrid[tg_mast].nodeCoor.find(node_tj);
					tempXico += sqrt(pow((iterNoco->second)[0], 2.0) 
						+ pow((iterNoco->second)[1], 2.0));
					tempEtac += atan2((iterNoco->second)[1], (iterNoco->second)[0]);
				}
				mastCoor[0].emplace_back(tempXico / 4.0);
				mastCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
			VECTOR2D slavCoor(2);
			for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
				double tempXico = 0.0;
				double tempEtac = 0.0;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = searCont[ts].slavSegm[ti][tj];
					auto iterNoco = multGrid[tg_slav].nodeCoor.find(node_tj);
					tempXico += sqrt(pow((iterNoco->second)[0], 2.0) 
						+ pow((iterNoco->second)[1], 2.0));
					tempEtac += atan2((iterNoco->second)[1], (iterNoco->second)[0]);
				}
				slavCoor[0].emplace_back(tempXico / 4.0);
				slavCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].CONTACT_SEARCH(slavCoor);
			searCont[ts].OUTPUT_INPO(ts);
		}
	}
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_DD surfaces between wheel teeth");
	buckNumb = {dehwSurf.gridNumb[1][1] * (1 << (dehwSurf.globHomo)), 
		dehwSurf.gridNumb[1][4] / dehwSurf.gridNumb[1][6] 
		* (1 << (dehwSurf.globInho + dehwSurf.globHomo - 1))
	};
	#pragma omp parallel for
	for(long tv_0 = 0; tv_0 < dehwSurf.gridNumb[1][5] - 1; tv_0 ++){
		for(long tv_1 = 0; tv_1 < dehwSurf.gridNumb[1][6]; tv_1 ++){
			//
			CURVEDS tempCusu = whdeAucu[tv_1];
			Eigen::Matrix3d rotaMatr;
			double tempAngl = 2.0 * PI / dehwSurf.z[1] * (double)tv_0;
			rotaMatr << cos(tempAngl), -sin(tempAngl), 0.0,
				sin(tempAngl), cos(tempAngl), 0.0,
				0.0, 0.0, 1.0;
			Eigen::Vector3d tranVect = Eigen::Vector3d::Zero(3,1);
			tempCusu.RIGI_ROTR(rotaMatr, tranVect);
			//
			long ts = searCoun + dehwSurf.gridNumb[0][6] - 1 
				+ dehwSurf.gridNumb[0][6] - dehwSurf.circNumb 
				+ dehwSurf.gridNumb[1][5] * (dehwSurf.gridNumb[1][6] - 1) 
				+ tv_0 * dehwSurf.gridNumb[1][6] + tv_1;
			long tg_mast = contBody[ts][0];
			long tg_slav = contBody[ts][1];
			EFACE_SURFACE iterEfsu_0(&(multGrid[tg_mast]), &tempCusu);
			while(iterEfsu_0.INCREMENT() == 1){
				searCont[ts].mastSegm.emplace_back(iterEfsu_0.currNode);
			}
			EFACE_SURFACE iterEfsu_1(&(multGrid[tg_slav]), &tempCusu);
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
					auto iterNoco = multGrid[tg_mast].nodeCoor.find(node_tj);
					tempXico += sqrt(pow((iterNoco->second)[0], 2.0) 
						+ pow((iterNoco->second)[1], 2.0));
					tempEtac += (iterNoco->second)[2];
				}
				mastCoor[0].emplace_back(tempXico / 4.0);
				mastCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
			VECTOR2D slavCoor(2);
			for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
				double tempXico = 0.0;
				double tempEtac = 0.0;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = searCont[ts].slavSegm[ti][tj];
					auto iterNoco = multGrid[tg_slav].nodeCoor.find(node_tj);
					tempXico += sqrt(pow((iterNoco->second)[0], 2.0) 
						+ pow((iterNoco->second)[1], 2.0));
					tempEtac += (iterNoco->second)[2];
				}
				slavCoor[0].emplace_back(tempXico / 4.0);
				slavCoor[1].emplace_back(tempEtac / 4.0);
			}
			searCont[ts].CONTACT_SEARCH(slavCoor);
			searCont[ts].OUTPUT_INPO(ts);
		}
	}
	return 1;
}

long DEHW_1::CONT_INTE_NODD(){
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_NODD local mesh refinement");
	auto CART_CURV = [=](COOR tempCoor, double &tempXico, double &tempEtac){
		tempXico += tempCoor[1];
		tempEtac += sqrt(pow(tempCoor[0] + (dehwSurf.a_h2 + centErro), 2.0)
			+ pow(tempCoor[2], 2.0));
	};
	std::vector<std::vector<CURVEDS>> cusuTabl_0(4, std::vector<CURVEDS>(2));
	for(long tt = 0; tt < 4; tt ++){
		Eigen::Matrix3d wormRota_0;
		wormRota_0 << cos(analAngl[0]), -sin(analAngl[0]), 0.0, 
			sin(analAngl[0]), cos(analAngl[0]), 0.0, 
			0.0, 0.0, 1.0;
		Eigen::Matrix3d wormRota_1;
		wormRota_1 << 1.0, 0.0, 0.0, 
			0.0, 0.0, 1.0, 
			0.0, -1.0, 0.0;
		wormRota_1 = wormRota_1 * wormRota_0;
		Eigen::Vector3d wormTran;
		wormTran << - (dehwSurf.a_h2 + centErro), 0.0, 0.0;
		CURVEDS mastSurf = dehwSurf.wormTosu;
		mastSurf.RIGI_ROTR(wormRota_1, wormTran);
		double wheeAngl_1 = analAngl[1] + 2.0 * PI / dehwSurf.z[1] * (1.0 + tt);
		Eigen::Matrix3d wheeRota;
		wheeRota << cos(wheeAngl_1), -sin(wheeAngl_1), 0.0, 
			sin(wheeAngl_1), cos(wheeAngl_1), 0.0, 
			0.0, 0.0, 1.0;
		Eigen::Vector3d wheeTran;
		wheeTran << 0.0, 0.0, 0.0;
		CURVEDS slavSurf = dehwSurf.wheeTosu;
		slavSurf.RIGI_ROTR(wheeRota, wheeTran);
		cusuTabl_0[tt][0] = mastSurf;
		cusuTabl_0[tt][1] = slavSurf;
		//
		for(long tr = 0; tr < dehwSurf.locaLeve; tr ++){
			long buckFact = (1 << 
				(dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve - 1 
				- (dehwSurf.locaLeve - tr))
			);
			std::vector<long> buckNumb = {
				dehwSurf.gridNumb[0][4] * dehwSurf.gridNumb[0][6] * buckFact, 
				dehwSurf.gridNumb[0][3] * buckFact
			};
			CSEARCH searCont_0;
			bool tempFlag;
			searCont_0.ADAPTIVE_REFINE(
				&(multGrid[0]), &(multGrid[1]), tempFlag, &mastSurf, &slavSurf, 
				dehwSurf.globInho + dehwSurf.globHomo + tr, distCrit[tr], buckNumb, CART_CURV
			);
		}
	}
	//
	loadIncr.resize(2);
	SUBR_COLO_WORM(-1);
	SUBR_COLO_WHEE(-1);
	#pragma omp parallel for
	for(long tw = 0; tw < multGrid.size(); tw ++){
		if(tw < 1){
			SUBR_COLO_WORM(tw);
		}
		else{
			SUBR_COLO_WHEE(tw);
		}
	}
	#pragma omp parallel for
	for(long tw = 0; tw < multGrid.size(); tw ++){
		multGrid[tw].OUTPUT_ELEMENT(tw);
	}
	//****************************************************************************************
	long totaSear = 4;
	std::vector<std::vector<CURVEDS>> cusuTabl(totaSear, std::vector<CURVEDS>(2));
	searCont.resize(totaSear);
	penaFact_n.resize(totaSear);
	penaFact_f.resize(totaSear);
	fricCoef.resize(totaSear);
	contBody.resize(totaSear);
	long searCoun = 0;
	for(long tt = 0; tt < 4; tt ++){
		contBody[searCoun] = {0, 1};
		searCont[searCoun].mastGrid = &(multGrid[contBody[searCoun][0]]);
		searCont[searCoun].slavGrid = &(multGrid[contBody[searCoun][1]]);
		penaFact_n[searCoun] = (multGrid[contBody[searCoun][0]].mateElas 
			+ multGrid[contBody[searCoun][1]].mateElas) * 500.0;
		penaFact_f[searCoun] = (multGrid[contBody[searCoun][0]].mateElas 
			+ multGrid[contBody[searCoun][1]].mateElas) * 
			((coloSett == 1) ? 5.0 : 50.0);
		fricCoef[searCoun] = (coloSett == 1) ? 0.08 : 0.2;
		cusuTabl[searCoun][0] = cusuTabl_0[tt][0];
		cusuTabl[searCoun][1] = cusuTabl_0[tt][1];
		searCoun ++;
	}
	cusuTabl_0.clear();
	//****************************************************************************************
	OUTPUT_TIME("DEHW_1::CONT_INTE_NODD contact surface");
	long totaLeve = dehwSurf.globInho + dehwSurf.globHomo + dehwSurf.locaLeve;
	std::vector<long> buckNumb = {
		dehwSurf.gridNumb[0][4] * dehwSurf.gridNumb[0][6] * (1 << (totaLeve - 3)), 
		dehwSurf.gridNumb[0][3] * (1 << (totaLeve - 2))
	};
	#pragma omp parallel for
	for(long ts = 0; ts < searCoun; ts ++){
		EFACE_SURFACE iterEfsu_0(searCont[ts].mastGrid, &(cusuTabl[ts][0]));
		while(iterEfsu_0.INCREMENT() == 1){
			//NO NEED: if(multGrid[contBody[ts][0]].elemVect[iterEfsu_0.eid].level == totaLeve){
			searCont[ts].mastSegm.emplace_back(iterEfsu_0.currNode);
		}
		EFACE_SURFACE iterEfsu_1(searCont[ts].slavGrid, &(cusuTabl[ts][1]));
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
				auto iterNoco = multGrid[contBody[ts][0]].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			mastCoor[0].emplace_back(tempXico / 4.0);
			mastCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].BUCKET_SORT(mastCoor, buckNumb);
		VECTOR2D slavCoor(2);
		for(long ti = 0; ti < searCont[ts].slavSegm.size(); ti ++){
			double tempXico = 0.0;
			double tempEtac = 0.0;
			for(long tj = 0; tj < 4; tj ++){
				long node_tj = searCont[ts].slavSegm[ti][tj];
				auto iterNoco = multGrid[contBody[ts][1]].nodeCoor.find(node_tj);
				CART_CURV(iterNoco->second, tempXico, tempEtac);
			}
			slavCoor[0].emplace_back(tempXico / 4.0);
			slavCoor[1].emplace_back(tempEtac / 4.0);
		}
		searCont[ts].CONTACT_SEARCH(slavCoor, distCrit[dehwSurf.locaLeve - 1]);
		searCont[ts].OUTPUT_INPO(ts);
	}
	return 1;
}

long DEHW_1::SOLVE(long appsCont, long isnoDode){
	//****************************************************************************************
	dehwSurf.ESTABLISH();
	//
	if(appsCont == 0){
		muscSett = (1 << 1);
	}
	else{
		muscSett = whadCosp;
	}
	//
	distCrit = {55.0E-6, 35.0E-6, 15.0E-6};
	centErro = 0.0E-6;
	analAngl = {0.0, 0.0};//must be zero
	//
	if(isnoDode == 1){
		doleMcsc.assign(
			dehwSurf.gridNumb[0][6] + dehwSurf.gridNumb[1][5] * dehwSurf.gridNumb[1][6], 1
		);
		multGrid.resize(
			dehwSurf.gridNumb[0][6] + dehwSurf.gridNumb[1][5] * dehwSurf.gridNumb[1][6]
		);
		WORM_MESH(isnoDode);
		WHEE_MESH_DD();
		for(long tv = dehwSurf.gridNumb[0][6]; tv < multGrid.size(); tv ++){
			multGrid[tv].mateElas = 110.0E9;
		}
		CONT_INTE_DD();
	}
	else{
		doleMcsc.assign(2, 1);
		multGrid.resize(2);
		WORM_MESH(isnoDode);
		WHEE_MESH_NODD();
		multGrid[1].mateElas = 110.0E9;
		CONT_INTE_NODD();
	}
	// for(long tv = 0; tv < multGrid.size(); tv ++){
		// doleMcsc[tv] = multGrid[tv].mgpi.maxiLeve;
	// }
	//****************************************************************************************
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
