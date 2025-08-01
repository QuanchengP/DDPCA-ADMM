
#ifndef _MCONTACT_H
#define _MCONTACT_H

#include "CSEARCH.h"
#include "MGPIS.h"

//multibody contact
class MCONTACT{
public:
	/*********************************************************************************************/
	std::vector<MULTIGRID> multGrid;//multi-body contact
	std::vector<CSEARCH> searCont;//multi contact region
	VECTOR2L contBody;//contact body indices
	//fricCoef > 0: frictional contact
	//fricCoef = 0: normal contact
	//fricCoef < 0: perfect interface
	std::vector<double> fricCoef;//friction coefficient of each contact region
	//
	std::vector<double> penaFact_n;//normal penalty parameter of each contact region
	std::vector<double> penaFact_f;//tangential penalty parameter of each contact region
	long muscSett;//setting of coarse space correction
	std::vector<long> doleMcsc;//domain level of multiplicative coarse space correction
	/*******************************************Output********************************************/
	std::vector<Eigen::VectorXd> resuDisp;
	std::vector<std::vector<Eigen::VectorXd>> inteAuxi;
	std::vector<std::vector<Eigen::VectorXd>> inteLagr;
	/*********************************************************************************************/
	DIRE_SOLV mugrDiso[MAXI_DOMA_NUMB];//direct solver for subproblem of each subdomain
	DIRE_SOLV inteDiso[MAXI_INTE_NUMB][2];//direct solver for interface subproblem
	DIRE_SOLV inteDiso_pena[MAXI_INTE_NUMB][2];//penalty * inteDiso
	//node to contact number
	std::vector<std::vector<std::map<long,long>>> nodeCont;
	//penalty matrix for integral point
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> pemaInpo;
	std::vector<std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>>> pemaInpo_r;
	SPARMATS systMass;//system mass matrix
	SPARMATS systTran;//system transfer matrix
	SPARMATS systTran_pena;
	SPARMATS inteMass;//interface mass matrix
	SPARMATS inteMass_pena;
	SPARMATS inteMass_pere;
	SPARMATS inpoLagr;//Lagrange at integral point
	SPARMATS inpoDisp;//displacement at integral point
	SPARMATS inteInpo;//from integral point to iterface
	std::vector<Eigen::VectorXd> inpoNgap;//initial normal gap
	//output contact pressure/traction
	long OUTPUT_PRTR(const Eigen::VectorXd &inpoGamm, const Eigen::VectorXi &fricStat, long ts);
	//output auxiliary variable/lagrange multiplier
	long OUTPUT_AULA();
	//node coordinate system is changed to cylindrical coordinate
	long CONT_ROTA(const MULTIGRID &inpuMugr, const std::vector<long> &inpuCono, 
		Eigen::SparseMatrix<double,Eigen::RowMajor> &resuRota
	);
	long ESTABLISH();
	long CONTACT_ANALYSIS();
	//output monitor information at each step
	long MONITOR(const long tc, std::ofstream &tempOfst, VECTOR2D &moniReco, 
		const std::vector<Eigen::VectorXd> &resuDisp_0, 
		const std::vector<std::vector<Eigen::VectorXd>> &inteAuxi_0, 
		const std::vector<std::vector<Eigen::VectorXd>> &inteLagr_0
	);
	/*********************************************************************************************/
	//
	Eigen::SparseMatrix<double,Eigen::RowMajor> globCoup;
	SPARMATS globTran;
	SPARMATS globTran_pena;
	SPARMATS globTran_D;
	std::vector<long> baseReco;
	DIRE_SOLV coarSolv_D;
	COGR_SOLV coarSolv_C;
	MGPIS mgpi;
	long DOUBLE_M(VECTOR3L coarNode);//
	long MULTISCALE();
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> accuProl;//accumulated prolong
	/*********************************************************************************************/
	//!!!!!!!!!!only for perfect interface!!!!!!!!!!
	Eigen::SparseMatrix<double,Eigen::RowMajor> globCoup_1;
	SPARMATS globTran_1;
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> globTran_D_1;
	Eigen::VectorXd globForc_1;
	DIRE_SOLV coarSolv_D_1;
	COGR_SOLV coarSolv_C_1;
	MGPIS mgpi_1;
	long DOUBLE_M_1();
	long MULTISCALE_1();
	/*********************************************************************************************/
	long APPS_MPL();
	long APPS();//automatic penalty parameter selection
	/*********************************************************************************************/
	//lagrange multiplier method: the levels of all bodies are the same
	long LAGRANGE(long precType);//can not deal with cross corner problem
};

long MCONTACT::OUTPUT_PRTR(const Eigen::VectorXd &inpoGamm, 
	const Eigen::VectorXi &fricStat, long ts){
	//
	std::stringstream tempStre;
	tempStre << "resuCont_" << ts << ".txt";
	std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
	tempStre.str("");
	tempStre.clear();
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
		if(fricCoef[ts] == 0.0){
			tempOfst << std::setw(30) << inpoGamm(ti) << std::endl;
		}
		else{
			Eigen::Vector3d tempTang_0 = searCont[ts].intePoin[ti].basiVect[1].transpose();
			Eigen::Vector3d tempTang_1 = searCont[ts].intePoin[ti].basiVect[2].transpose();
			Eigen::Vector3d tangTrac = inpoGamm(3 * ti + 1) * tempTang_0 
				+ inpoGamm(3 * ti + 2) * tempTang_1;
			tempOfst << std::setw(30) << inpoGamm(3 * ti + 0) 
				<< std::setw(30) << tangTrac(0) << std::setw(30) << tangTrac(1) 
				<< std::setw(30) << tangTrac(2) 
				<< std::setw(10) << fricStat(3 * ti + 1) << std::endl;
		}
	}
	tempOfst.close();
	return 1;
}

long MCONTACT::OUTPUT_AULA(){
	#pragma omp parallel for
	for(long ts = 0; ts < inteAuxi.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			//
			std::stringstream tempStre;
			tempStre << "resuAula_" << ts << "_" << tv << ".txt";
			std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
			tempStre.str("");
			tempStre.clear();
			tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
			for(long ti = 0; ti < nodeCont[ts][tv].size(); ti ++){
				if(fricCoef[ts] == 0.0){
					tempOfst << std::setw(30) << inteAuxi[ts][tv](ti);
					tempOfst << std::setw(30) << inteLagr[ts][tv](ti);
				}
				else{
					tempOfst << std::setw(30) << inteAuxi[ts][tv](3 * ti + 0) 
						<< std::setw(30) << inteAuxi[ts][tv](3 * ti + 1) 
						<< std::setw(30) << inteAuxi[ts][tv](3 * ti + 2);
					tempOfst << std::setw(30) << inteLagr[ts][tv](3 * ti + 0) 
						<< std::setw(30) << inteLagr[ts][tv](3 * ti + 1) 
						<< std::setw(30) << inteLagr[ts][tv](3 * ti + 2);
				}
				tempOfst << std::endl;
			}
			tempOfst.close();
		}
	}
	return 1;
}

long MCONTACT::CONT_ROTA(const MULTIGRID &inpuMugr, const std::vector<long> &inpuCono, 
	Eigen::SparseMatrix<double,Eigen::RowMajor> &resuRota){
	std::vector<Eigen::Triplet<double>> tempList;
	tempList.reserve(inpuCono.size() * 3 * 3);
	for(long ti = 0; ti < inpuCono.size(); ti ++){
		auto iterNoro = inpuMugr.nodeRota.find(inpuCono[ti]);
		if(iterNoro == inpuMugr.nodeRota.end()){
			for(long tk = 0; tk < 3; tk ++){
				tempList.emplace_back(3 * ti + tk, 3 * ti + tk, 1.0);
			}
		}
		else{
			for(long tj = 0; tj < 3; tj ++){
				for(long tk = 0; tk < 3; tk ++){
					tempList.emplace_back(3 * ti + tj, 3 * ti + tk, (iterNoro->second)(tj,tk));
				}
			}
		}
	}
	resuRota.resize(3 * inpuCono.size(), 3 * inpuCono.size());
	resuRota.setFromTriplets(tempList.begin(), tempList.end());
	return 1;
}

long MCONTACT::ESTABLISH(){
	if(multGrid.size() > MAXI_DOMA_NUMB){
		OUTPUT_TIME("MCONTACT::ESTABLISH ERROR 1");
	}
	if(doleMcsc.size() == 0){
		doleMcsc.assign(multGrid.size(), 0);
	}
	//***************************************************************************************
	nodeCont.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		nodeCont[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			// std::stringstream tempStre;
			// tempStre << "resuNoco_" << ts << "_" << tv << ".txt";
			// std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
			// tempStre.str("");
			// tempStre.clear();
			for(const auto &iterInpo : searCont[ts].intePoin){
				for(long tk = 0; tk < 4; tk ++){
					long tempNode = iterInpo.node[tv][tk];
					auto iterNoco = nodeCont[ts][tv].find(tempNode);
					if(iterNoco == nodeCont[ts][tv].end()){
						long tempSize = nodeCont[ts][tv].size();
						nodeCont[ts][tv].emplace(tempNode, tempSize);
						// tempOfst << std::setw(10) << tempNode << std::endl;
					}
				}
			}
			// tempOfst.close();
		}
	}
	pemaInpo.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		std::vector<Eigen::Triplet<double>> tempList;
		for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
			if(fricCoef[ts] == 0.0){
				tempList.emplace_back(ti, ti, penaFact_n[ts]);
			}
			else{
				tempList.emplace_back(3 * ti + 0, 3 * ti + 0, penaFact_n[ts]);
				tempList.emplace_back(3 * ti + 1, 3 * ti + 1, penaFact_f[ts]);
				tempList.emplace_back(3 * ti + 2, 3 * ti + 2, penaFact_f[ts]);
			}
		}
		if(fricCoef[ts] == 0.0){
			pemaInpo[ts].resize(
				searCont[ts].intePoin.size(), searCont[ts].intePoin.size()
			);
		}
		else{
			pemaInpo[ts].resize(
				3 * searCont[ts].intePoin.size(), 3 * searCont[ts].intePoin.size()
			);
		}
		pemaInpo[ts].setFromTriplets(tempList.begin(), tempList.end());
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::ESTABLISH systMass");
	systMass.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		systMass[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> normList;
				normList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,3> tempNorm = iterInpo.basiVect[0].transpose();
					std::vector<double> M_e = iterInpo.shapFunc[tv];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 0.0, 
						0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 
						0.0, 0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3];
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig * 
						N_e.transpose() * tempNorm.transpose() * penaFact_n[ts] * tempNorm * N_e;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								long node_tm = iterInpo.node[tv][tm];
								for(long tn = 0; tn < 3; tn ++){
									long free_tn = 3 * node_tm + tn;
									normList.emplace_back(free_tk, free_tn, 
										matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				systMass[ts][tv].resize(tempNumb, tempNumb);
				systMass[ts][tv].setFromTriplets(normList.begin(), normList.end());
			}
			else{
				std::vector<Eigen::Triplet<double>> fricList;
				fricList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					std::vector<double> M_e = iterInpo.shapFunc[tv];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 0.0, 
						0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 
						0.0, 0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3];
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					//tempTang.transpose() * tempTang = I
					Eigen::MatrixXd matr_1 = iterInpo.quadWeig * 
						N_e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								long node_tm = iterInpo.node[tv][tm];
								for(long tn = 0; tn < 3; tn ++){
									long free_tn = 3 * node_tm + tn;
									fricList.emplace_back(free_tk, free_tn, 
										matr_1(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				systMass[ts][tv].resize(tempNumb, tempNumb);
				systMass[ts][tv].setFromTriplets(fricList.begin(), fricList.end());
			}
		}
	}
	//
	OUTPUT_TIME("MCONTACT::ESTABLISH systTran");
	systTran.resize(searCont.size());
	systTran_pena.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		systTran[ts].resize(2);
		systTran_pena[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> normList;
				normList.reserve(searCont[ts].intePoin.size() * 12 * 4);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,3> tempNorm = iterInpo.basiVect[0].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig 
						* N_e.transpose() * tempNorm.transpose() * M_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_0 = tempRota.transpose() * matr_0;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iterNoco = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								normList.emplace_back(free_tk, iterNoco->second, 
									matr_0(3 * tj + tk, tm)
								);
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				systTran[ts][tv].resize(tempNumb, nodeCont[ts][tv].size());
				systTran[ts][tv].setFromTriplets(normList.begin(), normList.end());
				systTran_pena[ts][tv] = penaFact_n[ts] * systTran[ts][tv];
			}
			else{
				std::vector<Eigen::Triplet<double>> fricList, fricList_1;
				fricList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				fricList_1.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					//tempTang.transpose() * tempTang = I
					Eigen::MatrixXd matr_1 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempTang * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_1 = tempRota.transpose() * matr_1;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iterNoco = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								for(long tn = 0; tn < 3; tn ++){
									long col_tn = 3 * iterNoco->second + tn;
									fricList.emplace_back(free_tk, col_tn, 
										matr_1(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					matr_1 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					matr_1 = tempRota.transpose() * matr_1;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iterNoco = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								for(long tn = 0; tn < 3; tn ++){
									long col_tn = 3 * iterNoco->second + tn;
									fricList_1.emplace_back(free_tk, col_tn, 
										matr_1(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				systTran[ts][tv].resize(tempNumb, 3 * nodeCont[ts][tv].size());
				systTran[ts][tv].setFromTriplets(fricList.begin(), fricList.end());
				systTran_pena[ts][tv].resize(tempNumb, 3 * nodeCont[ts][tv].size());
				systTran_pena[ts][tv].setFromTriplets(fricList_1.begin(), fricList_1.end());
			}
		}
	}
	//
	OUTPUT_TIME("MCONTACT::ESTABLISH inteMass");
	inteMass.resize(searCont.size());
	inteMass_pere.resize(searCont.size());
	inteMass_pena.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		inteMass[ts].resize(2);
		inteMass_pere[ts].resize(2);
		inteMass_pena[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> noinList;
				noinList.reserve(searCont[ts].intePoin.size() * 4 * 4);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig * M_e.transpose() * M_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][tv].find(iterInpo.node[tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][tv].find(iterInpo.node[tv][tk]);
							long col_tk = iter_Ck->second;
							noinList.emplace_back(row_tj, col_tk, matr_0(tj, tk));
						}
					}
				}
				inteMass[ts][tv].resize(nodeCont[ts][tv].size(), nodeCont[ts][tv].size());
				inteMass[ts][tv].setFromTriplets(noinList.begin(), noinList.end());
				inteMass_pere[ts][tv] = inteMass[ts][tv] / penaFact_n[ts];
				inteMass_pena[ts][tv] = inteMass[ts][tv] * penaFact_n[ts];
			}
			else{
				std::vector<Eigen::Triplet<double>> frinList, frinList_1, frinList_2;
				frinList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				frinList_1.reserve(searCont[ts].intePoin.size() * 12 * 12);
				frinList_2.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					//tempTang.transpose() * tempTang = I
					Eigen::MatrixXd matr_1 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][tv].find(iterInpo.node[tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][tv].find(iterInpo.node[tv][tk]);
							long col_tk = iter_Ck->second;
							for(long tl = 0; tl < 3; tl ++){
								for(long tm = 0; tm < 3; tm ++){
									frinList.emplace_back(3 * row_tj + tl, 3 * col_tk + tm, 
										matr_1(3 * tj + tl, 3 * tk + tm)
									);
								}
							}
						}
					}
					Eigen::Matrix3d tempPena;
					tempPena << 1.0 / penaFact_n[ts], 0.0, 0.0, 
						0.0, 1.0 / penaFact_f[ts], 0.0, 
						0.0, 0.0, 1.0 / penaFact_f[ts];
					matr_1 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][tv].find(iterInpo.node[tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][tv].find(iterInpo.node[tv][tk]);
							long col_tk = iter_Ck->second;
							for(long tl = 0; tl < 3; tl ++){
								for(long tm = 0; tm < 3; tm ++){
									frinList_1.emplace_back(3 * row_tj + tl, 3 * col_tk + tm, 
										matr_1(3 * tj + tl, 3 * tk + tm)
									);
								}
							}
						}
					}
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					matr_1 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][tv].find(iterInpo.node[tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][tv].find(iterInpo.node[tv][tk]);
							long col_tk = iter_Ck->second;
							for(long tl = 0; tl < 3; tl ++){
								for(long tm = 0; tm < 3; tm ++){
									frinList_2.emplace_back(3 * row_tj + tl, 3 * col_tk + tm, 
										matr_1(3 * tj + tl, 3 * tk + tm)
									);
								}
							}
						}
					}
				}
				inteMass[ts][tv].resize(3 * nodeCont[ts][tv].size(), 3 * nodeCont[ts][tv].size());
				inteMass[ts][tv].setFromTriplets(frinList.begin(), frinList.end());
				inteMass_pere[ts][tv].resize(
					3 * nodeCont[ts][tv].size(), 3 * nodeCont[ts][tv].size()
				);
				inteMass_pere[ts][tv].setFromTriplets(frinList_1.begin(), frinList_1.end());
				inteMass_pena[ts][tv].resize(
					3 * nodeCont[ts][tv].size(), 3 * nodeCont[ts][tv].size()
				);
				inteMass_pena[ts][tv].setFromTriplets(frinList_2.begin(), frinList_2.end());
			}
		}
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::ESTABLISH inpoLagr");
	inpoLagr.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		inpoLagr[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> normList;
				normList.reserve(searCont[ts].intePoin.size() * 4);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						auto iter_Cj = nodeCont[ts][tv].find(node_tj);
						normList.emplace_back(ti, iter_Cj->second, M_e(tj));
					}
				}
				inpoLagr[ts][tv].resize(searCont[ts].intePoin.size(), nodeCont[ts][tv].size());
				inpoLagr[ts][tv].setFromTriplets(normList.begin(), normList.end());
			}
			else{
				std::vector<Eigen::Triplet<double>> fricList;
				fricList.reserve(searCont[ts].intePoin.size() * 12);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = searCont[ts].intePoin[ti].basiVect[0].transpose();
					tempTang.block(1,0,1,3) = searCont[ts].intePoin[ti].basiVect[1].transpose();
					tempTang.block(2,0,1,3) = searCont[ts].intePoin[ti].basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_1 = tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						auto iter_Cj = nodeCont[ts][tv].find(node_tj);
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 3; tm ++){
								fricList.emplace_back(3 * ti + tm, 3 * iter_Cj->second + tk, 
									matr_1(tm, 3 * tj + tk)
								);
							}
						}
					}
				}
				inpoLagr[ts][tv].resize(3 * searCont[ts].intePoin.size(), 
					3 * nodeCont[ts][tv].size()
				);
				inpoLagr[ts][tv].setFromTriplets(fricList.begin(), fricList.end());
			}
		}
	}
	//
	OUTPUT_TIME("MCONTACT::ESTABLISH inpoDisp");
	inpoDisp.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		inpoDisp[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> normList;
				normList.reserve(searCont[ts].intePoin.size() * 1 * 12);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,1,3> tempNorm = 
						searCont[ts].intePoin[ti].basiVect[0].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_0 = tempNorm * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], 
						searCont[ts].intePoin[ti].node[tv], tempRota
					);
					matr_0 = matr_0 * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							normList.emplace_back(ti, 3 * node_tj + tk, 
								matr_0(0, 3 * tj + tk)
							);
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				inpoDisp[ts][tv].resize(searCont[ts].intePoin.size(), tempNumb);
				inpoDisp[ts][tv].setFromTriplets(normList.begin(), normList.end());
			}
			else{
				std::vector<Eigen::Triplet<double>> fricList;
				fricList.reserve(searCont[ts].intePoin.size() * 3 * 12);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = searCont[ts].intePoin[ti].basiVect[0].transpose();
					tempTang.block(1,0,1,3) = searCont[ts].intePoin[ti].basiVect[1].transpose();
					tempTang.block(2,0,1,3) = searCont[ts].intePoin[ti].basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_1 = tempTang * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], 
						searCont[ts].intePoin[ti].node[tv], tempRota
					);
					matr_1 = matr_1 * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 3; tm ++){
								fricList.emplace_back(3 * ti + tm, 3 * node_tj + tk, 
									matr_1(tm, 3 * tj + tk)
								);
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				inpoDisp[ts][tv].resize(3 * searCont[ts].intePoin.size(), tempNumb);
				inpoDisp[ts][tv].setFromTriplets(fricList.begin(), fricList.end());
			}
		}
	}
	//
	OUTPUT_TIME("MCONTACT::ESTABLISH inteInpo");
	VECT_RESI(inteInpo, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> normList;
				normList.reserve(searCont[ts].intePoin.size() * 4 * 1);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					Eigen::MatrixXd matr_0 = searCont[ts].intePoin[ti].quadWeig 
						* M_e.transpose();
					if(tv == 0){
						matr_0 = - matr_0;
					}
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						auto iter_Cj = nodeCont[ts][tv].find(node_tj);
						normList.emplace_back(iter_Cj->second, ti, matr_0(tj, 0));
					}
				}
				inteInpo[ts][tv].resize(nodeCont[ts][tv].size(), searCont[ts].intePoin.size());
				inteInpo[ts][tv].setFromTriplets(normList.begin(), normList.end());
			}
			else{
				std::vector<Eigen::Triplet<double>> fricList;
				fricList.reserve(searCont[ts].intePoin.size() * 12 * 3);
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = searCont[ts].intePoin[ti].basiVect[0].transpose();
					tempTang.block(1,0,1,3) = searCont[ts].intePoin[ti].basiVect[1].transpose();
					tempTang.block(2,0,1,3) = searCont[ts].intePoin[ti].basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << searCont[ts].intePoin[ti].shapFunc[tv][0], 
						searCont[ts].intePoin[ti].shapFunc[tv][1], 
						searCont[ts].intePoin[ti].shapFunc[tv][2], 
						searCont[ts].intePoin[ti].shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_1 = searCont[ts].intePoin[ti].quadWeig 
						* N_e.transpose() * tempTang.transpose();
					if(tv == 0){
						matr_1 = - matr_1;
					}
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = searCont[ts].intePoin[ti].node[tv][tj];
						auto iter_Cj = nodeCont[ts][tv].find(node_tj);
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 3; tm ++){
								fricList.emplace_back(3 * iter_Cj->second + tk, 3 * ti + tm, 
									matr_1(3 * tj + tk, tm)
								);
							}
						}
					}
				}
				inteInpo[ts][tv].resize(3 * nodeCont[ts][tv].size(), 
					3 * searCont[ts].intePoin.size()
				);
				inteInpo[ts][tv].setFromTriplets(fricList.begin(), fricList.end());
			}
		}
	}
	//
	OUTPUT_TIME("MCONTACT::ESTABLISH inpoNgap");
	inpoNgap.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		if(fricCoef[ts] == 0.0){
			inpoNgap[ts] = Eigen::VectorXd::Zero(searCont[ts].intePoin.size());
			for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
				inpoNgap[ts](ti) = searCont[ts].intePoin[ti].initNgap;
			}
		}
		else{
			inpoNgap[ts] = Eigen::VectorXd::Zero(3 * searCont[ts].intePoin.size());
			for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
				inpoNgap[ts](3 * ti + 0) = searCont[ts].intePoin[ti].initNgap;
			}
		}
	}
	//
	VECT_RESI(pemaInpo_r,searCont.size(),2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		pemaInpo_r[ts][0] = pemaInpo[ts] * inpoDisp[ts][0];
		pemaInpo_r[ts][1] = pemaInpo[ts] * inpoDisp[ts][1];
	}
	//***************************************************************************************
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		multGrid[tv].TRANSFER();
		multGrid[tv].STIF_MATR();
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(long ti = 0; ti < 2; ti ++){
				if(tv != contBody[ts][ti]){
					continue;
				}
				multGrid[tv].origStif[multGrid[tv].mgpi.maxiLeve + 1] += systMass[ts][ti];
			}
		}
		multGrid[tv].CONSTRAINT(1);
	}
	OUTPUT_TIME("Factorization of consStif");
	#pragma omp parallel for// if(multGrid.size() <= omp_get_max_threads() / 2)
	for(long tv = 0; tv < multGrid.size(); tv ++){
		const Eigen::SparseMatrix<double,Eigen::RowMajor> &tempCost = 
			multGrid[tv].mgpi.consStif[multGrid[tv].mgpi.maxiLeve];
		if(tempCost.rows() < DIRE_MAXI_SUBD){
			(mugrDiso[tv]).compute(tempCost);
			std::cout << " " << tv << " ";
		}
	}
	std::cout << std::endl;
	OUTPUT_TIME("Factorization of inteMass");
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			if(inteMass[ts][tv].rows() < DIRE_MAXI){
				(inteDiso[ts][tv]).compute(inteMass[ts][tv]);
				(inteDiso_pena[ts][tv]).compute(inteMass_pena[ts][tv]);
				std::cout << " " << ts << "-" << tv << " ";
			}
		}
	}
	std::cout << std::endl;
	//***************************************************************************************
	long baseNumb = 0;
	baseReco.resize(multGrid.size() + 1);
	for(long tv = 0; tv < multGrid.size(); tv ++){
		const auto &tempStif = multGrid[tv].mgpi.consStif[doleMcsc[tv]];
		baseReco[tv] = baseNumb;
		baseNumb += tempStif.rows();
	}
	baseReco[multGrid.size()] = baseNumb;
	if((muscSett >> 0) % 2 == 1){
		MULTISCALE();
	}
	if((muscSett >> 1) % 2 == 1){
		MULTISCALE_1();
	}
	accuProl.resize(multGrid.size());
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		accuProl[tv] = multGrid[tv].consOper[doleMcsc[tv]].transpose();
		for(long ti = doleMcsc[tv]; ti < multGrid[tv].mgpi.maxiLeve; ti ++){
			accuProl[tv] = multGrid[tv].prolOper[ti] * accuProl[tv];
		}
		accuProl[tv] = multGrid[tv].consOper[multGrid[tv].mgpi.maxiLeve] * accuProl[tv];
	}
	//***************************************************************************************
	//initialization
	resuDisp.resize(multGrid.size());
	inteAuxi.resize(searCont.size());
	inteLagr.resize(searCont.size());
	for(long tv = 0; tv < multGrid.size(); tv ++){//for MONITOR
		resuDisp[tv] = Eigen::VectorXd::Zero(3 * multGrid[tv].nodeCoor.size());
	}
	for(long ts = 0; ts < searCont.size(); ts ++){//necessary
		inteAuxi[ts].resize(2);
		inteLagr[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				inteAuxi[ts][tv] = Eigen::VectorXd::Zero(nodeCont[ts][tv].size());
				inteLagr[ts][tv] = Eigen::VectorXd::Zero(nodeCont[ts][tv].size());
			}
			else{
				inteAuxi[ts][tv] = Eigen::VectorXd::Zero(3 * nodeCont[ts][tv].size());
				inteLagr[ts][tv] = Eigen::VectorXd::Zero(3 * nodeCont[ts][tv].size());
			}
		}
	}
	return 1;
}

long MCONTACT::MULTISCALE(){
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE ficoCotr");
	SPARMATS ficoCotr;//contact node transfer from coarset level to finest level
	ficoCotr.resize(searCont.size());
	VECTOR3L coarNode;
	VECT_RESI(coarNode, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		ficoCotr[ts].resize(2);
		for(long tv = 0; tv < 2; tv ++){
			//
			const auto &tempGrid =  multGrid[contBody[ts][tv]];
			Eigen::SparseMatrix<double,Eigen::RowMajor> tempFico = tempGrid.scalEarl;
			for(long ti = tempGrid.mgpi.maxiLeve; ti >= doleMcsc[contBody[ts][tv]]; ti --){
				tempFico = tempFico * tempGrid.scalProl[ti];
			}
			//!!!!!no need: tempFico * consOper[doleMcsc[contBody[ts][tv]]].transpose()
			//from total to contact
			Eigen::SparseMatrix<double,Eigen::RowMajor> contTota;
			std::vector<Eigen::Triplet<double>> tempList;
			for(const auto &iterNoco : nodeCont[ts][tv]){
				tempList.emplace_back(iterNoco.second, iterNoco.first, 1.0);
			}
			contTota.resize(nodeCont[ts][tv].size(), tempGrid.nodeCoor.size());
			contTota.setFromTriplets(tempList.begin(), tempList.end());
			tempFico = contTota * tempFico;
			//
			Eigen::SparseMatrix<double,Eigen::ColMajor> coluMatr = tempFico;
			long tempColu = 0;
			if(fricCoef[ts] == 0.0){
				std::vector<Eigen::Triplet<double>> tempList_f;
				for(long ti = 0; ti < coluMatr.cols(); ti ++){
					long coluTota = 0;
					for(CSPA_INNE iterStif(coluMatr, ti); iterStif; ++ iterStif){
						tempList_f.emplace_back(iterStif.row(), tempColu, iterStif.value());
						coluTota ++;
					}
					if(coluTota > 0){
						coarNode[ts][tv].emplace_back(ti);
						tempColu ++;
					}
				}
				ficoCotr[ts][tv].resize(nodeCont[ts][tv].size(), tempColu);
				ficoCotr[ts][tv].setFromTriplets(tempList_f.begin(), tempList_f.end());
			}
			else{
				std::vector<Eigen::Triplet<double>> tempList_f;
				for(long ti = 0; ti < coluMatr.cols(); ti ++){
					long coluTota = 0;
					for(CSPA_INNE iterStif(coluMatr, ti); iterStif; ++ iterStif){
						for(long tj = 0; tj < 3; tj ++){
							tempList_f.emplace_back(3 * iterStif.row() + tj, 3 * tempColu + tj, 
								iterStif.value()
							);
						}
						coluTota ++;
					}
					if(coluTota > 0){
						coarNode[ts][tv].emplace_back(ti);
						tempColu ++;
					}
				}
				ficoCotr[ts][tv].resize(3 * nodeCont[ts][tv].size(), 3 * tempColu);
				ficoCotr[ts][tv].setFromTriplets(tempList_f.begin(), tempList_f.end());
			}
		}
	}
	long contNumb = 0;
	std::vector<long> contReco = std::vector<long>(searCont.size(), 0);
	for(long ts = 0; ts < searCont.size(); ts ++){
		long acti_tv = 0;
		contReco[ts] = contNumb;
		contNumb += ficoCotr[ts][acti_tv].cols();
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE globCoup");
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> domaCoup(multGrid.size());
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		std::vector<Eigen::Triplet<double>> globList;
		//subtraction of the same thing
		//no need of special treatment for constrained DOF: F_F - K_C * U_C 
		//(subscript: F: free, C: constraint)
		const auto &tempStif = multGrid[tv].mgpi.consStif[doleMcsc[tv]];
		for(long ti = 0; ti < tempStif.rows(); ti ++){
			for(RSPA_INNE iterStif(tempStif, ti); iterStif; ++ iterStif){
				globList.emplace_back(baseReco[tv] + ti, 
					baseReco[tv] + iterStif.col(), iterStif.value()
				);
			}
		}
		domaCoup[tv].resize(
			baseReco[multGrid.size()] + contNumb, baseReco[multGrid.size()] + contNumb
		);
		domaCoup[tv].setFromTriplets(globList.begin(), globList.end());
		std::cout << " " << tv << " ";
	}
	globCoup = domaCoup[0];
	for(long tv = 1; tv < multGrid.size(); tv ++){
		globCoup += domaCoup[tv];
	}
	//
	SPARMATS inteCoup;
	VECT_RESI(inteCoup, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		long acti_tv = 0;
		//
		for(long tv = 0; tv < 2; tv ++){
			std::vector<Eigen::Triplet<double>> globList;
			if(fricCoef[ts] == 0.0){
				//
				std::vector<Eigen::Triplet<double>> dispList, unbaList;
				dispList.reserve(searCont[ts].intePoin.size() * 12 * 4);
				unbaList.reserve(searCont[ts].intePoin.size() * 4 * 4);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,3> tempNorm = iterInpo.basiVect[0].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig * penaFact_n[ts] 
						* N_e.transpose() * tempNorm.transpose() * M_0e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_0 = tempRota.transpose() * matr_0;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iterNoco = 
									nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tm]);
								dispList.emplace_back(free_tk, iterNoco->second, 
									matr_0(3 * tj + tk, tm)
								);
							}
						}
					}
					matr_0 = iterInpo.quadWeig * penaFact_n[ts] * M_0e.transpose() * M_0e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tk]);
							long col_tk = iter_Ck->second;
							unbaList.emplace_back(row_tj, col_tk, matr_0(tj, tk));
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> dispUnba;
				dispUnba.resize(tempNumb, nodeCont[ts][acti_tv].size());
				dispUnba.setFromTriplets(dispList.begin(), dispList.end());
				Eigen::SparseMatrix<double,Eigen::RowMajor> unbaMatr;
				unbaMatr.resize(nodeCont[ts][acti_tv].size(), nodeCont[ts][acti_tv].size());
				unbaMatr.setFromTriplets(unbaList.begin(), unbaList.end());
				//
				const auto &tempGrid = multGrid[contBody[ts][tv]];
				dispUnba = 
					tempGrid.earlTran.transpose() * dispUnba * ficoCotr[ts][acti_tv];
				for(long ti = tempGrid.mgpi.maxiLeve; ti >= doleMcsc[contBody[ts][tv]]; ti --){
					dispUnba = tempGrid.prolOper[ti].transpose() * dispUnba;
				}
				dispUnba = tempGrid.consOper[doleMcsc[contBody[ts][tv]]] * dispUnba;
				for(long ti = 0; ti < dispUnba.rows(); ti ++){
					for(RSPA_INNE iterStif(dispUnba, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							- iterStif.value()
						);
						globList.emplace_back(
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							baseReco[contBody[ts][tv]] + ti, 
							- iterStif.value()
						);
					}
				}
				//
				unbaMatr = ficoCotr[ts][acti_tv].transpose() * unbaMatr * ficoCotr[ts][acti_tv];
				for(long ti = 0; ti < unbaMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(unbaMatr, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							iterStif.value()
						);
					}
				}
			}
			else{
				//
				std::vector<Eigen::Triplet<double>> dispList, unbaList;
				dispList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				unbaList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::Matrix<double,3,12> N_0e;
					N_0e << 
					M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 0.0, 
					0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 
					0.0, 0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3);
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig 
						* N_e.transpose() * tempTang.transpose() * tempPena * tempTang * N_0e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_0 = tempRota.transpose() * matr_0;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iterNoco = 
									nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tm]);
								for(long tn = 0; tn < 3; tn ++){
									dispList.emplace_back(free_tk, 
										3 * iterNoco->second + tn, 
										matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
					matr_0 = iterInpo.quadWeig 
						* N_0e.transpose() * tempTang.transpose() * tempPena * tempTang * N_0e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk  < 3; tk ++){
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = 
									nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tm]);
								long col_tm = iter_Cm->second;
								for(long tn = 0; tn < 3; tn ++){
									unbaList.emplace_back(3 * row_tj + tk, 3 * col_tm + tn, 
										matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> dispUnba;
				dispUnba.resize(tempNumb, 3 * nodeCont[ts][acti_tv].size());
				dispUnba.setFromTriplets(dispList.begin(), dispList.end());
				Eigen::SparseMatrix<double,Eigen::RowMajor> unbaMatr;
				unbaMatr.resize(3 * nodeCont[ts][acti_tv].size(), 
					3 * nodeCont[ts][acti_tv].size()
				);
				unbaMatr.setFromTriplets(unbaList.begin(), unbaList.end());
				//
				const auto &tempGrid = multGrid[contBody[ts][tv]];
				dispUnba = 
					tempGrid.earlTran.transpose() * dispUnba * ficoCotr[ts][acti_tv];
				for(long ti = tempGrid.mgpi.maxiLeve; ti >= doleMcsc[contBody[ts][tv]]; ti --){
					dispUnba = tempGrid.prolOper[ti].transpose() * dispUnba;
				}
				dispUnba = tempGrid.consOper[doleMcsc[contBody[ts][tv]]] * dispUnba;
				for(long ti = 0; ti < dispUnba.rows(); ti ++){
					for(RSPA_INNE iterStif(dispUnba, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							- iterStif.value()
						);
						globList.emplace_back(
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							baseReco[contBody[ts][tv]] + ti, 
							- iterStif.value()
						);
					}
				}
				//
				unbaMatr = ficoCotr[ts][acti_tv].transpose() * unbaMatr * ficoCotr[ts][acti_tv];
				for(long ti = 0; ti < unbaMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(unbaMatr, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							baseReco[multGrid.size()] + contReco[ts] + iterStif.col(), 
							iterStif.value()
						);
					}
				}
			}
			inteCoup[ts][tv].resize(
				baseReco[multGrid.size()] + contNumb, baseReco[multGrid.size()] + contNumb
			);
			inteCoup[ts][tv].setFromTriplets(globList.begin(), globList.end());
			std::cout << " " << ts << "-" << tv << " ";
		}
	}
	std::cout << std::endl;
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			globCoup += inteCoup[ts][tv];
		}
	}
	std::cout << "Dimension of globCoup: " << globCoup.rows();
	OUTPUT_TIME("");
	if(globCoup.rows() < DIRE_MAXI){
		coarSolv_D.compute(globCoup);
	}
	else if(globCoup.rows() < COGR_MAXI){
		coarSolv_C.compute(globCoup);
	}
	else{
		DOUBLE_M(coarNode);
		mgpi.ESTABLISH();
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE globTran");
	globTran.resize(searCont.size());
	globTran_pena.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		globTran[ts].resize(2);
		globTran_pena[ts].resize(2);
		long acti_tv = 0;
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				//
				std::vector<Eigen::Triplet<double>> tranList;
				tranList.reserve(searCont[ts].intePoin.size() * 4 * 4);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig * M_0e.transpose() * M_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							auto iter_Ck = nodeCont[ts][tv].find(iterInpo.node[tv][tk]);
							long col_tk = iter_Ck->second;
							tranList.emplace_back(row_tj, col_tk, matr_0(tj, tk));
						}
					}
				}
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr;
				tempMatr.resize(nodeCont[ts][acti_tv].size(), nodeCont[ts][tv].size());
				tempMatr.setFromTriplets(tranList.begin(), tranList.end());
				tempMatr = ficoCotr[ts][acti_tv].transpose() * tempMatr;
				//
				tranList.clear();
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(tempMatr, ti); iterStif; ++ iterStif){
						tranList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				globTran[ts][tv].resize(
					baseReco[multGrid.size()] + contNumb, nodeCont[ts][tv].size()
				);
				globTran[ts][tv].setFromTriplets(tranList.begin(), tranList.end());
				globTran_pena[ts][tv] = penaFact_n[ts] * globTran[ts][tv];
			}
			else{
				//
				std::vector<Eigen::Triplet<double>> tranList, tranList_1;
				tranList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				tranList_1.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::Matrix<double,3,12> N_0e;
					N_0e << 
					M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 0.0, 
					0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 
					0.0, 0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3);
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig 
						* N_0e.transpose() * tempTang.transpose() * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								long col_tm = iter_Cm->second;
								for(long tn = 0; tn < 3; tn ++){
									tranList.emplace_back(3 * row_tj + tk, 
										3 * col_tm + tn, matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					matr_0 = iterInpo.quadWeig 
						* N_0e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								long col_tm = iter_Cm->second;
								for(long tn = 0; tn < 3; tn ++){
									tranList_1.emplace_back(3 * row_tj + tk, 
										3 * col_tm + tn, matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr;
				tempMatr.resize(3 * nodeCont[ts][acti_tv].size(), 3 * nodeCont[ts][tv].size());
				tempMatr.setFromTriplets(tranList.begin(), tranList.end());
				tempMatr = ficoCotr[ts][acti_tv].transpose() * tempMatr;
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr_1;
				tempMatr_1.resize(3 * nodeCont[ts][acti_tv].size(), 3 * nodeCont[ts][tv].size());
				tempMatr_1.setFromTriplets(tranList_1.begin(), tranList_1.end());
				tempMatr_1 = ficoCotr[ts][acti_tv].transpose() * tempMatr_1;
				//
				tranList.clear();
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(tempMatr, ti); iterStif; ++ iterStif){
						tranList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				tranList_1.clear();
				for(long ti = 0; ti < tempMatr_1.rows(); ti ++){
					for(RSPA_INNE iterStif(tempMatr_1, ti); iterStif; ++ iterStif){
						tranList_1.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				globTran[ts][tv].resize(
					baseReco[multGrid.size()] + contNumb, 3 * nodeCont[ts][tv].size()
				);
				globTran[ts][tv].setFromTriplets(tranList.begin(), tranList.end());
				globTran_pena[ts][tv].resize(
					baseReco[multGrid.size()] + contNumb, 3 * nodeCont[ts][tv].size()
				);
				globTran_pena[ts][tv].setFromTriplets(tranList_1.begin(), tranList_1.end());
			}
		}
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE globTran_D");
	globTran_D.resize(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		globTran_D[ts].resize(2);
		long acti_tv = 0;
		for(long tv = 0; tv < 2; tv ++){
			if(fricCoef[ts] == 0.0){
				//
				std::vector<Eigen::Triplet<double>> tranList;
				tranList.reserve(searCont[ts].intePoin.size() * 4 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,1,3> tempNorm = iterInpo.basiVect[0].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig 
						* penaFact_n[ts] * M_0e.transpose() * tempNorm * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_0 = matr_0 * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 4; tk ++){
							long node_tk = iterInpo.node[tv][tk];
							for(long tm = 0; tm < 3; tm ++){
								long free_tm = 3 * node_tk + tm;
								tranList.emplace_back(row_tj, free_tm, matr_0(tj, 3 * tk + tm));
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr;
				tempMatr.resize(nodeCont[ts][acti_tv].size(), tempNumb);
				tempMatr.setFromTriplets(tranList.begin(), tranList.end());
				tempMatr = ficoCotr[ts][acti_tv].transpose() * tempMatr;
				//
				tranList.clear();
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(tempMatr, ti); iterStif; ++ iterStif){
						tranList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				globTran_D[ts][tv].resize(baseReco[multGrid.size()] + contNumb, tempNumb);
				globTran_D[ts][tv].setFromTriplets(tranList.begin(), tranList.end());
			}
			else{
				//
				std::vector<Eigen::Triplet<double>> tranList;
				tranList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::Matrix<double,3,3> tempTang;
					tempTang.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempTang.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempTang.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,1,4> M_0e;
					M_0e << iterInpo.shapFunc[acti_tv][0], 
						iterInpo.shapFunc[acti_tv][1], 
						iterInpo.shapFunc[acti_tv][2], 
						iterInpo.shapFunc[acti_tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::Matrix<double,3,12> N_0e;
					N_0e << 
					M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 0.0, 
					0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3), 0.0, 
					0.0, 0.0, M_0e(0), 0.0, 0.0, M_0e(1), 0.0, 0.0, M_0e(2), 0.0, 0.0, M_0e(3);
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					Eigen::MatrixXd matr_0 = iterInpo.quadWeig 
						* N_0e.transpose() * tempTang.transpose() * tempPena * tempTang * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_0 = matr_0 * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						auto iter_Cj = nodeCont[ts][acti_tv].find(iterInpo.node[acti_tv][tj]);
						long row_tj = iter_Cj->second;
						for(long tk = 0; tk < 3; tk ++){
							for(long tm = 0; tm < 4; tm ++){
								long node_tm = iterInpo.node[tv][tm];
								for(long tn = 0; tn < 3; tn ++){
									long free_tn = 3 * node_tm + tn;
									tranList.emplace_back(3 * row_tj + tk, free_tn, 
										matr_0(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr;
				tempMatr.resize(3 * nodeCont[ts][acti_tv].size(), tempNumb);
				tempMatr.setFromTriplets(tranList.begin(), tranList.end());
				tempMatr = ficoCotr[ts][acti_tv].transpose() * tempMatr;
				//
				tranList.clear();
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(tempMatr, ti); iterStif; ++ iterStif){
						tranList.emplace_back(baseReco[multGrid.size()] + contReco[ts] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				globTran_D[ts][tv].resize(baseReco[multGrid.size()] + contNumb, tempNumb);
				globTran_D[ts][tv].setFromTriplets(tranList.begin(), tranList.end());
			}
		}
	}
	return 1;
}

long MCONTACT::DOUBLE_M(VECTOR3L coarNode){
	mgpi.maxiLeve = *(std::max_element(doleMcsc.begin(), doleMcsc.end()));
	mgpi.consStif.resize(mgpi.maxiLeve + 1);
	mgpi.realProl.resize(mgpi.maxiLeve);
	mgpi.consStif[mgpi.maxiLeve] = globCoup;
	for(long tl = 1; tl <= mgpi.maxiLeve; tl ++){
		SPARMATS ficoCotr;
		VECT_RESI(ficoCotr, searCont.size(), 2);
		#pragma omp parallel for
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(long tv = 0; tv < 2; tv ++){
				//
				const auto &tempGrid = multGrid[contBody[ts][tv]];
				long currLeve = doleMcsc[contBody[ts][tv]] - tl;
				if(currLeve >= 0){
					const auto &origProl = tempGrid.scalProl[currLeve];
					std::vector<Eigen::Triplet<double>> tempList;
					for(long ti = 0; ti < coarNode[ts][tv].size(); ti ++){
						long temp_ti = coarNode[ts][tv][ti];
						for(RSPA_INNE iterMatr(origProl, temp_ti); iterMatr; ++ iterMatr){
							tempList.emplace_back(ti, iterMatr.col(), iterMatr.value());
						}
					}
					Eigen::SparseMatrix<double,Eigen::ColMajor> coluMatr;
					coluMatr.resize(coarNode[ts][tv].size(), origProl.cols());
					coluMatr.setFromTriplets(tempList.begin(), tempList.end());
					coarNode[ts][tv].clear();
					long tempColu = 0;
					std::vector<Eigen::Triplet<double>> tempList_f;
					if(fricCoef[ts] == 0.0){
						for(long ti = 0; ti < coluMatr.cols(); ti ++){
							long coluTota = 0;
							for(CSPA_INNE iterStif(coluMatr, ti); iterStif; ++ iterStif){
								tempList_f.emplace_back(iterStif.row(), 
									tempColu, iterStif.value()
								);
								coluTota ++;
							}
							if(coluTota > 0){
								coarNode[ts][tv].emplace_back(ti);
								tempColu ++;
							}
						}
						ficoCotr[ts][tv].resize(coluMatr.rows(), tempColu);
					}
					else{
						for(long ti = 0; ti < coluMatr.cols(); ti ++){
							long coluTota = 0;
							for(CSPA_INNE iterStif(coluMatr, ti); iterStif; ++ iterStif){
								for(long tj = 0; tj < 3; tj ++){
									tempList_f.emplace_back(3 * iterStif.row() + tj, 
										3 * tempColu + tj, iterStif.value()
									);
								}
								coluTota ++;
							}
							if(coluTota > 0){
								coarNode[ts][tv].emplace_back(ti);
								tempColu ++;
							}
						}
						ficoCotr[ts][tv].resize(3 * coluMatr.rows(), 3 * tempColu);
					}
					ficoCotr[ts][tv].setFromTriplets(tempList_f.begin(), tempList_f.end());
				}
				else{
					std::vector<Eigen::Triplet<double>> tempList_f;
					if(fricCoef[ts] == 0.0){
						for(long ti = 0; ti < coarNode[ts][tv].size(); ti ++){
							tempList_f.emplace_back(ti, ti, 1.0);
						}
						ficoCotr[ts][tv].resize(coarNode[ts][tv].size(), 
							coarNode[ts][tv].size()
						);
					}
					else{
						for(long ti = 0; ti < 3 * coarNode[ts][tv].size(); ti ++){
							tempList_f.emplace_back(ti, ti, 1.0);
						}
						ficoCotr[ts][tv].resize(3 * coarNode[ts][tv].size(), 
							3 * coarNode[ts][tv].size()
						);
					}
					ficoCotr[ts][tv].setFromTriplets(tempList_f.begin(), tempList_f.end());
				}
			}
		}
		//
		std::vector<Eigen::Triplet<double>> tempList;
		long freeNumb_0 = 0;
		long freeNumb_1 = 0;
		for(long tv = 0; tv < multGrid.size(); tv ++){
			long tempLeve = doleMcsc[tv] - tl;
			if(tempLeve >= 0){
				const auto &tempMatr = multGrid[tv].mgpi.realProl[tempLeve];
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterMatr(tempMatr, ti); iterMatr; ++ iterMatr){
						tempList.emplace_back(
							freeNumb_0 + ti, freeNumb_1 + iterMatr.col(), iterMatr.value()
						);
					}
				}
				freeNumb_0 += multGrid[tv].mgpi.consStif[tempLeve + 1].rows();
			}
			else{
				tempLeve = 0;
				for(long ti = 0; ti < multGrid[tv].mgpi.consStif[tempLeve].rows(); ti ++){
					tempList.emplace_back(freeNumb_0 + ti, freeNumb_1 + ti, 1.0);
				}
				freeNumb_0 += multGrid[tv].mgpi.consStif[tempLeve].rows();
			}
			freeNumb_1 += multGrid[tv].mgpi.consStif[tempLeve].rows();
		}
		for(long ts = 0; ts < searCont.size(); ts ++){
			long acti_tv = 0;
			for(long ti = 0; ti < ficoCotr[ts][acti_tv].rows(); ti ++){
				for(RSPA_INNE iterMatr(ficoCotr[ts][acti_tv], ti); iterMatr; ++ iterMatr){
					tempList.emplace_back(
						freeNumb_0 + ti, freeNumb_1 + iterMatr.col(), iterMatr.value()
					);
				}
			}
			freeNumb_0 += ficoCotr[ts][acti_tv].rows();
			freeNumb_1 += ficoCotr[ts][acti_tv].cols();
		}
		long currLeve = mgpi.maxiLeve - tl;
		mgpi.realProl[currLeve].resize(freeNumb_0, freeNumb_1);
		mgpi.realProl[currLeve].setFromTriplets(tempList.begin(), tempList.end());
		mgpi.consStif[currLeve] = mgpi.realProl[currLeve].transpose() 
			* mgpi.consStif[currLeve + 1] * mgpi.realProl[currLeve];
	}
	return 1;
}

long MCONTACT::MULTISCALE_1(){
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE_1 globCoup_1");
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> domaCoup(multGrid.size());
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		std::vector<Eigen::Triplet<double>> globList;
		const auto &tempStif = multGrid[tv].mgpi.consStif[doleMcsc[tv]];
		for(long ti = 0; ti < tempStif.rows(); ti ++){
			for(RSPA_INNE iterStif(tempStif, ti); iterStif; ++ iterStif){
				globList.emplace_back(baseReco[tv] + ti, 
					baseReco[tv] + iterStif.col(), iterStif.value()
				);
			}
		}
		domaCoup[tv].resize(baseReco[multGrid.size()], baseReco[multGrid.size()]);
		domaCoup[tv].setFromTriplets(globList.begin(), globList.end());
		std::cout << " " << tv << " ";
	}
	globCoup_1 = domaCoup[0];
	for(long tv = 1; tv < multGrid.size(); tv ++){
		globCoup_1 += domaCoup[tv];
	}
	//
	SPARMATS inteCoup;
	VECT_RESI(inteCoup, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			std::vector<Eigen::Triplet<double>> globList;
			std::vector<Eigen::Triplet<double>> selfList, mateList;
			selfList.reserve(searCont[ts].intePoin.size() * 12 * 12);
			mateList.reserve(searCont[ts].intePoin.size() * 12 * 12);
			for(const auto &iterInpo : searCont[ts].intePoin){
				Eigen::MatrixXd tempNorm;
				if(fricCoef[ts] == 0.0){
					tempNorm = iterInpo.basiVect[0].transpose();
				}
				else{
					tempNorm.resize(3,3);
					tempNorm.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempNorm.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempNorm.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
				}
				Eigen::Matrix<double,1,4> M_e, M_m;
				M_e << iterInpo.shapFunc[tv][0], 
					iterInpo.shapFunc[tv][1], 
					iterInpo.shapFunc[tv][2], 
					iterInpo.shapFunc[tv][3];
				M_m << iterInpo.shapFunc[1 - tv][0], 
					iterInpo.shapFunc[1 - tv][1], 
					iterInpo.shapFunc[1 - tv][2], 
					iterInpo.shapFunc[1 - tv][3];
				Eigen::Matrix<double,3,12> N_e, N_m;
				N_e << 
					M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
					0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
					0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
				N_m << 
					M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 0.0, 
					0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 
					0.0, 0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3);
				Eigen::MatrixXd matr_e;
				if(fricCoef[ts] == 0.0){
					matr_e = - 0.5 * iterInpo.quadWeig * penaFact_n[ts] 
						* N_e.transpose() * tempNorm.transpose() * tempNorm * N_e;
				}
				else{
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					matr_e = - 0.5 * iterInpo.quadWeig 
						* N_e.transpose() * tempNorm.transpose() * tempPena * tempNorm * N_e;
				}
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
				CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
				matr_e = tempRota.transpose() * matr_e * tempRota;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = iterInpo.node[tv][tj];
					for(long tk = 0; tk < 3; tk ++){
						long free_tk = 3 * node_tj + tk;
						for(long tm = 0; tm < 4; tm ++){
							long node_tm = iterInpo.node[tv][tm];
							for(long tn = 0; tn < 3; tn ++){
								long free_tn = 3 * node_tm + tn;
								selfList.emplace_back(free_tk, free_tn, 
									matr_e(3 * tj + tk, 3 * tm + tn)
								);
							}
						}
					}
				}
				Eigen::MatrixXd matr_m;
				if(fricCoef[ts] == 0.0){
					matr_m = - 0.5 * iterInpo.quadWeig * penaFact_n[ts] 
						* N_e.transpose() * tempNorm.transpose() * tempNorm * N_m;
				}
				else{
					Eigen::Matrix3d tempPena;
					tempPena << penaFact_n[ts], 0.0, 0.0, 
						0.0, penaFact_f[ts], 0.0, 
						0.0, 0.0, penaFact_f[ts];
					matr_m = - 0.5 * iterInpo.quadWeig 
						* N_e.transpose() * tempNorm.transpose() * tempPena * tempNorm * N_m;
				}
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota_m;
				CONT_ROTA(multGrid[contBody[ts][1 - tv]], iterInpo.node[1 - tv], tempRota_m);
				matr_m = tempRota.transpose() * matr_m * tempRota_m;
				for(long tj = 0; tj < 4; tj ++){
					long node_tj = iterInpo.node[tv][tj];
					for(long tk = 0; tk < 3; tk ++){
						long free_tk = 3 * node_tj + tk;
						for(long tm = 0; tm < 4; tm ++){
							long node_tm = iterInpo.node[1 - tv][tm];
							for(long tn = 0; tn < 3; tn ++){
								long free_tn = 3 * node_tm + tn;
								mateList.emplace_back(free_tk, free_tn, 
									matr_m(3 * tj + tk, 3 * tm + tn)
								);
							}
						}
					}
				}
			}
			long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
			Eigen::SparseMatrix<double,Eigen::RowMajor> selfMatr;
			selfMatr.resize(tempNumb, tempNumb);
			selfMatr.setFromTriplets(selfList.begin(), selfList.end());
			long tempNumb_m = 3 * multGrid[contBody[ts][1 - tv]].nodeCoor.size();
			Eigen::SparseMatrix<double,Eigen::RowMajor> mateMatr;
			mateMatr.resize(tempNumb, tempNumb_m);
			mateMatr.setFromTriplets(mateList.begin(), mateList.end());
			//
			const auto &tempGrid = multGrid[contBody[ts][tv]];
			long tempLeve = doleMcsc[contBody[ts][tv]];
			selfMatr = 
				tempGrid.earlTran.transpose() * selfMatr * tempGrid.earlTran;
			mateMatr = 
				tempGrid.earlTran.transpose() * mateMatr;
			for(long ti = tempGrid.mgpi.maxiLeve; ti >= tempLeve; ti --){
				selfMatr = tempGrid.prolOper[ti].transpose() * selfMatr 
					* tempGrid.prolOper[ti];
				mateMatr = tempGrid.prolOper[ti].transpose() * mateMatr;
			}
			selfMatr = tempGrid.consOper[tempLeve] * selfMatr 
				* tempGrid.consOper[tempLeve].transpose();
			mateMatr = tempGrid.consOper[tempLeve] * mateMatr;
			for(long ti = 0; ti < selfMatr.rows(); ti ++){
				for(RSPA_INNE iterStif(selfMatr, ti); iterStif; ++ iterStif){
					globList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
						baseReco[contBody[ts][tv]] + iterStif.col(), iterStif.value()
					);
				}
			}
			//
			const auto &tempGrid_m = multGrid[contBody[ts][1 - tv]];
			long tempLeve_m = doleMcsc[contBody[ts][1 - tv]];
			mateMatr = mateMatr * tempGrid_m.earlTran;
			for(long ti = tempGrid_m.mgpi.maxiLeve; ti >= tempLeve_m; ti --){
				mateMatr = mateMatr * tempGrid_m.prolOper[ti];
			}
			mateMatr = mateMatr * tempGrid_m.consOper[tempLeve_m].transpose();
			for(long ti = 0; ti < mateMatr.rows(); ti ++){
				for(RSPA_INNE iterStif(mateMatr, ti); iterStif; ++ iterStif){
					globList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
						baseReco[contBody[ts][1 - tv]] + iterStif.col(), iterStif.value()
					);
				}
			}
			inteCoup[ts][tv].resize(baseReco[multGrid.size()], baseReco[multGrid.size()]);
			inteCoup[ts][tv].setFromTriplets(globList.begin(), globList.end());
			std::cout << " " << ts << "-" << tv << " ";
		}
	}
	std::cout << std::endl;
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			globCoup_1 += inteCoup[ts][tv];
		}
	}
	// globCoup_1.resize(baseReco[multGrid.size()], baseReco[multGrid.size()]);
	// globCoup_1.setFromTriplets(globList.begin(), globList.end());
	std::cout << "Dimension of globCoup_1: " << globCoup_1.rows();
	OUTPUT_TIME("");
	if(globCoup_1.rows() < DIRE_MAXI){
		coarSolv_D_1.compute(globCoup_1);
	}
	else if(globCoup_1.rows() < COGR_MAXI){
		coarSolv_C_1.compute(globCoup_1);
	}
	else{
		DOUBLE_M_1();
		mgpi_1.ESTABLISH();
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE_1 globTran_D_1");
	globTran_D_1.resize(multGrid.size());
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> neceTran(multGrid.size());
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		const auto &tempGrid = multGrid[tv];
		neceTran[tv] = tempGrid.earlTran.transpose();
		const auto &tempProl = tempGrid.prolOper[tempGrid.mgpi.maxiLeve];
		long tempCols = neceTran[tv].cols();
		neceTran[tv] = 
			tempProl.block(0,0,tempProl.cols(),tempProl.cols()).transpose()
			* neceTran[tv].block(0,0,tempProl.cols(),tempCols).eval();
		neceTran[tv] = tempGrid.consOper[tempGrid.mgpi.maxiLeve] * neceTran[tv].eval();
		Eigen::SparseMatrix<double,Eigen::RowMajor> tempStif = 
			multGrid[tv].mgpi.consStif[tempGrid.mgpi.maxiLeve] * neceTran[tv];
		for(long ti = tempGrid.mgpi.maxiLeve - 1; ti >= doleMcsc[tv]; ti --){
			tempStif = tempGrid.mgpi.realProl[ti].transpose() * tempStif.eval();
		}
		std::vector<Eigen::Triplet<double>> globList;
		for(long ti = 0; ti < tempStif.rows(); ti ++){
			for(RSPA_INNE iterStif(tempStif, ti); iterStif; ++ iterStif){
				globList.emplace_back(baseReco[tv] + ti, iterStif.col(), iterStif.value());
			}
		}
		long tempNumb = 3 * tempGrid.nodeCoor.size();
		globTran_D_1[tv].resize(baseReco[multGrid.size()], tempNumb);
		globTran_D_1[tv].setFromTriplets(globList.begin(), globList.end());
	}
	//
	#pragma omp parallel for
	for(long tv_real = 0; tv_real < multGrid.size(); tv_real ++){
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(long tv = 0; tv < 2; tv ++){
				if(tv_real != contBody[ts][tv]){
					continue;
				}
				std::vector<Eigen::Triplet<double>> selfList, mateList;
				selfList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				mateList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				for(const auto &iterInpo : searCont[ts].intePoin){
					Eigen::MatrixXd tempNorm;
					if(fricCoef[ts] == 0.0){
						tempNorm = iterInpo.basiVect[0].transpose();
					}
					else{
						tempNorm.resize(3,3);
						tempNorm.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
						tempNorm.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
						tempNorm.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
					}
					Eigen::Matrix<double,1,4> M_e, M_m;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					M_m << iterInpo.shapFunc[1 - tv][0], 
						iterInpo.shapFunc[1 - tv][1], 
						iterInpo.shapFunc[1 - tv][2], 
						iterInpo.shapFunc[1 - tv][3];
					Eigen::Matrix<double,3,12> N_e, N_m;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					N_m << 
						M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 0.0, 
						0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 
						0.0, 0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3);
					Eigen::MatrixXd matr_e;
					if(fricCoef[ts] == 0.0){
						matr_e = - 0.5 * iterInpo.quadWeig * penaFact_n[ts] 
							* N_e.transpose() * tempNorm.transpose() * tempNorm * N_e;
					}
					else{
						Eigen::Matrix3d tempPena;
						tempPena << penaFact_n[ts], 0.0, 0.0, 
							0.0, penaFact_f[ts], 0.0, 
							0.0, 0.0, penaFact_f[ts];
						matr_e = - 0.5 * iterInpo.quadWeig 
							* N_e.transpose() * tempNorm.transpose() * tempPena * tempNorm * N_e;
					}
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_e = tempRota.transpose() * matr_e * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								long node_tm = iterInpo.node[tv][tm];
								for(long tn = 0; tn < 3; tn ++){
									long free_tn = 3 * node_tm + tn;
									selfList.emplace_back(free_tk, free_tn, 
										matr_e(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
					Eigen::MatrixXd matr_m;
					if(fricCoef[ts] == 0.0){
						matr_m = - 0.5 * iterInpo.quadWeig * penaFact_n[ts] 
							* N_m.transpose() * tempNorm.transpose() * tempNorm * N_e;
					}
					else{
						Eigen::Matrix3d tempPena;
						tempPena << penaFact_n[ts], 0.0, 0.0, 
							0.0, penaFact_f[ts], 0.0, 
							0.0, 0.0, penaFact_f[ts];
						matr_m = - 0.5 * iterInpo.quadWeig 
							* N_m.transpose() * tempNorm.transpose() * tempPena * tempNorm * N_e;
					}
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota_m;
					CONT_ROTA(multGrid[contBody[ts][1 - tv]], iterInpo.node[1 - tv], tempRota_m);
					matr_m = tempRota_m.transpose() * matr_m * tempRota;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[1 - tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								long node_tm = iterInpo.node[tv][tm];
								for(long tn = 0; tn < 3; tn ++){
									long free_tn = 3 * node_tm + tn;
									mateList.emplace_back(free_tk, free_tn, 
										matr_m(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
				long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> selfMatr;
				selfMatr.resize(tempNumb, tempNumb);
				selfMatr.setFromTriplets(selfList.begin(), selfList.end());
				long tempNumb_m = 3 * multGrid[contBody[ts][1 - tv]].nodeCoor.size();
				Eigen::SparseMatrix<double,Eigen::RowMajor> mateMatr;
				mateMatr.resize(tempNumb_m, tempNumb);
				mateMatr.setFromTriplets(mateList.begin(), mateList.end());
				//
				const auto &tempGrid = multGrid[contBody[ts][tv]];
				long tempLeve = doleMcsc[contBody[ts][tv]];
				selfMatr = tempGrid.earlTran.transpose() * selfMatr.eval() * tempGrid.earlTran;
				selfMatr = 
					tempGrid.prolOper[tempGrid.mgpi.maxiLeve].transpose() * selfMatr.eval() 
					* tempGrid.prolOper[tempGrid.mgpi.maxiLeve];
				selfMatr = tempGrid.consOper[tempGrid.mgpi.maxiLeve] * selfMatr.eval() 
					* tempGrid.consOper[tempGrid.mgpi.maxiLeve].transpose();
				selfMatr = selfMatr.eval() * neceTran[contBody[ts][tv]];
				for(long ti = tempGrid.mgpi.maxiLeve - 1; ti >= tempLeve; ti --){
					selfMatr = tempGrid.mgpi.realProl[ti].transpose() * selfMatr.eval();
				}
				std::vector<Eigen::Triplet<double>> globList;
				for(long ti = 0; ti < selfMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(selfMatr, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				//
				const auto &tempGrid_m = multGrid[contBody[ts][1 - tv]];
				long tempLeve_m = doleMcsc[contBody[ts][1 - tv]];
				mateMatr = tempGrid_m.earlTran.transpose() * mateMatr.eval() * tempGrid.earlTran;
				mateMatr = 
					tempGrid_m.prolOper[tempGrid_m.mgpi.maxiLeve].transpose() * mateMatr.eval() 
					* tempGrid.prolOper[tempGrid.mgpi.maxiLeve];
				mateMatr = tempGrid_m.consOper[tempGrid_m.mgpi.maxiLeve] * mateMatr.eval() 
					* tempGrid.consOper[tempGrid.mgpi.maxiLeve].transpose();
				mateMatr = mateMatr.eval() * neceTran[contBody[ts][tv]];
				for(long ti = tempGrid_m.mgpi.maxiLeve - 1; ti >= tempLeve_m; ti --){
					mateMatr = tempGrid_m.mgpi.realProl[ti].transpose() * mateMatr.eval();
				}
				for(long ti = 0; ti < mateMatr.rows(); ti ++){
					for(RSPA_INNE iterStif(mateMatr, ti); iterStif; ++ iterStif){
						globList.emplace_back(baseReco[contBody[ts][1 - tv]] + ti, 
							iterStif.col(), iterStif.value()
						);
					}
				}
				//
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr;
				tempMatr.resize(baseReco[multGrid.size()], tempNumb);
				tempMatr.setFromTriplets(globList.begin(), globList.end());
				globTran_D_1[contBody[ts][tv]] += tempMatr;
			}
		}
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE_1 globForc_1");
	globForc_1.resize(baseReco[multGrid.size()]);
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		//
		Eigen::VectorXd tempCofo = multGrid[tv].consForc;
		for(long ti = multGrid[tv].mgpi.maxiLeve - 1; ti >= doleMcsc[tv]; ti --){
			tempCofo = multGrid[tv].mgpi.realProl[ti].transpose() * tempCofo.eval();
		}
		globForc_1.block(baseReco[tv],0,tempCofo.rows(),1) = tempCofo;
	}
	#pragma omp parallel for
	for(long tv_real = 0; tv_real < multGrid.size(); tv_real ++){
		for(long ts = 0; ts < searCont.size(); ts ++){
			// if(fricCoef[ts] < 0.0){
				// continue;//for perfect interface, the initial gap is always 0!
			// }
			for(long tv = 0; tv < 2; tv ++){
				if(tv_real != contBody[ts][tv]){
					continue;
				}
				Eigen::VectorXd tempGafo = Eigen::VectorXd::Zero(
					3 * multGrid[contBody[ts][tv]].nodeCoor.size()
				);
				for(const auto &iterInpo : searCont[ts].intePoin){
					//only normal gap, no tangential gap
					Eigen::MatrixXd tempNorm = iterInpo.basiVect[0].transpose();
					Eigen::Matrix<double,1,4> M_e;
					M_e << iterInpo.shapFunc[tv][0], 
						iterInpo.shapFunc[tv][1], 
						iterInpo.shapFunc[tv][2], 
						iterInpo.shapFunc[tv][3];
					Eigen::Matrix<double,3,12> N_e;
					N_e << 
						M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
						0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
						0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
					Eigen::MatrixXd matr_e = 0.5 * iterInpo.quadWeig * penaFact_n[ts] 
						* N_e.transpose() * tempNorm.transpose() * iterInpo.initNgap;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_e = tempRota.transpose() * matr_e;
					if(tv == 1){
						matr_e = - matr_e;
					}
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							tempGafo(free_tk) += matr_e(3 * tj + tk, 0);
						}
					}
				}
				//
				const auto &tempGrid = multGrid[contBody[ts][tv]];
				long tempLeve = doleMcsc[contBody[ts][tv]];
				tempGafo = tempGrid.earlTran.transpose() * tempGafo;
				tempGafo = tempGrid.prolOper[tempGrid.mgpi.maxiLeve].transpose() * tempGafo;
				tempGafo = tempGrid.consOper[tempGrid.mgpi.maxiLeve] * tempGafo;
				for(long ti = tempGrid.mgpi.maxiLeve - 1; ti >= tempLeve; ti --){
					tempGafo = tempGrid.mgpi.realProl[ti].transpose() * tempGafo;
				}
				globForc_1.block(baseReco[contBody[ts][tv]],0,tempGafo.rows(),1) += tempGafo;
			}
		}
	}
	//***************************************************************************************
	OUTPUT_TIME("MCONTACT::MULTISCALE_1 globTran_1");
	VECT_RESI(globTran_1, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			//
			std::vector<Eigen::Triplet<double>> selfList, mateList;
			if(fricCoef[ts] == 0.0){
				selfList.reserve(searCont[ts].intePoin.size() * 12 * 4);
				mateList.reserve(searCont[ts].intePoin.size() * 12 * 4);
			}
			else{
				selfList.reserve(searCont[ts].intePoin.size() * 12 * 12);
				mateList.reserve(searCont[ts].intePoin.size() * 12 * 12);
			}
			for(const auto &iterInpo : searCont[ts].intePoin){
				Eigen::MatrixXd tempNorm;
				if(fricCoef[ts] == 0.0){
					tempNorm = iterInpo.basiVect[0].transpose();
				}
				else{
					tempNorm.resize(3,3);
					tempNorm.block(0,0,1,3) = iterInpo.basiVect[0].transpose();
					tempNorm.block(1,0,1,3) = iterInpo.basiVect[1].transpose();
					tempNorm.block(2,0,1,3) = iterInpo.basiVect[2].transpose();
				}
				Eigen::Matrix<double,1,4> M_e;
				M_e << iterInpo.shapFunc[tv][0], 
					iterInpo.shapFunc[tv][1], 
					iterInpo.shapFunc[tv][2], 
					iterInpo.shapFunc[tv][3];
				Eigen::Matrix<double,1,4> M_m;
				M_m << iterInpo.shapFunc[1 - tv][0], 
					iterInpo.shapFunc[1 - tv][1], 
					iterInpo.shapFunc[1 - tv][2], 
					iterInpo.shapFunc[1 - tv][3];
				Eigen::Matrix<double,3,12> N_e, N_m;
				N_e << 
					M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 0.0, 
					0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3), 0.0, 
					0.0, 0.0, M_e(0), 0.0, 0.0, M_e(1), 0.0, 0.0, M_e(2), 0.0, 0.0, M_e(3);
				N_m << 
					M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 0.0, 
					0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3), 0.0, 
					0.0, 0.0, M_m(0), 0.0, 0.0, M_m(1), 0.0, 0.0, M_m(2), 0.0, 0.0, M_m(3);
				if(fricCoef[ts] == 0.0){
					Eigen::MatrixXd matr_e = - 0.5 * iterInpo.quadWeig 
						* N_e.transpose() * tempNorm.transpose() * M_e;
					Eigen::MatrixXd matr_m = 0.5 * iterInpo.quadWeig 
						* N_m.transpose() * tempNorm.transpose() * M_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_e = tempRota.transpose() * matr_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota_m;
					CONT_ROTA(multGrid[contBody[ts][1 - tv]], iterInpo.node[1 - tv], tempRota_m);
					matr_m = tempRota_m.transpose() * matr_m;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								selfList.emplace_back(free_tk, iter_Cm->second, 
									matr_e(3 * tj + tk, tm)
								);
							}
						}
					}
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[1 - tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								mateList.emplace_back(free_tk, iter_Cm->second, 
									matr_m(3 * tj + tk, tm)
								);
							}
						}
					}
				}
				else{
					Eigen::MatrixXd matr_e = - 0.5 * iterInpo.quadWeig 
						* N_e.transpose() * tempNorm.transpose() * tempNorm * N_e;
					Eigen::MatrixXd matr_m = 0.5 * iterInpo.quadWeig 
						* N_m.transpose() * tempNorm.transpose() * tempNorm * N_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota;
					CONT_ROTA(multGrid[contBody[ts][tv]], iterInpo.node[tv], tempRota);
					matr_e = tempRota.transpose() * matr_e;
					Eigen::SparseMatrix<double,Eigen::RowMajor> tempRota_m;
					CONT_ROTA(multGrid[contBody[ts][1 - tv]], iterInpo.node[1 - tv], tempRota_m);
					matr_m = tempRota_m.transpose() * matr_m;
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								for(long tn = 0; tn < 3; tn ++){
									selfList.emplace_back(free_tk, 3 * iter_Cm->second + tn, 
										matr_e(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
					for(long tj = 0; tj < 4; tj ++){
						long node_tj = iterInpo.node[1 - tv][tj];
						for(long tk = 0; tk < 3; tk ++){
							long free_tk = 3 * node_tj + tk;
							for(long tm = 0; tm < 4; tm ++){
								auto iter_Cm = nodeCont[ts][tv].find(iterInpo.node[tv][tm]);
								for(long tn = 0; tn < 3; tn ++){
									mateList.emplace_back(free_tk, 3 * iter_Cm->second + tn, 
										matr_m(3 * tj + tk, 3 * tm + tn)
									);
								}
							}
						}
					}
				}
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> selfMatr, mateMatr;
			long tempNumb = 3 * multGrid[contBody[ts][tv]].nodeCoor.size();
			long tempNumb_m = 3 * multGrid[contBody[ts][1 - tv]].nodeCoor.size();
			if(fricCoef[ts] == 0.0){
				selfMatr.resize(tempNumb, nodeCont[ts][tv].size());
				mateMatr.resize(tempNumb_m, nodeCont[ts][tv].size());
				globTran_1[ts][tv].resize(baseReco[multGrid.size()], nodeCont[ts][tv].size());
			}
			else{
				selfMatr.resize(tempNumb, 3 * nodeCont[ts][tv].size());
				mateMatr.resize(tempNumb_m, 3 * nodeCont[ts][tv].size());
				globTran_1[ts][tv].resize(baseReco[multGrid.size()], 
					3 * nodeCont[ts][tv].size()
				);
			}
			selfMatr.setFromTriplets(selfList.begin(), selfList.end());
			mateMatr.setFromTriplets(mateList.begin(), mateList.end());
			//
			const auto &tempGrid = multGrid[contBody[ts][tv]];
			long tempLeve = doleMcsc[contBody[ts][tv]];
			selfMatr = tempGrid.earlTran.transpose() * selfMatr.eval();
			selfMatr = tempGrid.prolOper[tempGrid.mgpi.maxiLeve].transpose() * selfMatr.eval();
			selfMatr = tempGrid.consOper[tempGrid.mgpi.maxiLeve] * selfMatr.eval();
			for(long ti = tempGrid.mgpi.maxiLeve - 1; ti >= tempLeve; ti --){
				selfMatr = tempGrid.mgpi.realProl[ti].transpose() * selfMatr.eval();
			}
			std::vector<Eigen::Triplet<double>> tranList;
			for(long ti = 0; ti < selfMatr.rows(); ti ++){
				for(RSPA_INNE iterStif(selfMatr, ti); iterStif; ++ iterStif){
					tranList.emplace_back(baseReco[contBody[ts][tv]] + ti, 
						iterStif.col(), iterStif.value()
					);
				}
			}
			//
			const auto &tempGrid_m = multGrid[contBody[ts][1 - tv]];
			long tempLeve_m = doleMcsc[contBody[ts][1 - tv]];
			mateMatr = tempGrid_m.earlTran.transpose() * mateMatr.eval();
			mateMatr = 
				tempGrid_m.prolOper[tempGrid_m.mgpi.maxiLeve].transpose() * mateMatr.eval();
			mateMatr = tempGrid_m.consOper[tempGrid_m.mgpi.maxiLeve] * mateMatr.eval();
			for(long ti = tempGrid_m.mgpi.maxiLeve - 1; ti >= tempLeve_m; ti --){
				mateMatr = tempGrid_m.mgpi.realProl[ti].transpose() * mateMatr.eval();
			}
			for(long ti = 0; ti < mateMatr.rows(); ti ++){
				for(RSPA_INNE iterStif(mateMatr, ti); iterStif; ++ iterStif){
					tranList.emplace_back(baseReco[contBody[ts][1 - tv]] + ti, 
						iterStif.col(), iterStif.value()
					);
				}
			}
			globTran_1[ts][tv].setFromTriplets(tranList.begin(), tranList.end());
		}
	}
	return 1;
}

long MCONTACT::DOUBLE_M_1(){
	mgpi_1.maxiLeve = *(std::max_element(doleMcsc.begin(), doleMcsc.end()));
	mgpi_1.consStif.resize(mgpi_1.maxiLeve + 1);
	mgpi_1.realProl.resize(mgpi_1.maxiLeve);
	mgpi_1.consStif[mgpi_1.maxiLeve] = globCoup_1;
	for(long tl = 1; tl <= mgpi_1.maxiLeve; tl ++){
		std::vector<Eigen::Triplet<double>> tempList;
		long freeNumb_0 = 0;
		long freeNumb_1 = 0;
		for(long tv = 0; tv < multGrid.size(); tv ++){
			long tempLeve = doleMcsc[tv] - tl;
			if(tempLeve >= 0){
				const auto &tempMatr = multGrid[tv].mgpi.realProl[tempLeve];
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterMatr(tempMatr, ti); iterMatr; ++ iterMatr){
						tempList.emplace_back(
							freeNumb_0 + ti, freeNumb_1 + iterMatr.col(), iterMatr.value()
						);
					}
				}
				freeNumb_0 += multGrid[tv].mgpi.consStif[tempLeve + 1].rows();
			}
			else{
				tempLeve = 0;
				for(long ti = 0; ti < multGrid[tv].mgpi.consStif[tempLeve].rows(); ti ++){
					tempList.emplace_back(freeNumb_0 + ti, freeNumb_1 + ti, 1.0);
				}
				freeNumb_0 += multGrid[tv].mgpi.consStif[tempLeve].rows();
			}
			freeNumb_1 += multGrid[tv].mgpi.consStif[tempLeve].rows();
		}
		long currLeve = mgpi_1.maxiLeve - tl;
		mgpi_1.realProl[currLeve].resize(freeNumb_0, freeNumb_1);
		mgpi_1.realProl[currLeve].setFromTriplets(tempList.begin(), tempList.end());
		mgpi_1.consStif[currLeve] = mgpi_1.realProl[currLeve].transpose() 
			* mgpi_1.consStif[currLeve + 1] * mgpi_1.realProl[currLeve];
	}
	return 1;
}

long MCONTACT::APPS(){
	OUTPUT_TIME("MCONTACT::APPS");
	if(globCoup_1.rows() <= 0){
		OUTPUT_TIME("MCONTACT::APPS ERROR 1");
		return -1;
	}
	//
	long freqNumb = 10;
	typedef Spectra::SparseSymMatProd<double,Eigen::Lower,Eigen::RowMajor> SSMP;
	SSMP oper(globCoup_1);
	int tempNumb = globCoup_1.rows();
	Spectra::SymEigsSolver<SSMP> eigs_0(oper, freqNumb, std::min(20 * (int)freqNumb, tempNumb));
	eigs_0.init();
	//SmallestMagn,LargestMagn
	int ncon = eigs_0.compute(Spectra::SortRule::SmallestMagn, 
		1000000, 1.0E-6, Spectra::SortRule::SmallestMagn
	);
	if(eigs_0.info() != Spectra::CompInfo::Successful){
		OUTPUT_TIME("MCONTACT::APPS ERROR 2");
		return -1;
	}
	Eigen::VectorXd tempValu = eigs_0.eigenvalues();
	Eigen::MatrixXd tempVect = eigs_0.eigenvectors();
	//
	Eigen::VectorXd tempForc = globForc_1.normalized();
	std::ofstream tempOfst(DIRECTORY("resuFreq.txt"), std::ios::out);
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long tf = 0; tf < freqNumb; tf ++){
		Eigen::VectorXd globSolu = tempVect.col(tf);
		double tempCorr = globSolu.dot(tempForc);
		tempOfst << std::setw(30) << tempValu[tf] << std::setw(30) << tempCorr  << std::endl;
		std::cout << "Oder: " << tf + 1 << ", freq: " << tempValu[tf] 
			<< ", corr: " << tempCorr << std::endl;
	}
	tempOfst.close();
	for(long tf = 0; tf < freqNumb; tf ++){
		Eigen::VectorXd globSolu = tempVect.col(tf);
		std::cout << "Normalization of " << tf << "-th mode: " 
			<< globSolu.dot(globSolu) << std::endl;
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			Eigen::VectorXd resuSolu = globSolu.block(
				baseReco[tv],0,multGrid[tv].mgpi.consStif[doleMcsc[tv]].rows(),1
			);
			resuSolu = multGrid[tv].consOper[doleMcsc[tv]].transpose() * resuSolu;
			for(long ti = doleMcsc[tv]; ti < multGrid[tv].mgpi.maxiLeve; ti ++){
				resuSolu = multGrid[tv].prolOper[ti] * resuSolu;
			}
			resuSolu = multGrid[tv].consOper[multGrid[tv].mgpi.maxiLeve] * resuSolu;
			Eigen::VectorXd outpDisp;
			multGrid[tv].OUTP_SUB1(resuSolu, outpDisp);
			std::stringstream tempStre;
			tempStre << tf + 1 << "-" << tv;
			multGrid[tv].OUTP_SUB2(outpDisp, tempStre.str());
			tempStre.str("");
			tempStre.clear();
		}
	}
	OUTPUT_TIME("MCONTACT::APPS finished.");
	return 1;
}

long MCONTACT::APPS_MPL(){
	OUTPUT_TIME("MCONTACT::APPS_MPL");
	if(globCoup.rows() <= 0){
		OUTPUT_TIME("MCONTACT::APPS_MPL ERROR 1");
		return -1;
	}
	//
	long freqNumb = 10;
	typedef Spectra::SparseSymMatProd<double,Eigen::Lower,Eigen::RowMajor> SSMP;
	SSMP oper(globCoup);
	int tempNumb = globCoup.rows();
	Spectra::SymEigsSolver<SSMP> eigs_0(oper, freqNumb, std::min(20 * (int)freqNumb, tempNumb));
	eigs_0.init();
	//SmallestMagn,LargestMagn
	int ncon = eigs_0.compute(Spectra::SortRule::SmallestMagn, 
		1000000, 1.0E-6, Spectra::SortRule::SmallestMagn
	);
	if(eigs_0.info() != Spectra::CompInfo::Successful){
		OUTPUT_TIME("MCONTACT::APPS_MPL ERROR 2");
		return -1;
	}
	Eigen::VectorXd tempValu = eigs_0.eigenvalues();
	Eigen::MatrixXd tempVect = eigs_0.eigenvectors();
	//
	globForc_1.resize(baseReco[multGrid.size()]);
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		//
		Eigen::VectorXd tempCofo = multGrid[tv].consForc;
		for(long ti = multGrid[tv].mgpi.maxiLeve - 1; ti >= doleMcsc[tv]; ti --){
			tempCofo = multGrid[tv].mgpi.realProl[ti].transpose() * tempCofo.eval();
		}
		globForc_1.block(baseReco[tv],0,tempCofo.rows(),1) = tempCofo;
	}
	Eigen::VectorXd tempForc = globForc_1.normalized();
	std::ofstream tempOfst(DIRECTORY("resuFreq.txt"), std::ios::out);
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(long tf = 0; tf < freqNumb; tf ++){
		Eigen::VectorXd globSolu = tempVect.col(tf);
		Eigen::VectorXd tempModa = globSolu.head(tempForc.size());
		double tempCorr = tempModa.dot(tempForc);
		tempOfst << std::setw(30) << tempValu[tf] << std::setw(30) << tempCorr  << std::endl;
		std::cout << "Oder: " << tf + 1 << ", freq: " << tempValu[tf] 
			<< ", corr: " << tempCorr << std::endl;
	}
	tempOfst.close();
	for(long tf = 0; tf < freqNumb; tf ++){
		Eigen::VectorXd globSolu = tempVect.col(tf);
		std::cout << "Normalization of " << tf << "-th mode: " 
			<< globSolu.dot(globSolu) << std::endl;
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			Eigen::VectorXd resuSolu = globSolu.block(
				baseReco[tv],0,multGrid[tv].mgpi.consStif[doleMcsc[tv]].rows(),1
			);
			resuSolu = multGrid[tv].consOper[doleMcsc[tv]].transpose() * resuSolu;
			for(long ti = doleMcsc[tv]; ti < multGrid[tv].mgpi.maxiLeve; ti ++){
				resuSolu = multGrid[tv].prolOper[ti] * resuSolu;
			}
			resuSolu = multGrid[tv].consOper[multGrid[tv].mgpi.maxiLeve] * resuSolu;
			Eigen::VectorXd outpDisp;
			multGrid[tv].OUTP_SUB1(resuSolu, outpDisp);
			std::stringstream tempStre;
			tempStre << tf + 1 << "-" << tv;
			multGrid[tv].OUTP_SUB2(outpDisp, tempStre.str());
			tempStre.str("");
			tempStre.clear();
		}
	}
	OUTPUT_TIME("MCONTACT::APPS_MPL finished.");
	return 1;
}

long MCONTACT::CONTACT_ANALYSIS(){
	long moniCycl = 10;
	VECTOR2D moniReco(multGrid.size() + 4 * searCont.size());
	for(long ti = 0; ti < moniReco.size(); ti ++){
		moniReco[ti].resize(moniCycl);
	}
	//
	long tc;
	long maxiIter = 1000;
	std::ofstream tempOfst(DIRECTORY("resuMoni.txt"), std::ios::out);
	tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
	for(tc = 0; tc < maxiIter; tc ++){
		std::cout << "The " << tc << "-th iteration";
		OUTPUT_TIME("");
		std::vector<Eigen::VectorXd> resuDisp_0 = resuDisp;
		std::vector<std::vector<Eigen::VectorXd>> inteAuxi_0 = inteAuxi;
		std::vector<std::vector<Eigen::VectorXd>> inteLagr_0 = inteLagr;
		//*************************************body blance***********************************
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			//
			Eigen::VectorXd addiForc = Eigen::VectorXd::Zero(3 * multGrid[tv].nodeCoor.size());
			for(long ts = 0; ts < searCont.size(); ts ++){
				for(long ti = 0; ti < 2; ti ++){
					if(tv != contBody[ts][ti]){
						continue;
					}
					addiForc += systTran_pena[ts][ti] * inteAuxi[ts][ti] 
						- systTran[ts][ti] * inteLagr[ts][ti];
				}
			}
			multGrid[tv].ADDITIONAL_FORCE(addiForc);
			//solving
			Eigen::VectorXd resuSolu;
			if(addiForc.rows() < DIRE_MAXI_SUBD){
				resuSolu = (mugrDiso[tv]).solve(multGrid[tv].consForc + addiForc);
			}
			else{
				multGrid[tv].mgpi.CG_SOLV(1, multGrid[tv].consForc + addiForc, resuSolu);
			}
			multGrid[tv].OUTP_SUB1(resuSolu, resuDisp[tv]);
			if(((muscSett >> 0) % 2 == 1 || (muscSett >> 1) % 2 == 1) && tc <= MULT_MAXI){
				continue;
			}
			multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
		}
		//**************************************macroscopic***********************************
		if((muscSett >> 0) % 2 == 1 && tc <= MULT_MAXI){
			Eigen::VectorXd globForc = Eigen::VectorXd::Zero(globCoup.rows());
			for(long ts = 0; ts < searCont.size(); ts ++){
				for(long tv = 0; tv < 2; tv ++){
					globForc += 
						globTran[ts][tv] * inteLagr[ts][tv] 
						- globTran_pena[ts][tv] * inteAuxi[ts][tv] 
						+ globTran_D[ts][tv] * resuDisp[contBody[ts][tv]];
				}
			}
			Eigen::VectorXd globSolu;
			OUTPUT_TIME("MCONTACT::CONTACT_ANALYSIS: coupling solver");
			if(globCoup.rows() < DIRE_MAXI){
				globSolu = coarSolv_D.solve(globForc);
			}
			else if(globCoup.rows() < COGR_MAXI){
				globSolu = coarSolv_C.solve(globForc);
			}
			else{
				mgpi.CG_SOLV(1, globForc, globSolu);
			}
			OUTPUT_TIME("MCONTACT::CONTACT_ANALYSIS: coupling solver done");
			#pragma omp parallel for
			for(long tv = 0; tv < multGrid.size(); tv ++){
				Eigen::VectorXd resuSolu = globSolu.block(
					baseReco[tv],0,multGrid[tv].mgpi.consStif[doleMcsc[tv]].rows(),1
				);
				resuSolu = accuProl[tv] * resuSolu;
				Eigen::VectorXd outpDisp;
				multGrid[tv].OUTP_SUB1(resuSolu, outpDisp);
				resuDisp[tv] += outpDisp;
				multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
			}
		}
		//**********************************interface-eliminated*****************************
		if((muscSett >> 1) % 2 == 1 && tc <= MULT_MAXI){
			Eigen::VectorXd globForc = globForc_1;
			for(long ts = 0; ts < searCont.size(); ts ++){
				for(long tv = 0; tv < 2; tv ++){
					globForc += globTran_1[ts][tv] * inteLagr[ts][tv];
				}
			}
			for(long tv = 0; tv < multGrid.size(); tv ++){
				globForc -= (globTran_D_1[tv] * resuDisp[tv]);
			}
			Eigen::VectorXd globSolu;
			OUTPUT_TIME("MCONTACT::CONTACT_ANALYSIS: coupling solver 1");
			if(globCoup_1.rows() < DIRE_MAXI){
				globSolu = coarSolv_D_1.solve(globForc);
			}
			else if(globCoup_1.rows() < COGR_MAXI){
				globSolu = coarSolv_C_1.solve(globForc);
			}
			else{
				mgpi_1.CG_SOLV(1, globForc, globSolu);
			}
			OUTPUT_TIME("MCONTACT::CONTACT_ANALYSIS: coupling solver 1 done");
			std::vector<Eigen::VectorXd> outpDisp(multGrid.size());
			#pragma omp parallel for
			for(long tv = 0; tv < multGrid.size(); tv ++){
				Eigen::VectorXd resuSolu = globSolu.block(
					baseReco[tv],0,multGrid[tv].mgpi.consStif[doleMcsc[tv]].rows(),1
				);
				resuSolu = accuProl[tv] * resuSolu;
				multGrid[tv].OUTP_SUB1(resuSolu, outpDisp[tv]);
				resuDisp[tv] += outpDisp[tv];
				multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
			}
			// #pragma omp parallel for
			// for(long ts = 0; ts < searCont.size(); ts ++){
				// for(long tv = 0; tv < 2; tv ++){
					// Eigen::VectorXd inteForc = 
						// inteInte[ts][2 * tv] * outpDisp[contBody[ts][tv]] 
						// + inteInte[ts][2 * tv + 1] * outpDisp[contBody[ts][1 - tv]];
					// if(inteMass[ts][tv].rows() < DIRE_MAXI){
						// inteAuxi[ts][tv] += (inteDiso[ts][tv]).solve(inteForc);
					// }
					// else{
						// COGR_SOLV solv_0;
						// solv_0.compute(inteMass[ts][tv]);
						// inteAuxi[ts][tv] += solv_0.solve(inteForc);
					// }
				// }
			// }
		}
		//********************************interface balance**********************************
		std::cout << "Auxiliary solver";
		OUTPUT_TIME("");
		std::vector<Eigen::VectorXd> inpoGamm(searCont.size());
		#pragma omp parallel for
		for(long ts = 0; ts < searCont.size(); ts ++){
			//
			inpoGamm[ts] = 0.5 * (inpoLagr[ts][0] * inteLagr[ts][0] 
				- inpoLagr[ts][1] * inteLagr[ts][1] 
				+ pemaInpo_r[ts][0] * resuDisp[contBody[ts][0]] 
				- pemaInpo_r[ts][1] * resuDisp[contBody[ts][1]] 
				- pemaInpo[ts] * inpoNgap[ts]);
			if(fricCoef[ts] >= 0.0){
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					if(fricCoef[ts] == 0.0){
						inpoGamm[ts](ti) = std::max(0.0, inpoGamm[ts](ti));
					}
					else{
						inpoGamm[ts](3 * ti + 0) = std::max(0.0, inpoGamm[ts](3 * ti + 0));
					}
				}
			}
			Eigen::VectorXi fricStat = Eigen::VectorXi::Zero(inpoGamm[ts].rows());
			if(fricCoef[ts] > 0.0){
				for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
					if(inpoGamm[ts](3 * ti + 0) > 0.0){
						double tempSlid = fricCoef[ts] * inpoGamm[ts](3 * ti + 0);
						Eigen::Vector2d tempOrig = inpoGamm[ts].block(3 * ti + 1, 0, 2, 1);
						if(tempOrig.norm() >= tempSlid){
							inpoGamm[ts].block(3 * ti + 1, 0, 2, 1) = tempSlid 
								/ tempOrig.norm() * tempOrig;
							fricStat[3 * ti + 1] = 1;
						}
						else{
							fricStat[3 * ti + 1] = 2;
						}
					}
					else{
						inpoGamm[ts](3 * ti + 1) = 0.0;
						inpoGamm[ts](3 * ti + 2) = 0.0;
						fricStat[3 * ti + 1] = 0;
					}
				}
			}
			OUTPUT_PRTR(inpoGamm[ts], fricStat, ts);
			//
			for(long tv = 0; tv < 2; tv ++){
				Eigen::VectorXd inteForc = 
					systTran_pena[ts][tv].transpose() * resuDisp[contBody[ts][tv]] 
					+ inteMass[ts][tv] * inteLagr[ts][tv] 
					+ inteInpo[ts][tv] * inpoGamm[ts];
				if(inteMass_pena[ts][tv].rows() < DIRE_MAXI){
					inteAuxi[ts][tv] = (inteDiso_pena[ts][tv]).solve(inteForc);
				}
				else{
					COGR_SOLV solv_0;
					solv_0.compute(inteMass_pena[ts][tv]);
					inteAuxi[ts][tv] = solv_0.solve(inteForc);
				}
			}
		}
		//*******************************Lagrange multiplier*********************************
		std::cout << "Lagrange multiplier solver";
		OUTPUT_TIME("");
		#pragma omp parallel for
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(long tv = 0; tv < 2; tv ++){
				Eigen::VectorXd inteForc = 
					systTran_pena[ts][tv].transpose() * resuDisp[contBody[ts][tv]] 
					- inteMass_pena[ts][tv] * inteAuxi[ts][tv];
				if(inteMass[ts][tv].rows() < DIRE_MAXI){
					inteLagr[ts][tv] += (inteDiso[ts][tv]).solve(inteForc);
				}
				else{
					COGR_SOLV solv_0;
					solv_0.compute(inteMass[ts][tv]);
					inteLagr[ts][tv] += solv_0.solve(inteForc);
				}
			}
		}
		if(debuMode >= 1){
			OUTPUT_AULA();
		}
		//********************************stopping criterion*********************************
		if(MONITOR(tc, tempOfst, moniReco, resuDisp_0, inteAuxi_0, inteLagr_0) == 1){
			break;
		}
	}
	tempOfst.close();
	if(tc >= maxiIter){
		std::cout << "Nonconvergence in MCONTACT::CONTACT_ANALYSIS";
	}
	else{
		std::cout << "Converge after " << tc << "-th iterations";
	}
	OUTPUT_TIME("");
	return 1;
}

long MCONTACT::MONITOR(const long tc, std::ofstream &tempOfst, VECTOR2D &moniReco, 
	const std::vector<Eigen::VectorXd> &resuDisp_0, 
	const std::vector<std::vector<Eigen::VectorXd>> &inteAuxi_0, 
	const std::vector<std::vector<Eigen::VectorXd>> &inteLagr_0){
	long moniCycl = moniReco[0].size();
	bool tempFlag_0 = (tc >= moniCycl) ? true : false;
	bool tempFlag_1 = true;
	double critRati_0 = 0.1;
	double critRati_1 = 1.0E-12;
	double critRati_2 = 1.0E-10;
	double convValu = 0.0;
	double convCrit = 0.0;
	for(long tv = 0; tv < resuDisp.size(); tv ++){
		moniReco[tv][tc % moniCycl] = (resuDisp[tv] - resuDisp_0[tv]).squaredNorm();
		double dispAllo = resuDisp[tv].squaredNorm();
		convValu += moniReco[tv][tc % moniCycl];
		convCrit += dispAllo;
		tempOfst << std::setw(30) << moniReco[tv][tc % moniCycl] 
			<< std::setw(30) << dispAllo;
		if(debuMode >= 1){
			std::cout << "tv = " << tv << ": " << "convValu_D = " 
				<< moniReco[tv][tc % moniCycl];
			OUTPUT_TIME("");
		}
		if(tc >= moniCycl){
			double dispMedi, dispOsci;
			VECT_MEDI_OSCI(moniReco[tv], dispMedi, dispOsci);
			if(dispOsci > critRati_0 * dispMedi){
				if(debuMode >= 1){
					std::cout << "tv = " << tv << ": " << "dispMedi = " 
						<< dispMedi << ", dispOsci = " << dispOsci << std::endl;
				}
				tempFlag_0 = false;
			}
		}
		if(moniReco[tv][tc % moniCycl] > critRati_1 * dispAllo){
			if(debuMode >= 1){
				std::cout << "tv = " << tv << ": " 
					<< "dispAllo = " << dispAllo << std::endl;
			}
			tempFlag_1 = false;
		}
	}
	for(long ts = 0; ts < inteAuxi.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			//
			long tempIndi = resuDisp.size() + 4 * ts + 2 * tv + 0;
			moniReco[tempIndi][tc % moniCycl] = 
				(inteAuxi[ts][tv] - inteAuxi_0[ts][tv]).squaredNorm();
			double auxiAllo = inteAuxi[ts][tv].squaredNorm();
			convValu += moniReco[tempIndi][tc % moniCycl];
			convCrit += auxiAllo;
			tempOfst << std::setw(30) << moniReco[tempIndi][tc % moniCycl] 
				 << std::setw(30) << auxiAllo;
			if(debuMode >= 1){
				std::cout << "ts = " << ts << ",tv = " << tv << ": " 
					<< "convValu_A = " << moniReco[tempIndi][tc % moniCycl];
				OUTPUT_TIME("");
			}
			if(tc >= moniCycl){
				double auxiMedi, auxiOsci;
				VECT_MEDI_OSCI(moniReco[tempIndi], auxiMedi, auxiOsci);
				if(auxiOsci > critRati_0 * auxiMedi){
					if(debuMode >= 1){
						std::cout << "ts = " << ts << ",tv = " << tv << ": " << "auxiMedi = " 
							<< auxiMedi << ", auxiOsci = " << auxiOsci << std::endl;
					}
					tempFlag_0 = false;
				}
			}
			if(moniReco[tempIndi][tc % moniCycl] > critRati_1 * auxiAllo){
				if(debuMode >= 1){
					std::cout << "ts = " << ts << ",tv = " << tv << ": " 
						<< "auxiAllo = " << auxiAllo << std::endl;
				}
				tempFlag_1 = false;
			}
			//
			tempIndi = tempIndi + 1;
			moniReco[tempIndi][tc % moniCycl] = 
				(inteLagr[ts][tv] - inteLagr_0[ts][tv]).squaredNorm();
			double lagrAllo = inteLagr[ts][tv].squaredNorm();
			tempOfst << std::setw(30) << moniReco[tempIndi][tc % moniCycl] 
				<< std::setw(30) << lagrAllo;
			if(debuMode >= 1){
				std::cout << "ts = " << ts << ",tv = " << tv << ": " 
					<< "convValu_L = " << moniReco[tempIndi][tc % moniCycl];
				OUTPUT_TIME("");
			}
			if(tc >= moniCycl){
				double lagrMedi, lagrOsci;
				VECT_MEDI_OSCI(moniReco[tempIndi], lagrMedi, lagrOsci);
				if(lagrOsci > critRati_0 * lagrMedi){
					if(debuMode >= 1){
						std::cout << "ts = " << ts << ",tv = " << tv << ": " << "lagrMedi = " 
							<< lagrMedi << ", lagrOsci = " << lagrOsci << std::endl;
					}
					// tempFlag_0 = false;//more strict convergence criterion
				}
			}
			if(moniReco[tempIndi][tc % moniCycl] > critRati_2 * lagrAllo){
				if(debuMode >= 1){
					std::cout << "ts = " << ts << ",tv = " << tv << ": " 
						<< "lagrAllo = " << lagrAllo << std::endl;
				}
				// tempFlag_1 = false;
			}
		}
	}
	std::cout << "Cvalu = " << convValu << ", Ccrit = " << convCrit << std::endl;
	tempOfst << std::setw(30) << convValu << std::setw(30) << convCrit;
	tempOfst << std::endl;
	//
	if(tempFlag_0 == true){
		MULT_MAXI = tc;
	}
	if(tempFlag_1 == true){
		return 1;
	}
	return -1;
}

long MCONTACT::LAGRANGE(long precType){
	resuDisp.resize(multGrid.size());
	//***************************************************************************************
	#pragma omp parallel for
	for(long tv = 0; tv < multGrid.size(); tv ++){
		multGrid[tv].TRANSFER();
		multGrid[tv].STIF_MATR();
		if(precType == 1){
			multGrid[tv].CONSTRAINT(1);
		}
		else{
			multGrid[tv].CONSTRAINT(-1);
		}
	}
	long baseNumb = 0;
	baseReco.resize(multGrid.size() + 1, 0);
	for(long tv = 0; tv < multGrid.size(); tv ++){
		baseReco[tv] = baseNumb;
		const auto &tempStif = multGrid[tv].mgpi.consStif[multGrid[tv].mgpi.maxiLeve];
		baseNumb += tempStif.rows();
	}
	baseReco[multGrid.size()] = baseNumb;
	//***************************************************************************************
	//non-mortar node: can not be hanging node
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(auto iterInpo = searCont[ts].intePoin.begin(); 
			iterInpo != searCont[ts].intePoin.end();){
			bool tempFlag = true;
			long tv = 0;
			for(long tk = 0; tk < 4; tk ++){
				long tempNode = (iterInpo->node)[tv][tk];
				long tempMult = contBody[ts][tv];
				if(multGrid[tempMult].nodeLepo[tempNode][0] 
					== multGrid[tempMult].mgpi.maxiLeve + 1){
					tempFlag = false;
					break;
				}
			}
			if(tempFlag == false){
				iterInpo = searCont[ts].intePoin.erase(iterInpo);
			}
			else{
				iterInpo ++;
			}
		}
	}
	//dual mortar, boundary consistent treatment
	OUTPUT_TIME("Boundary consistent treatment");
	//segment: INTEGRAL_POINT
	std::vector<std::map<std::vector<long>,VECTOR2L>> segmInpo(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		long tv = 0;
		for(long ti = 0; ti < searCont[ts].intePoin.size(); ti ++){
			std::vector<long> tempSegm = searCont[ts].intePoin[ti].node[tv];
			// sort(tempSegm.begin(), tempSegm.end());//!!!!!NO NEED!!!!!
			std::vector<long> tempInpo = {ts, ti};
			auto iterSein = segmInpo[ts].find(tempSegm);
			if(iterSein == segmInpo[ts].end()){
				VECTOR2L tempVect = {tempInpo};
				segmInpo[ts].emplace(tempSegm, tempVect);
			}
			else{
				(iterSein->second).emplace_back(tempInpo);
			}
		}
	}
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		long tv = 0;
		for(const auto &iterSein : segmInpo[ts]){
			//
			Eigen::MatrixXd tempD = Eigen::MatrixXd::Zero(4,4);
			Eigen::MatrixXd tempM = Eigen::MatrixXd::Zero(4,4);
			for(const auto &iterScip : iterSein.second){//Scip: searCont, intePoin
				const INTEGRAL_POINT &tempInpo = searCont[iterScip[0]].intePoin[iterScip[1]];
				Eigen::MatrixXd tempD_tj = Eigen::MatrixXd::Zero(4,4);
				Eigen::MatrixXd tempM_tj = Eigen::MatrixXd::Zero(4,4);
				for(long tk = 0; tk < 4; tk ++){
					tempD_tj(tk,tk) = tempInpo.shapFunc[tv][tk];
					for(long tm = 0; tm < 4; tm ++){
						tempM_tj(tk,tm) = tempInpo.shapFunc[tv][tk] * tempInpo.shapFunc[tv][tm];
					}
				}
				tempD_tj = tempInpo.quadWeig * tempD_tj;
				tempM_tj = tempInpo.quadWeig * tempM_tj;
				tempD = tempD + tempD_tj;
				tempM = tempM + tempM_tj;
			}
			Eigen::MatrixXd tempA = tempD * tempM.inverse();
			//
			for(const auto &iterScip : iterSein.second){
				INTEGRAL_POINT &tempInpo = searCont[iterScip[0]].intePoin[iterScip[1]];
				Eigen::Vector<double,4> tempShfu;
				tempShfu << tempInpo.shapFunc[tv][0], tempInpo.shapFunc[tv][1], 
					tempInpo.shapFunc[tv][2], tempInpo.shapFunc[tv][3];
				tempInpo.dual = tempA * tempShfu;
			}
		}
	}
	segmInpo.clear();
	//***************************************************************************************
	//a non-mortar node: 0 - inactive, 1 - active/sliding, 2 - active/sticking
	std::vector<std::map<std::vector<long>,long>> acinStsl(searCont.size());
	//a non-mortar node: a number
	std::vector<std::map<std::vector<long>,long>> acinNumb(searCont.size());
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		long tv = 0;
		long tempStat = (fricCoef[ts] == 0.0) ? 1 : 2;
		for(const auto &iterInpo : searCont[ts].intePoin){
			for(long tk = 0; tk < 4; tk ++){
				std::vector<long> tempNode = {contBody[ts][tv], iterInpo.node[tv][tk]};
				//the same new key will be abandoned!!!!!
				acinStsl[ts].emplace(tempNode, tempStat);
				long tempSize = acinNumb[ts].size();
				acinNumb[ts].emplace(tempNode, tempSize);
			}
		}
	}
	//all non-mortar nodes: normal vector, tangential vector
	std::map<std::vector<long>,Eigen::Matrix3d> nmnoNota;
	std::map<std::vector<long>,double> notaWeig;
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(const auto &iterAcst : acinNumb[ts]){
			nmnoNota.emplace(iterAcst.first, Eigen::Matrix3d::Zero());
			notaWeig.emplace(iterAcst.first, 0.0);
		}
	}
	for(long ts = 0; ts < searCont.size(); ts ++){
		long tv = 0;
		for(const auto &iterInpo : searCont[ts].intePoin){
			for(long tk = 0; tk < 4; tk ++){
				std::vector<long> tempNode = {contBody[ts][tv], iterInpo.node[tv][tk]};
				auto iterNmno = nmnoNota.find(tempNode);
				(iterNmno->second).col(0) = (iterNmno->second).col(0).eval()
					+ iterInpo.quadWeig * iterInpo.basiVect[0];
				auto iterNowe = notaWeig.find(tempNode);
				(iterNowe->second) += iterInpo.quadWeig;
			}
		}
	}
	for(auto &iterNmno : nmnoNota){
		auto iterNowe = notaWeig.find(iterNmno.first);
		(iterNmno.second).col(0) = (iterNmno.second).col(0) / (iterNowe->second);
		double d1x = (iterNmno.second)(0,0);
		double d1y = (iterNmno.second)(1,0);
		double d1z = (iterNmno.second)(2,0);
		double EPSI = 1.0E-14;
		if(abs(d1y) < EPSI){
			(iterNmno.second)(1,1) = 1.0;
			if(abs(d1x) < EPSI){
				(iterNmno.second)(0,2) = 1.0;
			}
			else{
				if(abs(d1z) < EPSI){
					(iterNmno.second)(2,2) = 1.0;
				}
				else{
					(iterNmno.second)(0,2) = d1z / sqrt(pow(d1z, 2.0) + pow(d1x, 2.0));
					(iterNmno.second)(2,2) = - d1x / sqrt(pow(d1z, 2.0) + pow(d1x, 2.0));
				}
			}
		}
		else if(abs(d1z) < EPSI){
			(iterNmno.second)(2,1) = 1.0;
			if(abs(d1x) < EPSI){
				(iterNmno.second)(0,2) = 1.0;
			}
			else{
				(iterNmno.second)(0,2) = d1y / sqrt(pow(d1y, 2.0) + pow(d1x, 2.0));
				(iterNmno.second)(1,2) = - d1x / sqrt(pow(d1y, 2.0) + pow(d1x, 2.0));
			}
		}
		else{
			double b2c2 = pow(d1y, 2.0) + pow(d1z, 2.0);
			(iterNmno.second)(1,1) = d1z / sqrt(b2c2);
			(iterNmno.second)(2,1) = - d1y / sqrt(b2c2);
			double abc = sqrt(pow(b2c2, 2.0) + pow(d1x * d1y, 2.0) + pow(d1x * d1z, 2.0));
			(iterNmno.second)(0,2) = b2c2 / abc;
			(iterNmno.second)(1,2) = (- d1x * d1y) / abc;
			(iterNmno.second)(2,2) = (- d1x * d1z) / abc;
		}
		Eigen::Vector3d tempVect_0 = (iterNmno.second).col(0);
		Eigen::Vector3d tempVect_1 = (iterNmno.second).col(1);
		Eigen::Vector3d tempVect_2 = (iterNmno.second).col(2);
		if((tempVect_0.cross(tempVect_1)).dot(tempVect_2) < 0.0){
			(iterNmno.second).col(2) = - tempVect_2;
		}
	}
	notaWeig.clear();
	//***************************************************************************************
	OUTPUT_TIME("Mortar coupling matrix");
	SPARMATS notaMoco;//normal tangential mortar coupling
	VECT_RESI(notaMoco, searCont.size(), 2);
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			std::vector<Eigen::Triplet<double>> mocoList;
			for(const auto &iterInpo : searCont[ts].intePoin){
				Eigen::VectorXd U_e = iterInpo.dual;
				Eigen::Matrix<double,3,12> D_e;
				D_e << 
					U_e(0), 0.0, 0.0, U_e(1), 0.0, 0.0, U_e(2), 0.0, 0.0, U_e(3), 0.0, 0.0, 
					0.0, U_e(0), 0.0, 0.0, U_e(1), 0.0, 0.0, U_e(2), 0.0, 0.0, U_e(3), 0.0, 
					0.0, 0.0, U_e(0), 0.0, 0.0, U_e(1), 0.0, 0.0, U_e(2), 0.0, 0.0, U_e(3);
				std::vector<double> M_e = iterInpo.shapFunc[tv];
				Eigen::Matrix<double,3,12> N_e;
				N_e << 
					M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 0.0, 
					0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3], 0.0, 
					0.0, 0.0, M_e[0], 0.0, 0.0, M_e[1], 0.0, 0.0, M_e[2], 0.0, 0.0, M_e[3];
				Eigen::MatrixXd matr_0 = iterInpo.quadWeig * D_e.transpose()* N_e;
				for(long tj = 0; tj < 4; tj ++){
					std::vector<long> tempNode = {contBody[ts][0], iterInpo.node[0][tj]};
					auto iterAinu = acinNumb[ts].find(tempNode);
					for(long tk = 0; tk < 3; tk ++){
						long row_tj = 3 * (iterAinu -> second) + tk;
						for(long tm = 0; tm < 4; tm ++){
							long node_tm = iterInpo.node[tv][tm];
							for(long tn = 0; tn < 3; tn ++){
								long col_tn = 3 * node_tm + tn;
								if(tv == 0 && 3 * tj + tk != 3 * tm + tn){
									continue;
								}
								mocoList.emplace_back(
									row_tj, col_tn, matr_0(3 * tj + tk, 3 * tm + tn)
								);
							}
						}
					}
				}
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> mocoMatr;
			const auto &tempGrid = multGrid[contBody[ts][tv]];
			long tempNumb = 3 * tempGrid.nodeCoor.size();
			mocoMatr.resize(3 * acinNumb[ts].size(), tempNumb);
			mocoMatr.setFromTriplets(mocoList.begin(), mocoList.end());
			mocoMatr = mocoMatr * tempGrid.earlTran;
			mocoMatr = mocoMatr * tempGrid.prolOper[tempGrid.mgpi.maxiLeve];
			mocoMatr = mocoMatr * tempGrid.consOper[tempGrid.mgpi.maxiLeve].transpose();
			if(tv == 1){
				mocoMatr = - mocoMatr;
			}
			std::vector<Eigen::Triplet<double>> notaList;
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterNmno = nmnoNota.find(iterAinu.first);
				for(long tj = 0; tj < 3; tj ++){
					for(long tk = 0; tk < 3; tk ++){
						notaList.emplace_back(
							3 * iterAinu.second + tj, 3 * iterAinu.second + tk, 
							(iterNmno->second)(tk,tj)//!!!Order!!!
						);
					}
				}
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> notaMatr;
			notaMatr.resize(3 * acinNumb[ts].size(), 3 * acinNumb[ts].size());
			notaMatr.setFromTriplets(notaList.begin(), notaList.end());
			notaMoco[ts][tv] = notaMatr * mocoMatr;
		}
	}
	//***************************************************************************************
	std::vector<Eigen::VectorXd> initGapd(searCont.size());//initial gap distance
	#pragma omp parallel for
	for(long ts = 0; ts < searCont.size(); ts ++){
		initGapd[ts] = Eigen::VectorXd::Zero(3 * acinNumb[ts].size());
		for(const auto &iterInpo : searCont[ts].intePoin){
			Eigen::VectorXd U_e = iterInpo.dual;
			Eigen::MatrixXd matr_0 = iterInpo.quadWeig * U_e.transpose()* iterInpo.initNgap;
			for(long tj = 0; tj < 4; tj ++){
				std::vector<long> tempNode = {contBody[ts][0], iterInpo.node[0][tj]};
				auto iterAinu = acinNumb[ts].find(tempNode);
				initGapd[ts](3 * (iterAinu -> second) + 0) += matr_0(tj);
			}
		}
	}
	//***************************************************************************************
	OUTPUT_TIME("Original coupling matrix");
	std::vector<Eigen::Triplet<double>> origList;
	for(long tv = 0; tv < multGrid.size(); tv ++){
		const auto &tempGrid = multGrid[tv];
		const auto &tempStif = tempGrid.mgpi.consStif[tempGrid.mgpi.maxiLeve];
		for(long ti = 0; ti < tempStif.rows(); ti ++){
			for(RSPA_INNE iterStif(tempStif, ti); iterStif; ++ iterStif){
				origList.emplace_back(
					baseReco[tv] + ti, baseReco[tv] + iterStif.col(), iterStif.value()
				);
			}
		}
	}
	long contNumb = 0;
	std::vector<long> acinReco(searCont.size() + 1, 0);
	for(long ts = 0; ts < searCont.size(); ts ++){
		acinReco[ts] = contNumb;
		contNumb += acinNumb[ts].size();
	}
	acinReco[searCont.size()] = contNumb;
	for(long ts = 0; ts < searCont.size(); ts ++){
		for(long tv = 0; tv < 2; tv ++){
			const auto &tempMatr = notaMoco[ts][tv];
			for(long ti = 0; ti < tempMatr.rows(); ti ++){
				for(RSPA_INNE iterMatr(tempMatr, ti); iterMatr; ++ iterMatr){
					origList.emplace_back(
						baseReco[multGrid.size()] + 3 * acinReco[ts] + ti, 
						baseReco[contBody[ts][tv]] + iterMatr.col(), 
						iterMatr.value()
					);
					origList.emplace_back(
						baseReco[contBody[ts][tv]] + iterMatr.col(), 
						baseReco[multGrid.size()] + 3 * acinReco[ts] + ti, 
						iterMatr.value()
					);
				}
			}
		}
	}
	long origNumb = baseReco[multGrid.size()] + 3 * acinReco[searCont.size()];
	// Eigen::SparseMatrix<double,Eigen::RowMajor> origCoup;
	// origCoup.resize(origNumb, origNumb);
	// origCoup.setFromTriplets(origList.begin(), origList.end());
	Eigen::VectorXd origForc = Eigen::VectorXd::Zero(origNumb);
	for(long tv = 0; tv < multGrid.size(); tv ++){
		const auto &tempGrid = multGrid[tv];
		origForc.block(baseReco[tv],0,tempGrid.consForc.rows(),1) = tempGrid.consForc;
	}
	for(long ts = 0; ts < searCont.size(); ts ++){
		origForc.block(
			baseReco[multGrid.size()] + 3 * acinReco[ts],0,3 * acinNumb[ts].size(),1
		) = initGapd[ts];
	}
	//***************************************************************************************
	//relative displacement of non-mortar side relative to mortar side, no need to initialize!
	std::vector<Eigen::VectorXd> nmnoWedi(searCont.size());
	std::vector<std::map<std::vector<long>,long>> aissHist = acinStsl;
	std::vector<Eigen::VectorXd> resuLagr(searCont.size());
	for(long tc = 0; true; tc ++){
		//***************************************************************************************
		OUTPUT_TIME("SlidCoup generation");
		std::vector<Eigen::Triplet<double>> slidList = origList;
		for(long ts = 0; ts < searCont.size(); ts ++){
			if(fricCoef[ts] <= 0.0){
				continue;
			}
			std::vector<Eigen::Triplet<double>> tempList;
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				//for tc = 0, fricCoef[ts] > 0.0 AND iterAcst->second == 1 
				//will never be simultaneously satisfied
				if(iterAcst->second == 1){
					auto iterAihi = aissHist[ts].find(iterAinu.first);
					double tangComp_0, tangComp_1;
					if(iterAihi->second == 1 || iterAihi->second == 0){
						tangComp_0 = (nmnoWedi[ts])(3 * iterAinu.second + 1);
						tangComp_1 = (nmnoWedi[ts])(3 * iterAinu.second + 2);
					}
					else{//iterAihi->second == 2
						//- * - = +, the left side of work equation -> move to the right side
						tangComp_0 = (resuLagr[ts])(3 * iterAinu.second + 1);
						tangComp_1 = (resuLagr[ts])(3 * iterAinu.second + 2);
					}
					double tangTota = sqrt(pow(tangComp_0, 2.0) + pow(tangComp_1, 2.0));
					tangComp_0 = tangComp_0 / tangTota;
					tangComp_1 = tangComp_1 / tangTota;
					tempList.emplace_back(
						3 * iterAinu.second + 0, 3 * iterAinu.second + 1, tangComp_0
					);
					tempList.emplace_back(
						3 * iterAinu.second + 0, 3 * iterAinu.second + 2, tangComp_1
					);
				}
			}
			Eigen::SparseMatrix<double,Eigen::RowMajor> slidMatr;
			slidMatr.resize(3 * acinNumb[ts].size(), 3 * acinNumb[ts].size());
			slidMatr.setFromTriplets(tempList.begin(), tempList.end());
			for(long tv = 0; tv < 2; tv ++){
				Eigen::SparseMatrix<double,Eigen::RowMajor> tempMatr = 
					fricCoef[ts] * slidMatr * notaMoco[ts][tv];
				for(long ti = 0; ti < tempMatr.rows(); ti ++){
					for(RSPA_INNE iterMatr(tempMatr, ti); iterMatr; ++ iterMatr){
						slidList.emplace_back(
							baseReco[contBody[ts][tv]] + iterMatr.col(), 
							baseReco[multGrid.size()] + 3 * acinReco[ts] + ti, 
							iterMatr.value()
						);
					}
				}
			}
		}
		Eigen::SparseMatrix<double,Eigen::RowMajor> slidCoup;
		slidCoup.resize(origNumb, origNumb);
		slidCoup.setFromTriplets(slidList.begin(), slidList.end());
		slidList.clear();
		//***************************************************************************************
		OUTPUT_TIME("RealCoup generation");
		std::vector<Eigen::Triplet<double>> realList;
		for(long ti = 0; ti < baseReco[multGrid.size()]; ti ++){
			realList.emplace_back(ti, ti, 1.0);
		}
		long consNumb = 0;
		for(long ts = 0; ts < searCont.size(); ts ++){
			long tempBase = baseReco[multGrid.size()] + 3 * acinReco[ts];
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				if(iterAcst->second == 1){
					realList.emplace_back(
						baseReco[multGrid.size()] + consNumb, 
						tempBase + 3 * iterAinu.second + 0, 
						1.0
					);
					consNumb ++;
				}
				else if(iterAcst->second == 2){
					for(long td = 0; td < 3; td ++){
						realList.emplace_back(
							baseReco[multGrid.size()] + consNumb, 
							tempBase + 3 * iterAinu.second + td, 
							1.0
						);
						consNumb ++;
					}
				}
			}
		}
		long totaNumb = baseReco[multGrid.size()] + consNumb;//total number of DOF
		Eigen::SparseMatrix<double,Eigen::RowMajor> realMatr;
		realMatr.resize(totaNumb, origNumb);
		realMatr.setFromTriplets(realList.begin(), realList.end());
		realList.clear();
		Eigen::SparseMatrix<double,Eigen::RowMajor> realCoup = 
			realMatr * slidCoup * realMatr.transpose();
		Eigen::VectorXd realForc = realMatr * origForc;
		//***************************************************************************************
		OUTPUT_TIME("SoluCoup generation");
		//condensed DOF numbering
		std::vector<std::vector<std::set<long>>> condDofn(searCont.size());
		consNumb = 0;
		for(long ts = 0; ts < searCont.size(); ts ++){
			condDofn[ts].resize(acinNumb[ts].size());
			long starColn = baseReco[contBody[ts][0]];
			long termColn = baseReco[contBody[ts][0] + 1];
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				if(iterAcst->second == 1){
					long row_temp = baseReco[multGrid.size()] + consNumb;
					double tempMaxi = -1.0;
					long tempCond = -1;
					for(RSPA_INNE iterMatr(realCoup, row_temp); iterMatr; ++ iterMatr){
						if(starColn <= iterMatr.col() && iterMatr.col() < termColn 
							&& abs(iterMatr.value()) > tempMaxi){
							tempMaxi = abs(iterMatr.value());
							tempCond = iterMatr.col();
						}
					}
					condDofn[ts][iterAinu.second].emplace(tempCond);
					consNumb ++;
					if(tempCond == -1){
						std::cout << "ERROR 1 in MCONTACT::LAGRANGE!" << std::endl;
					}
				}
				else if(iterAcst->second == 2){
					for(long td = 0; td < 3; td ++){
						long row_temp = baseReco[multGrid.size()] + consNumb;
						for(RSPA_INNE iterMatr(realCoup, row_temp); iterMatr; ++ iterMatr){
							if(starColn <= iterMatr.col() && iterMatr.col() < termColn 
								/*&& abs(iterMatr.value()) > 1.0E-12*/){
								condDofn[ts][iterAinu.second].emplace(iterMatr.col());
							}
						}
						consNumb ++;
					}
					if(condDofn[ts][iterAinu.second].size() != 3){
						std::cout << "ERROR 2 in MCONTACT::LAGRANGE!" << std::endl;
					}
				}
			}
		}
		std::vector<Eigen::Triplet<double>> soluList;
		std::vector<bool> condFlag(totaNumb, false);
		consNumb = 0;
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(const auto &iterAinu : acinNumb[ts]){
				for(const auto &iterCodo : condDofn[ts][iterAinu.second]){
					soluList.emplace_back(iterCodo, consNumb, 1.0);
					condFlag[iterCodo] = true;
					consNumb ++;
				}
			}
		}
		for(long ti = 0; ti < totaNumb; ti ++){
			if(condFlag[ti] == false){
				soluList.emplace_back(ti, consNumb, 1.0);
				consNumb ++;
			}
		}
		Eigen::SparseMatrix<double,Eigen::RowMajor> soluMatr;
		soluMatr.resize(totaNumb, totaNumb);
		soluMatr.setFromTriplets(soluList.begin(), soluList.end());
		soluList.clear();
		Eigen::SparseMatrix<double,Eigen::RowMajor> soluCoup = 
			soluMatr.transpose() * realCoup * soluMatr;
		Eigen::VectorXd soluForc = soluMatr.transpose() * realForc;
		//***************************************************************************************
		consNumb = totaNumb - baseReco[multGrid.size()];
		Eigen::SparseMatrix<double,Eigen::RowMajor> K_00 = 
			soluCoup.block(0,0,consNumb,consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> K_01 = 
			soluCoup.block(0,consNumb,consNumb,totaNumb - 2 * consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> K_10 = 
			soluCoup.block(consNumb,0,totaNumb - 2 * consNumb,consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> K_11 = 
			soluCoup.block(consNumb,consNumb,totaNumb - 2 * consNumb,totaNumb - 2 * consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> T_0 = 
			soluCoup.block(baseReco[multGrid.size()],0,consNumb,consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> T_1 = 
			soluCoup.block(baseReco[multGrid.size()],consNumb,consNumb,totaNumb - 2 * consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> T_0f = 
			soluCoup.block(0,baseReco[multGrid.size()],consNumb,consNumb);
		Eigen::SparseMatrix<double,Eigen::RowMajor> T_1f = 
			soluCoup.block(consNumb,baseReco[multGrid.size()],totaNumb - 2 * consNumb,consNumb);
		Eigen::VectorXd F_0 = soluForc.block(0,0,consNumb,1);
		Eigen::VectorXd F_1 = soluForc.block(consNumb,0,baseReco[multGrid.size()] - consNumb,1);
		Eigen::VectorXd g_0 = soluForc.block(baseReco[multGrid.size()],0,consNumb,1);
		//
		std::vector<Eigen::Triplet<double>> inveList_0, inveList_1;
		consNumb = 0;
		for(long ts = 0; ts < searCont.size(); ts ++){
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				if(iterAcst->second == 1){
					inveList_0.emplace_back(
						consNumb, consNumb, 1.0 / T_0.coeff(consNumb, consNumb)
					);
					inveList_1.emplace_back(
						consNumb, consNumb, 1.0 / T_0f.coeff(consNumb, consNumb)
					);
					consNumb ++;
				}
				else if(iterAcst->second == 2){
					Eigen::Matrix3d tempMatr_0 = T_0.block(consNumb, consNumb, 3, 3);
					Eigen::Matrix3d tempInve_0 = tempMatr_0.inverse();
					Eigen::Matrix3d tempMatr_1 = T_0f.block(consNumb, consNumb, 3, 3);
					Eigen::Matrix3d tempInve_1 = tempMatr_1.inverse();
					for(long ti = 0; ti < 3; ti ++){
						for(long tj = 0; tj < 3; tj ++){
							inveList_0.emplace_back(
								consNumb + ti, consNumb + tj, tempInve_0(ti, tj)
							);
							inveList_1.emplace_back(
								consNumb + ti, consNumb + tj, tempInve_1(ti, tj)
							);
						}
					}
					consNumb += 3;
				}
			}
		}
		Eigen::SparseMatrix<double,Eigen::RowMajor> inveT_0, inveT_0f;
		inveT_0.resize(consNumb, consNumb);
		inveT_0.setFromTriplets(inveList_0.begin(), inveList_0.end());
		inveList_0.clear();
		inveT_0f.resize(consNumb, consNumb);
		inveT_0f.setFromTriplets(inveList_1.begin(), inveList_1.end());
		inveList_1.clear();
		//
		Eigen::SparseMatrix<double,Eigen::RowMajor> K = K_11 - K_10 * inveT_0 * T_1 
			- T_1f * inveT_0f * K_01 + T_1f * inveT_0f * K_00 * inveT_0 * T_1;
		Eigen::VectorXd F = F_1 - K_10 * inveT_0 * g_0 
			- T_1f * inveT_0f * F_0 + T_1f * inveT_0f * K_00 * inveT_0 * g_0;
		Eigen::VectorXd U_1;
		if(precType == 1){
			OUTPUT_TIME("Prolongation operator");
			MGPIS mgpi;
			mgpi.maxiLeve = multGrid[0].mgpi.maxiLeve;
			mgpi.realProl.resize(mgpi.maxiLeve);
			mgpi.consStif.resize(mgpi.maxiLeve + 1);
			mgpi.consStif[mgpi.maxiLeve] = K;
			//
			std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> origProl;
			origProl.resize(mgpi.maxiLeve);
			#pragma omp parallel for
			for(long tl = mgpi.maxiLeve - 1; tl >= 0; tl --){
				std::vector<Eigen::Triplet<double>> tempList;
				long rowsSumm = 0;
				long coluSumm = 0;
				for(long tv = 0; tv < multGrid.size(); tv ++){
					const Eigen::SparseMatrix<double,Eigen::RowMajor> &tempRepr = 
						multGrid[tv].mgpi.realProl[tl];
					for(long ti = 0; ti < tempRepr.rows(); ti ++){
						for(RSPA_INNE iterMatr(tempRepr, ti); iterMatr; ++ iterMatr){
							tempList.emplace_back(
								rowsSumm + ti, coluSumm + iterMatr.col(), iterMatr.value()
							);
						}
					}
					rowsSumm += tempRepr.rows();
					coluSumm += tempRepr.cols();
				}
				origProl[tl].resize(rowsSumm, coluSumm);
				origProl[tl].setFromTriplets(tempList.begin(), tempList.end());
			}
			//
			Eigen::SparseMatrix<double,Eigen::RowMajor> reseMaxi = 
				soluMatr.transpose().block(
					consNumb,0,baseReco[multGrid.size()] - consNumb,baseReco[multGrid.size()]
				);
			std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> condProl;
			condProl.resize(mgpi.maxiLeve + 1);
			condProl[mgpi.maxiLeve] = 
				soluMatr.block(0,0,baseReco[multGrid.size()],consNumb) 
				* (- inveT_0 * T_1 * reseMaxi);
			for(long tl = mgpi.maxiLeve - 1; tl >= 0; tl --){
				condProl[tl] = condProl[tl + 1] * origProl[tl];
			}
			std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> reseProl;
			reseProl.resize(mgpi.maxiLeve);
			#pragma omp parallel for
			for(long tl = mgpi.maxiLeve - 1; tl >= 0; tl --){
				//condensed dof before/after
				long coluSumm = 0;
				long tempNumb_c = 0;
				std::vector<long> codoBeaf_c(origProl[tl].cols(), -1);
				for(long tv = 0; tv < multGrid.size(); tv ++){
					const Eigen::SparseMatrix<double,Eigen::RowMajor> &tempRepr = 
						multGrid[tv].mgpi.realProl[tl];
					for(long ti = 0; ti < tempRepr.cols(); ti ++){
						long ti_real = baseReco[tv] + ti;
						if(condFlag[ti_real] == false){
							codoBeaf_c[coluSumm + ti] = tempNumb_c;
							tempNumb_c ++;
						}
					}
					coluSumm += tempRepr.cols();
				}
				//
				std::vector<Eigen::Triplet<double>> tempList;
				coluSumm = 0;
				for(long tv = 0; tv < multGrid.size(); tv ++){
					const Eigen::SparseMatrix<double,Eigen::RowMajor> &tempRepr = 
						multGrid[tv].mgpi.realProl[tl];
					for(long ti = 0; ti < tempRepr.cols(); ti ++){
						long ti_real = baseReco[tv] + ti;
						if(condFlag[ti_real] == false){
							tempList.emplace_back(coluSumm + ti, coluSumm + ti, 1.0);
						}
						else{
							for(RSPA_INNE iterMatr(condProl[tl], ti_real); 
								iterMatr; ++ iterMatr){
								tempList.emplace_back(
									coluSumm + ti, iterMatr.col(), iterMatr.value()
								);
								if(codoBeaf_c[iterMatr.col()] == -1){
									std::cout << "ERROR 3 in MCONTACT::LAGRANGE!" << std::endl;
								}
							}
						}
					}
					coluSumm += tempRepr.cols();
				}
				reseProl[tl].resize(coluSumm, coluSumm);
				reseProl[tl].setFromTriplets(tempList.begin(), tempList.end());
			}
			//
			#pragma omp parallel for
			for(long tl = mgpi.maxiLeve - 1; tl >= 0; tl --){
				//condensed dof before/after
				long rowsSumm = 0;
				long coluSumm = 0;
				long tempNumb_r = 0;
				long tempNumb_c = 0;
				std::vector<long> codoBeaf_r(origProl[tl].rows(), -1);
				std::vector<long> codoBeaf_c(origProl[tl].cols(), -1);
				for(long tv = 0; tv < multGrid.size(); tv ++){
					const Eigen::SparseMatrix<double,Eigen::RowMajor> &tempRepr = 
						multGrid[tv].mgpi.realProl[tl];
					for(long ti = 0; ti < tempRepr.rows(); ti ++){
						long ti_real = baseReco[tv] + ti;
						if(condFlag[ti_real] == false){
							codoBeaf_r[rowsSumm + ti] = tempNumb_r;
							tempNumb_r ++;
							if(ti < tempRepr.cols()){
								codoBeaf_c[coluSumm + ti] = tempNumb_c;
								tempNumb_c ++;
							}
						}
					}
					rowsSumm += tempRepr.rows();
					coluSumm += tempRepr.cols();
				}
				//
				Eigen::SparseMatrix<double,Eigen::RowMajor> realProl_0 = 
					origProl[tl] * reseProl[tl];
				std::vector<Eigen::Triplet<double>> tempList;
				for(long ti = 0; ti < realProl_0.rows(); ti ++){
					if(codoBeaf_r[ti] == -1){
						continue;
					}
					for(RSPA_INNE iterMatr(realProl_0, ti); iterMatr; ++ iterMatr){
						if(codoBeaf_c[iterMatr.col()] == -1){
							continue;
						}
						tempList.emplace_back(
							codoBeaf_r[ti], codoBeaf_c[iterMatr.col()], iterMatr.value()
						);
					}
				}
				mgpi.realProl[tl].resize(tempNumb_r, tempNumb_c);
				mgpi.realProl[tl].setFromTriplets(tempList.begin(), tempList.end());
			}
			for(long tl = mgpi.maxiLeve - 1; tl >= 0; tl --){
				mgpi.consStif[tl] = mgpi.realProl[tl].transpose() 
					* mgpi.consStif[tl + 1] * mgpi.realProl[tl];
			}
			mgpi.ESTABLISH();
			mgpi.BiCGSTAB_SOLV(1, F, U_1);
			// mgpi.GMRES_SOLV(1, F, U_1);
		}
		else if(precType == 2){//BiCGSTAB in Eigen
			Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>> solv;
			solv.setTolerance(1.0E-14);
			solv.setMaxIterations(K.rows());
			solv.compute(K);
			OUTPUT_TIME("Eigen::BiCGSTAB solver");
			U_1 = solv.solve(F);
			if(solv.info() != Eigen::Success){
				std::cout << "Solution failed!" << std::endl;
			}
			std::cout << "#Iterations:      " << solv.iterations() << std::endl;
			std::cout << "#Estimated error: " << solv.error();
			OUTPUT_TIME("");
		}
		//
		Eigen::VectorXd U_0 = inveT_0 * g_0 - inveT_0 * T_1 * U_1;
		Eigen::VectorXd lamb = inveT_0f * F_0 - inveT_0f * K_00 * inveT_0 * g_0 
			- inveT_0f * K_01 * U_1 + inveT_0f * K_00 * inveT_0 * T_1 * U_1;
		//***************************************************************************************
		Eigen::VectorXd soluDisp = Eigen::VectorXd::Zero(totaNumb);
		soluDisp.block(0,0,consNumb,1) = U_0;
		soluDisp.block(consNumb,0,baseReco[multGrid.size()] - consNumb,1) = U_1;
		soluDisp.block(baseReco[multGrid.size()],0,consNumb,1) = lamb;
		Eigen::VectorXd realDisp = soluMatr * soluDisp;
		Eigen::VectorXd slidDisp = realMatr.transpose() * realDisp;
		for(long ts = 0; ts < searCont.size(); ts ++){
			nmnoWedi[ts] = - initGapd[ts];
			for(long tv = 0; tv < 2; tv ++){
				nmnoWedi[ts] += notaMoco[ts][tv] * slidDisp.block(
					baseReco[contBody[ts][tv]],0,
					baseReco[contBody[ts][tv] + 1] - baseReco[contBody[ts][tv]],1
				);
			}
			resuLagr[ts] = slidDisp.block(
				baseReco[multGrid.size()] + 3 * acinReco[ts],0,
				3 * acinReco[ts + 1] - 3 * acinReco[ts],1
			);
			//no need to obtain resuLagr of sliding node
		}
		#pragma omp parallel for
		for(long tv = 0; tv < multGrid.size(); tv ++){
			Eigen::VectorXd resuSolu = slidDisp.block(
				baseReco[tv],0,baseReco[tv + 1] - baseReco[tv],1
			);
			multGrid[tv].OUTP_SUB1(resuSolu, resuDisp[tv]);
			multGrid[tv].OUTP_SUB2(resuDisp[tv], tv);
		}
		#pragma omp parallel for
		for(long ts = 0; ts < searCont.size(); ts ++){
			std::stringstream tempStre;
			tempStre << "resuLagr_" << ts << ".txt";
			std::ofstream tempOfst(DIRECTORY(tempStre.str()), std::ios::out);
			tempStre.str("");
			tempStre.clear();
			tempOfst << std::setiosflags(std::ios::scientific) << std::setprecision(20);
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				tempOfst << std::setw(10) << (iterAinu.first)[1] 
					<< std::setw(10) << iterAcst->second 
					<< std::setw(30) << resuLagr[ts](3 * iterAinu.second + 0);
				if(iterAcst->second != 1){
					tempOfst << std::setw(30) << resuLagr[ts](3 * iterAinu.second + 1) 
						<< std::setw(30) << resuLagr[ts](3 * iterAinu.second + 2) << std::endl;
				}
				else{
					tempOfst 
						<< std::setw(30) << fricCoef[ts] * resuLagr[ts](3 * iterAinu.second + 0) 
						<< std::setw(30) << 0.0 << std::endl;
				}
			}
			tempOfst.close();
		}
		aissHist = acinStsl;
		long seneNumb = 0;
		for(long ts = 0; ts < searCont.size(); ts ++){
			if(fricCoef[ts] < 0.0){
				continue;
			}
			for(const auto &iterAinu : acinNumb[ts]){
				auto iterAcst = acinStsl[ts].find(iterAinu.first);
				//semi-smooth Newton, normal
				double seneNorm = resuLagr[ts](3 * iterAinu.second + 0) 
					+ 210.0E9 * nmnoWedi[ts](3 * iterAinu.second + 0);
				if(seneNorm <= 0.0){
					if(iterAcst->second != 0){
						seneNumb ++;
					}
					iterAcst->second = 0;
					continue;
				}
				if(fricCoef[ts] == 0.0){
					if(iterAcst->second != 1){
						seneNumb ++;
					}
					iterAcst->second = 1;
					continue;
				}
				//semi-smooth Newton, tangential
				double seneTang;
				if(iterAcst->second == 2){
					//for sticking, relative displacement is always 0;
					seneTang = sqrt(pow(resuLagr[ts](3 * iterAinu.second + 1), 2.0) 
						+ pow(resuLagr[ts](3 * iterAinu.second + 2), 2.0));
				}
				else{//iterAcst->second == 0/1
					//for sliding: the relative sliding velocity direction 
					//is the same (or only changes slightly) in two consecutive iterations
					seneTang = fricCoef[ts] * resuLagr[ts](3 * iterAinu.second + 0) 
						+ 210.0E9 * sqrt(pow(nmnoWedi[ts](3 * iterAinu.second + 1), 2.0) 
						+ pow(nmnoWedi[ts](3 * iterAinu.second + 2), 2.0));
				}
				if(seneTang >= fricCoef[ts] * seneNorm){
					if(iterAcst->second != 1){
						seneNumb ++;
					}
					iterAcst->second = 1;
				}
				else{
					if(iterAcst->second != 2){
						seneNumb ++;
					}
					iterAcst->second = 2;
				}
			}
		}
		if(seneNumb == 0){
			std::cout << "Converge after " << tc << "-th iteration";
			OUTPUT_TIME("");
			break;
		}
		else{
			std::cout << "There are " << seneNumb << " unconverged constraints at the " 
				<< tc << "-th iteration." << std::endl;
		}
	}
	return 1;
}

#endif
