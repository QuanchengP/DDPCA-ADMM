
#ifndef _MGPIS_H
#define _MGPIS_H

#include "PREP.h"

//multigrid preconditioned iterative solver
class MGPIS{
public:
	/*********************************************************************************************/
	//maximum level
	long maxiLeve;
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> realProl;//prolongation operator
	//consStif = realProl^T * origStif * realProl
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> consStif;
	long ESTABLISH();
	//multigrid solver
	long MULT_SOLV(const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu);
	//precSwit: preconditioner switch
	//0 - diagonally preconditioned conjugate gradient solver
	//1 - multigrid preconditioned conjugate gradient solver
	long CG_SOLV(long precSwit, const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu);
	//multigrid preconditioned GMRES
	long GMRES_SOLV(long precSwit, const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu);
	//multigrid preconditioned BiCGSTAB
	long BiCGSTAB_SOLV(long precSwit, const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu);
	/*********************************************************************************************/
	//lower triangular matrix with zeros on the diagonal
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> consLowe;
	//diagonal matrix
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> consDiag;
	//upper triangular matrix with zeros on the diagonal
	std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> consUppe;
	//V cycle
	long MULT_VCYC(long tempLeve, const Eigen::VectorXd &righHand, 
		Eigen::VectorXd &resuSolu, const DIRE_SOLV &direSolv
	);
};

long MGPIS::ESTABLISH(){
	consLowe.resize(maxiLeve + 1);
	consDiag.resize(maxiLeve + 1);
	consUppe.resize(maxiLeve + 1);
	for(long ti = 0; ti <= maxiLeve; ti ++){
		consLowe[ti] = consStif[ti].triangularView<Eigen::StrictlyLower>();
		consDiag[ti].resize(consStif[ti].rows(), consStif[ti].cols());
		for(long tj = 0; tj < consStif[ti].rows(); tj ++){
			consDiag[ti].insert(tj,tj) = consStif[ti].coeff(tj,tj);
		}
		consUppe[ti] = consStif[ti].triangularView<Eigen::StrictlyUpper>();
	}
	return 1;
}

long MGPIS::MULT_VCYC(long tempLeve, const Eigen::VectorXd &righHand, 
	Eigen::VectorXd &resuSolu, const DIRE_SOLV &direSolv){
	if(tempLeve == 0){
		resuSolu = direSolv.solve(righHand);
		return 1;
	}
	//Symmetric Gauss-Seidel for pre-smoothing
	//Yang, X, Li, S, Yuan, F et al. (2023) Optimizing Multi-grid Computation
	//and Parallelization on Multi-cores.
	Eigen::VectorXd p_1;
	for(long tc = 0; tc < 1; tc ++){
		Eigen::VectorXd p_0 = - consUppe[tempLeve] * resuSolu;
		Eigen::VectorXd bpp0 = righHand + p_0;
		for(long ti = 0; ti < consStif[tempLeve].rows(); ti ++){
			double tempSumm = consLowe[tempLeve].row(ti) * resuSolu;
			resuSolu(ti) = (bpp0(ti) - tempSumm) / consDiag[tempLeve].coeff(ti,ti);
		}
		p_1 = consDiag[tempLeve] * resuSolu - p_0;
		for(long ti = consStif[tempLeve].rows() - 1; ti >= 0; ti --){
			double tempSumm = consUppe[tempLeve].row(ti) * resuSolu;
			resuSolu(ti) = (p_1(ti) - tempSumm) / consDiag[tempLeve].coeff(ti,ti);
		}
	}
	// for(long tc = 0; tc < 1; tc ++){
		// for(long ti = 0; ti < consStif[tempLeve].rows(); ti ++){
			// double tempSumm = consStif[tempLeve].row(ti) * resuSolu;
				// // - consStif[tempLeve].coeff(ti,ti) * resuSolu(ti);
			// resuSolu(ti) = (righHand(ti) - tempSumm) / consStif[tempLeve].coeff(ti,ti) 
				// + resuSolu(ti);
		// }
		// for(long ti = consStif[tempLeve].rows() - 1; ti >= 0; ti --){
			// double tempSumm = consStif[tempLeve].row(ti) * resuSolu;
			// resuSolu(ti) = (righHand(ti) - tempSumm) / consStif[tempLeve].coeff(ti,ti) 
				// + resuSolu(ti);
		// }
	// }
	//coarse grid correction
	Eigen::VectorXd resiErro = righHand - (p_1 + consLowe[tempLeve] * resuSolu);
	Eigen::VectorXd resuSolu_1 = Eigen::VectorXd::Zero(consStif[tempLeve - 1].rows());
	MULT_VCYC(
		tempLeve - 1, 
		realProl[tempLeve - 1].transpose() * resiErro, 
		resuSolu_1, 
		direSolv
	);
	resuSolu = resuSolu + realProl[tempLeve - 1] * resuSolu_1;
	//Symmetric Gauss-Seidel for post-smoothing
	for(long tc = 0; tc < 1; tc ++){
		Eigen::VectorXd p_0 = - consUppe[tempLeve] * resuSolu;
		Eigen::VectorXd bpp0 = righHand + p_0;
		for(long ti = 0; ti < consStif[tempLeve].rows(); ti ++){
			double tempSumm = consLowe[tempLeve].row(ti) * resuSolu;
			resuSolu(ti) = (bpp0(ti) - tempSumm) / consDiag[tempLeve].coeff(ti,ti);
		}
		Eigen::VectorXd p_1 = consDiag[tempLeve] * resuSolu - p_0;
		for(long ti = consStif[tempLeve].rows() - 1; ti >= 0; ti --){
			double tempSumm = consUppe[tempLeve].row(ti) * resuSolu;
			resuSolu(ti) = (p_1(ti) - tempSumm) / consDiag[tempLeve].coeff(ti,ti);
		}
	}
	// for(long tc = 0; tc < 1; tc ++){
		// for(long ti = 0; ti < consStif[tempLeve].rows(); ti ++){
			// double tempSumm = consStif[tempLeve].row(ti) * resuSolu;
			// resuSolu(ti) = (righHand(ti) - tempSumm) / consStif[tempLeve].coeff(ti,ti) 
				// + resuSolu(ti);
		// }
		// for(long ti = consStif[tempLeve].rows() - 1; ti >= 0; ti --){
			// double tempSumm = consStif[tempLeve].row(ti) * resuSolu;
			// resuSolu(ti) = (righHand(ti) - tempSumm) / consStif[tempLeve].coeff(ti,ti) 
				// + resuSolu(ti);
		// }
	// }
	return 1;
}

long MGPIS::MULT_SOLV(const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu){
	OUTPUT_TIME("MGPIS::MULT_SOLV");
	//
	resuSolu = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
	long maxiNumb = 10000;
	double toleLimi = 1.0E-14 * totaForc.norm();
	DIRE_SOLV direSolv;
	direSolv.compute(consStif[0]);
	//
	long iterNumb = 0;
	Eigen::VectorXd resiErro = totaForc - consStif[maxiLeve] * resuSolu;
	std::vector<double> moniErro(5);
	while(iterNumb < maxiNumb){
		MULT_VCYC(maxiLeve, totaForc, resuSolu, direSolv);
		resiErro = totaForc - consStif[maxiLeve] * resuSolu;
		//
		moniErro[iterNumb % moniErro.size()] = resiErro.norm();
		if(iterNumb >= moniErro.size() - 1){
			double erroMedi, erroOsci;
			VECT_MEDI_OSCI(moniErro, erroMedi, erroOsci);
			if(erroOsci < 0.1 * erroMedi){
				break;
			}
		}
		iterNumb ++;
	}
	std::cout << "#Iteration: " << iterNumb << ", residual: " 
		<< moniErro[iterNumb % moniErro.size()] << "/" << toleLimi;
	OUTPUT_TIME(":");
	return 1;
}

//Jonathan Richard Shewchuk. An Introduction to the Conjugate Gradient Method Without the Agonizing Pain. 1994
long MGPIS::CG_SOLV(long precSwit, const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu){
	std::cout << "MGPIS::CG_SOLV";
	if(precSwit == 0){
		std::cout << " (diagonal preconditioner)";
	}
	else if(precSwit == 1){
		std::cout << " (multigrid preconditioner)";
	}
	OUTPUT_TIME("");
	//
	resuSolu = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
	long maxiNumb = resuSolu.rows();
	double toleLimi = 1.0E-14 * totaForc.norm();
	//
	//diagonal preconditioner
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> diagPrec(consStif[maxiLeve].rows());
	//multigrid preconditioner
	DIRE_SOLV direSolv;
	if(precSwit == 0){
		DIAG_PREC(consStif[maxiLeve], diagPrec);
	}
	else if(precSwit == 1){
		direSolv.compute(consStif[0]);
	}
	//
	long iterNumb = 0;
	Eigen::VectorXd resiErro = totaForc - consStif[maxiLeve] * resuSolu;
	Eigen::VectorXd searDire = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
	if(precSwit == 0){
		searDire = diagPrec * resiErro;
	}
	else if(precSwit == 1){
		MULT_VCYC(maxiLeve, resiErro, searDire, direSolv);
	}
	double delt_new = resiErro.transpose() * searDire;
	while(iterNumb < maxiNumb && resiErro.norm() > toleLimi){
		Eigen::VectorXd matrDire;
		matrDire = consStif[maxiLeve] * searDire;
		double alph = delt_new / (searDire.transpose() * matrDire);
		resuSolu = resuSolu + alph * searDire;
		resiErro = resiErro - alph * matrDire;
		Eigen::VectorXd precDire = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
		if(precSwit == 0){
			precDire = diagPrec * resiErro;
		}
		else if(precSwit == 1){
			MULT_VCYC(maxiLeve, resiErro, precDire, direSolv);
		}
		double delt_old = delt_new;
		delt_new = resiErro.transpose() * precDire;
		double beta = delt_new / delt_old;
		searDire = precDire + beta * searDire;
		if(iterNumb % 100 == 99){
			std::cout << "#Iteration: " << iterNumb << ", residual: " 
				<< resiErro.norm() << "/" << toleLimi << std::endl;
		}
		iterNumb ++;
	}
	std::cout << "#Iteration: " << iterNumb - 1 
		<< ", residual: " << resiErro.norm() << "/" << toleLimi;
	OUTPUT_TIME(":");
	return 1;
}

long MGPIS::GMRES_SOLV(long precSwit, const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu){
	std::cout << "MGPIS::GMRES_SOLV";
	if(precSwit == 0){
		std::cout << " (diagonal preconditioner)";
	}
	else if(precSwit == 1){
		std::cout << " (multigrid preconditioner)";
	}
	OUTPUT_TIME("");
	//diagonal preconditioner
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> diagPrec(consStif[maxiLeve].rows());
	//multigrid preconditioner
	DIRE_SOLV direSolv;
	if(precSwit == 0){
		DIAG_PREC(consStif[maxiLeve], diagPrec);
	}
	else if(precSwit == 1){
		std::cout << "DOF of level 0 = " << consStif[0].rows() << std::endl;
		direSolv.compute(consStif[0]);
	}
	//
	resuSolu = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
	long maxiNumb = resuSolu.rows();
	double toleLimi = 1.0E-12 * totaForc.norm();
	//
	long iterNumb = 0;
	long iterStag = 10;
	std::vector<double> moniErro(iterStag);//size <= iterStag
	Eigen::VectorXd resuSolu_0;
	Eigen::VectorXd origErro;
	double normR_0;
	Eigen::MatrixXd supeHess;
	Eigen::MatrixXd orthBasi;
	Eigen::MatrixXd Q, R;
	while(iterNumb < maxiNumb){
		//restart
		if(iterNumb % iterStag == 0){
			resuSolu_0 = resuSolu;
			origErro = totaForc - consStif[maxiLeve] * resuSolu_0;
			Eigen::VectorXd precErro = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
			if(precSwit == 0){
				precErro = diagPrec * origErro;
			}
			else if(precSwit == 1){
				MULT_VCYC(maxiLeve, origErro, precErro, direSolv);
			}
			normR_0 = precErro.norm();
			supeHess.conservativeResize(0, 0);
			orthBasi = precErro / normR_0;
		}
		//
		Eigen::VectorXd origV = consStif[maxiLeve] * orthBasi.col(iterNumb % iterStag);
		Eigen::VectorXd precV = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
		if(precSwit == 0){
			precV = diagPrec * origV;
		}
		else if(precSwit == 1){
			MULT_VCYC(maxiLeve, origV, precV, direSolv);
		}
		Eigen::VectorXd b_i = orthBasi.transpose() * precV;
		Eigen::VectorXd q_ip1 = precV - orthBasi * b_i;
		double normQip1 = q_ip1.norm();
		supeHess.conservativeResize(b_i.rows() + 1, supeHess.cols() + 1);
		supeHess.block(0,supeHess.cols() - 1,b_i.rows(),1) = b_i;
		supeHess.block(b_i.rows(),0,1,supeHess.cols() - 1) = 
			Eigen::MatrixXd::Zero(1,supeHess.cols() - 1);
		supeHess(b_i.rows(), supeHess.cols() - 1) = normQip1;
		q_ip1 = q_ip1 / normQip1;
		orthBasi.conservativeResize(orthBasi.rows(), orthBasi.cols() + 1);
		orthBasi.col(orthBasi.cols() - 1) = q_ip1;
		if(iterNumb % iterStag == 0){
			Q = supeHess / supeHess.norm();
			R.conservativeResize(1, 1);
			R(0,0) = supeHess.norm();
		}
		else{
			Q.conservativeResize(Q.rows() + 1, Q.cols() + 1);
			R.conservativeResize(R.rows() + 1, R.cols() + 1);
			Q.block(Q.rows() - 1,0,1,Q.cols() - 1) = 
				Eigen::MatrixXd::Zero(1,Q.cols() - 1);
			R.block(R.rows() - 1,0,1,R.cols() - 1) = 
				Eigen::MatrixXd::Zero(1,R.cols() - 1);
			R.block(0,R.cols() - 1,R.rows() - 1,1) = 
				Q.block(0,0,Q.rows(),Q.cols() - 1).transpose() 
				* supeHess.col(supeHess.cols() - 1);
			Q.col(Q.cols() - 1) = supeHess.col(supeHess.cols() - 1) 
				- Q.block(0,0,Q.rows(),Q.cols() - 1) * R.block(0,R.cols() - 1,R.rows() - 1,1);
			R(R.rows() - 1,R.cols() - 1) = Q.col(Q.cols() - 1).norm();
			Q.col(Q.cols() - 1) = Q.col(Q.cols() - 1) / R(R.rows() - 1,R.cols() - 1);
		}
		Eigen::VectorXd realRigh = normR_0 * Q.row(0).transpose();
		Eigen::VectorXd y = Eigen::VectorXd::Zero(R.rows());
		y(R.rows() - 1) = realRigh(R.rows() - 1) / R(R.rows() - 1, R.cols() - 1);
		for(long tj = R.rows() - 2; tj >= 0; tj --){
			double tempSumm = (R.block(tj, tj + 1, 1, R.cols() - 1 - tj) 
				* y.block(tj + 1, 0, R.cols() - 1 - tj, 1))(0,0);
			y(tj) = (realRigh(tj) - tempSumm) / R(tj,tj);
		}
		resuSolu = resuSolu_0 + orthBasi.block(0, 0, orthBasi.rows(), y.rows()) * y;
		origErro = totaForc - consStif[maxiLeve] * resuSolu;
		//
		moniErro[iterNumb % moniErro.size()] = origErro.norm();
		if(iterNumb % 100 == 99){
			std::cout << "#Iteration: " << iterNumb << ", residual: " 
				<< moniErro[iterNumb % moniErro.size()] << "/" << toleLimi << std::endl;
		}
		if(iterNumb >= moniErro.size() - 1){
			double erroMedi, erroOsci;
			VECT_MEDI_OSCI(moniErro, erroMedi, erroOsci);
			if(moniErro[iterNumb % moniErro.size()] <= toleLimi || 
				(moniErro[iterNumb % moniErro.size()] <= 1.0E2 * toleLimi 
				&& erroOsci < 0.1 * erroMedi)){
				break;
			}
		}
		iterNumb ++;
	}
	std::cout << "#Iteration: " << iterNumb << ", residual: " 
		<< moniErro[iterNumb % moniErro.size()] << "/" << toleLimi;
	OUTPUT_TIME(":");
	return 1;
}

long MGPIS::BiCGSTAB_SOLV(long precSwit, 
	const Eigen::VectorXd &totaForc, Eigen::VectorXd &resuSolu){
	std::cout << "MGPIS::BiCGSTAB_SOLV";
	if(precSwit == 0){
		std::cout << " (diagonal preconditioner)";
	}
	else if(precSwit == 1){
		std::cout << " (multigrid preconditioner)";
	}
	OUTPUT_TIME("");
	//
	resuSolu = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
	long maxiNumb = resuSolu.rows();
	double toleLimi = 1.0E-14 * totaForc.norm();
	//
	//diagonal preconditioner
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> diagPrec(consStif[maxiLeve].rows());
	//multigrid preconditioner
	DIRE_SOLV direSolv;
	if(precSwit == 0){
		DIAG_PREC(consStif[maxiLeve], diagPrec);
	}
	else if(precSwit == 1){
		direSolv.compute(consStif[0]);
	}
	//
	long iterNumb = 0;
	Eigen::VectorXd resiErro = totaForc - consStif[maxiLeve] * resuSolu;
	Eigen::VectorXd resiOver = resiErro;
	std::vector<double> rho(2);
	Eigen::VectorXd p, v;
	double alph, omeg;
	while(iterNumb < maxiNumb && resiErro.norm() > toleLimi){
		rho[(iterNumb - 1 + 2) % 2] = resiOver.dot(resiErro);
		if(abs(rho[(iterNumb - 1 + 2) % 2]) == 0.0){
			OUTPUT_TIME("MGPIS::BiCGSTAB_SOLV: ERROR 1");
			break;
		}
		if(iterNumb == 0){
			p = resiErro;
		}
		else{
			double beta = (rho[(iterNumb - 1 + 2) % 2] / rho[(iterNumb - 2 + 2) % 2]) 
				* (alph / omeg);
			p = resiErro + beta * (p - omeg * v);
		}
		Eigen::VectorXd prepP = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
		if(precSwit == 0){
			prepP = diagPrec * p;
		}
		else if(precSwit == 1){
			MULT_VCYC(maxiLeve, p, prepP, direSolv);
		}
		v = consStif[maxiLeve] * prepP;
		alph = rho[(iterNumb - 1 + 2) % 2] / resiOver.dot(v);
		Eigen::VectorXd s = resiErro - alph * v;
		if(s.norm() <= 0.0){
			resuSolu += alph * prepP;
			break;
		}
		Eigen::VectorXd prepS = Eigen::VectorXd::Zero(consStif[maxiLeve].rows());
		if(precSwit == 0){
			prepS = diagPrec * s;
		}
		else if(precSwit == 1){
			MULT_VCYC(maxiLeve, s, prepS, direSolv);
		}
		Eigen::VectorXd t = consStif[maxiLeve] * prepS;
		omeg = t.dot(s) / t.dot(t);
		resuSolu += alph * prepP + omeg * prepS;
		resiErro = s - omeg * t;
		//for continuation it is necessary that omeg != 0.0
		// if(iterNumb % 100 == 99){
			std::cout << "#Iteration: " << iterNumb << ", residual: " 
				<< resiErro.norm() << "/" << toleLimi << std::endl;
		// }
		iterNumb ++;
	}
	std::cout << "#Iteration: " << iterNumb - 1 
		<< ", residual: " << resiErro.norm() << "/" << toleLimi;
	OUTPUT_TIME(":");
	return 1;
}

#endif
