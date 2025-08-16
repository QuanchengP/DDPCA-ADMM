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
 * \file   PREP.h
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  Pre declaration of common variable, function, class.
 */

#ifndef _PREP_H
#define _PREP_H

/*************************************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <functional>//lambda expressions
#include <iterator>
#include <utility>
#include <ctime>
#include <sys/stat.h>//for windows without cygwin: <direct.h>
#include <omp.h>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/SparseLU"
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseSymMatProd.h"

/*************************************************************************************************/
const double PI = acos(-1.0);
//the maximum number of subdomains
const long MAXI_DOMA_NUMB = 1000;
//the maximum number of interfaces
const long MAXI_INTE_NUMB = 1000;
//when the DOF is larger than DIRE_MAXI, 
//use the Conjugate Gradient solver but not the direct solver
const long DIRE_MAXI = 120000;//macroscopic problem
const long DIRE_MAXI_SUBD = 50000;//subproblem
//when the DOF is larger than COGR_MAXI, 
//use the multigrid preconditioned Conjugate Gradient solver
const long COGR_MAXI = 100000;//almost never been used
//the maximum implementation number of CSC in contact analysis
long MULT_MAXI = 1000;//almost never been used
//control output of warning and error information: the larger, the more output
long debuMode = 0;

//output string with current time
long OUTPUT_TIME(std::string outpStri);

//different numerical examples, different output directories
std::string outpDire = "";//extern
//concatenate outpDire with fileName
std::string DIRECTORY(std::string fileName);

//output outpStri + current time
long OUTPUT_TIME(std::string outpStri){
	time_t now = time(0);
	char *cnow = ctime(&now);
	std::cout << outpStri << ": " << cnow << std::endl;
	return 1;
}

//add prefix "outpDire" to "fileName"
std::string DIRECTORY(std::string fileName){
	std::stringstream tempStre;
	tempStre << outpDire << fileName;
	return tempStre.str();
}

typedef std::vector<std::vector<long>> VECTOR2L;
typedef std::vector<std::vector<double>> VECTOR2D;
typedef std::vector<std::vector<std::vector<long>>> VECTOR3L;
typedef std::vector<std::vector<std::vector<double>>> VECTOR3D;
typedef std::vector<std::vector<std::vector<std::vector<long>>>> VECTOR4L;
typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<double,Eigen::RowMajor>> DIRE_SOLV;
typedef Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> DILU_SOLV;
typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, 
	Eigen::Lower|Eigen::Upper> COGR_SOLV;
typedef Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator RSPA_INNE;
typedef Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator CSPA_INNE;
//matrix about sparse matrix
typedef std::vector<std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>>> SPARMATS;

//resize of two dimensional std::vector
template<typename T>
long VECT_RESI(T &tempVect, long size_0, long size_1){
	tempVect.resize(size_0);
	for(long ti = 0; ti < size_0; ti ++){
		tempVect[ti].resize(size_1);
	}
	return 1;
}

//assignment of two dimensional std::vector
template<typename T_0, typename T_1>
long VECT_ASSI(T_0 &tempVect, long size_0, long size_1, T_1 valu_0){
	tempVect.resize(size_0);
	for(long ti = 0; ti < size_0; ti ++){
		tempVect[ti].assign(size_1, valu_0);
	}
	return 1;
}

//display row major sparse matrix
long DISP_RSPA(const Eigen::SparseMatrix<double,Eigen::RowMajor> &inpuSpar){
	for(long ti = 0; ti < inpuSpar.rows(); ti ++){
		for(RSPA_INNE iterStif(inpuSpar, ti); iterStif; ++ iterStif){
			std::cout << ti << "," << iterStif.col()
				<< "," << iterStif.value() << std::endl;
		}
	}
	return 1;
}

//median and oscilation of vector
long VECT_MEDI_OSCI(const std::vector<double> &inpuVect, double &resuMedi, double &resuOsci){
	auto resuMaxi = std::max_element(inpuVect.begin(), inpuVect.end());
	auto resuMini = std::min_element(inpuVect.begin(), inpuVect.end());
	resuMedi = (*resuMaxi + *resuMini) / 2.0;
	resuOsci = *resuMaxi - *resuMini;
	return 1;
}

/*************************************************************************************************/
//all lines of a hexahedron
VECTOR2L hexaLine = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0}, 
	{0, 4}, {1, 5}, {2, 6}, {3, 7}, 
	{4, 5}, {5, 6}, {6, 7}, {7, 4}
};
//all faces of a hexahedron, the normal directions to the outside
VECTOR2L hexaFace = {
	{0, 3, 2, 1}, {4, 5, 6, 7}, 
	{0, 4, 7, 3}, {1, 2, 6, 5}, 
	{0, 1, 5, 4}, {3, 7, 6, 2}
};

//3d coordinate, two coordinates can be compared
class COOR : public std::vector<double>{
public:
	COOR(){
	}
	COOR(double temp_x, double temp_y, double temp_z){
		(*this).resize(3);
		(*this)[0] = temp_x;
		(*this)[1] = temp_y;
		(*this)[2] = temp_z;
	}
	bool operator<(const COOR &coor_1) const{
		for(long ti = 0; ti < 3; ti ++){
			if((*this)[ti] < coor_1[ti] - 1.0E-10) return true;
			else if((*this)[ti] > coor_1[ti] + 1.0E-10) return false;
		}
		return false;
	}
	COOR operator+(const COOR &coor_1) const{
		COOR coor_2(0.0, 0.0, 0.0);
		for(long ti = 0; ti < 3; ti ++){
			coor_2[ti] = (*this)[ti] + coor_1[ti];
		}
		return coor_2;
	}
	COOR operator/(const double &divi) const{
		COOR coor_2(0.0, 0.0, 0.0);
		for(long ti = 0; ti < 3; ti ++){
			coor_2[ti] = (*this)[ti] / divi;
		}
		return coor_2;
	}
};

/*************************************************************************************************/
//octree for global/local mesh refinement
class TREE_ELEM{
public:
	long parent;
	std::vector<long> cornNode;
	long level;
	//refinement patter: 0 - xi,eta,zeta; 1 - xi,eta; 2 - eta,zeta; 3 - zeta,xi;
	//4 - xi; 5 - eta; 6 - zeta; 7 - not to be refined;
	long refiPatt;
	std::vector<long> children;
	TREE_ELEM(){
		parent = -1;
		cornNode.resize(8);
		level = 0;
		refiPatt = 7;
		children.resize(0);
	}
};

/*************************************************************************************************/
//numerical Gauss quadrature over trilinear element
class TRILINEAR_QUADRATURE{
public:
	long numbNgip;//number of numerical Gauss integration point
	VECTOR2D nacoNgip;//natural coordinate of Ngip
	std::vector<double> niwfNgip;//numerical integration weight factor of Ngip
	VECTOR2D nacoCorn;//natural coordinate of eight corner nodes
	//shape function of numerical Gauss integration point
	Eigen::MatrixXd shfuNgip;
	//partial derivative of N relative to xi/eta/zeta
	std::vector<Eigen::MatrixXd> pnpxNgip;
	TRILINEAR_QUADRATURE(){
		//{-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};
		std::vector<double> ngipList = {-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};
		//{1.0, 1.0};
		std::vector<double> niwfList = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
		nacoCorn = {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, 
			{1.0, 1.0, -1.0}, {-1.0, 1.0, -1.0}, 
			{-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0}, 
			{1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
		};
		numbNgip = ngipList.size() * ngipList.size() * ngipList.size();
		VECT_RESI(nacoNgip, numbNgip, 3);
		niwfNgip.resize(numbNgip);
		shfuNgip.resize(numbNgip,8);
		pnpxNgip.resize(numbNgip);
		for(long ti = 0; ti < numbNgip; ti ++){
			pnpxNgip[ti].resize(3,8);
		}
		for(long ti = 0; ti < ngipList.size(); ti ++){
			for(long tj = 0; tj < ngipList.size(); tj ++){
				for(long tk = 0; tk < ngipList.size(); tk ++){
					long tempIijk = ti * ngipList.size() * ngipList.size() 
						+ tj * ngipList.size() + tk;
					nacoNgip[tempIijk][0] = ngipList[ti];
					nacoNgip[tempIijk][1] = ngipList[tj];
					nacoNgip[tempIijk][2] = ngipList[tk];
					niwfNgip[tempIijk] = niwfList[ti] * niwfList[tj] * niwfList[tk];
				}
			}
		}
		for(long ti = 0; ti < numbNgip; ti ++){
			for(long tj = 0; tj < 8; tj ++){
				pnpxNgip[ti](0,tj) = nacoCorn[tj][0] 
					* (1.0 + nacoCorn[tj][1] * nacoNgip[ti][1]) 
					* (1.0 + nacoCorn[tj][2] * nacoNgip[ti][2]) / 8.0;
				pnpxNgip[ti](1,tj) = (1.0 + nacoCorn[tj][0] * nacoNgip[ti][0]) 
					* nacoCorn[tj][1] 
					* (1.0 + nacoCorn[tj][2] * nacoNgip[ti][2]) / 8.0;
				pnpxNgip[ti](2,tj) = (1.0 + nacoCorn[tj][0] * nacoNgip[ti][0]) 
					* (1.0 + nacoCorn[tj][1] * nacoNgip[ti][1]) 
					* nacoCorn[tj][2] / 8.0;
				shfuNgip(ti,tj) = (1.0 + nacoCorn[tj][0] * nacoNgip[ti][0]) 
					* (1.0 + nacoCorn[tj][1] * nacoNgip[ti][1]) 
					* (1.0 + nacoCorn[tj][2] * nacoNgip[ti][2]) / 8.0;
			}
		}
	}
}trilQuad;

//numerical Gauss quadrature over bilinear element
class BILINEAR_QUADRATURE{
public:
	long numbNgip;//number of numerical Gauss integration point
	VECTOR2D nacoNgip;//natural coordinate of Ngip
	std::vector<double> niwfNgip;//numerical integration weight factor of Ngip
	VECTOR2D nacoCorn;//natural coordinate of eight corner nodes
	//shape function of numerical Gauss integration point
	std::vector<Eigen::Matrix<double,3,12>> shfuNgip;
	//partial derivative of N relative to xi/eta/zeta
	std::vector<Eigen::MatrixXd> pnpxNgip;
	BILINEAR_QUADRATURE(){
		std::vector<double> ngipList = {-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};
		std::vector<double> niwfList = {1.0, 1.0};
		nacoCorn = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
		numbNgip = ngipList.size() * ngipList.size();
		VECT_RESI(nacoNgip, numbNgip, 2);
		niwfNgip.resize(numbNgip);
		shfuNgip.resize(numbNgip);
		pnpxNgip.resize(numbNgip);
		for(long ti = 0; ti < numbNgip; ti ++){
			pnpxNgip[ti].resize(2,4);
		}
		for(long ti = 0; ti < ngipList.size(); ti ++){
			for(long tj = 0; tj < ngipList.size(); tj ++){
				long tempIijk = ti * ngipList.size() + tj;
				nacoNgip[tempIijk][0] = ngipList[ti];
				nacoNgip[tempIijk][1] = ngipList[tj];
				niwfNgip[tempIijk] = niwfList[ti] * niwfList[tj];
			}
		}
		for(long ti = 0; ti < numbNgip; ti ++){
			Eigen::Matrix<double,3,12> shap_ti = Eigen::MatrixXd::Zero(3,12);
			for(long tj = 0; tj < 4; tj ++){
				pnpxNgip[ti](0,tj) = nacoCorn[tj][0] 
					* (1.0 + nacoCorn[tj][1] * nacoNgip[ti][1]) / 4.0;
				pnpxNgip[ti](1,tj) = (1.0 + nacoCorn[tj][0] * nacoNgip[ti][0]) 
					* nacoCorn[tj][1] / 4.0;
				double shap_tj = (1.0 + nacoCorn[tj][0] * nacoNgip[ti][0])
					* (1.0 + nacoCorn[tj][1] * nacoNgip[ti][1]) / 4.0;
				for(long tk = 0; tk < 3; tk ++){
					shap_ti(tk, 3 * tj + tk) = shap_tj;
				}
			}
			shfuNgip[ti] = shap_ti;
		}
	}
}biliQuad;

//intergral boundary transformation: from [-1, 1] * [-1, 1] to [0, 1] * [0, 1-x]
//H.T. Rathod, K.V. Nagaraja, B. Venkatesudu, N.L. Ramesh.
//Gauss Legendre quadrature over a triangle. IIS, 2004, 84: 183-188.
class TRIANGLE_QUADRATURE{
public:
	long numbNgip;//number of numerical Gauss integration point
	VECTOR2D nacoNgip;//natural coordinate of Ngip
	std::vector<double> niwfNgip;//numerical integration weight factor of Ngip
	TRIANGLE_QUADRATURE(){
		//quadrature point = - sqrt(3.0) / 3.0, quadrature weight = 1.0
		std::vector<double> ngipList = {-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};
		std::vector<double> niwfList = {1.0, 1.0};
		numbNgip = ngipList.size() * ngipList.size();
		VECT_RESI(nacoNgip, numbNgip, 3);
		niwfNgip.resize(numbNgip);
		for(long ti = 0; ti < ngipList.size(); ti ++){
			for(long tj = 0; tj < ngipList.size(); tj ++){
				nacoNgip[ti * ngipList.size() + tj][0] = (1.0 + ngipList[ti]) / 2.0;
				nacoNgip[ti * ngipList.size() + tj][1] = (1.0 - ngipList[ti]) 
					* (1.0 + ngipList[tj]) / 4.0;
				nacoNgip[ti * ngipList.size() + tj][2] = 
					1.0 - nacoNgip[ti * ngipList.size() + tj][0] 
					- nacoNgip[ti * ngipList.size() + tj][1];
				niwfNgip[ti * ngipList.size() + tj] = (1.0 - ngipList[ti]) 
					/ 8.0 * niwfList[ti] * niwfList[tj];
			}
		}
	}
}triaQuad;

//bilinear quadrature, coordinate transformation Jacobian
long BIQU_TRJA(const Eigen::Vector2d &tempXiet, 
	const std::vector<Eigen::Vector3d> &elemCoor, double &weigFact){
	Eigen::Matrix<double,3,2> PrmaPxie = Eigen::MatrixXd::Zero(3,2);
	for(long ti = 0; ti < 3; ti ++){
		for(long tj = 0; tj < 4; tj ++){
			PrmaPxie(ti,0) += elemCoor[tj](ti) * (biliQuad.nacoCorn[tj][0] / 4.0 
				+ biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] * tempXiet(1) / 4.0);
			PrmaPxie(ti,1) += elemCoor[tj](ti) * (biliQuad.nacoCorn[tj][1] / 4.0 
				+ biliQuad.nacoCorn[tj][0] * biliQuad.nacoCorn[tj][1] * tempXiet(0) / 4.0);
		}
	}
	weigFact = sqrt(
		pow(PrmaPxie(1,0) * PrmaPxie(2,1) - PrmaPxie(1,1) * PrmaPxie(2,0), 2.0) + 
		pow(PrmaPxie(2,0) * PrmaPxie(0,1) - PrmaPxie(2,1) * PrmaPxie(0,0), 2.0) + 
		pow(PrmaPxie(0,0) * PrmaPxie(1,1) - PrmaPxie(0,1) * PrmaPxie(1,0), 2.0)
	);
	return 1;
}

long BIQU_TRJA(const std::vector<double> &tempXiet, 
	const std::vector<Eigen::Vector3d> &elemCoor, double &weigFact){
	Eigen::Vector2d tempXiet_1;
	tempXiet_1 << tempXiet[0], tempXiet[1];
	BIQU_TRJA(tempXiet_1, elemCoor, weigFact);
	return 1;
}

/*************************************************************************************************/
//diagonal preconditioner
long DIAG_PREC(const Eigen::SparseMatrix<double,Eigen::RowMajor> &inpuMatr, 
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> &outpDipre){
	Eigen::VectorXd tempDiag = inpuMatr.diagonal();
	for(long ti = 0; ti < inpuMatr.rows(); ti ++){
		tempDiag(ti) = 1.0 / tempDiag(ti);
	}
	outpDipre.diagonal() = tempDiag;
	return 1;
}

#endif
