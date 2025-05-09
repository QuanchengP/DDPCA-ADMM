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
 * \file   DEHW.cpp
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  Implementation.
 */

#include "DEHW.h"
#include "DEHW_1.h"

long EXAM_DEHW();//example DEHW
long ISNO_EIGE();//is/not eigen analysis
long TANG_PEPA();//tangential penalty parameter
long ISNO_SELO();//is/not self-locking
long SHOW_ABBR();//show abbreviation
long SELE_COSP();//select multigrid strategy
long SELE_COSP_1();

int main(int argc, char **argv){
	omp_set_nested(1);//nested parallelism
	omp_set_dynamic(1);//dynamic thread adjustment
	//
	long isnoSelo = ISNO_SELO();
	std::cout << std::endl;
	std::cout << "******Double enveloping hourglass worm drive example****" << std::endl;
	std::cout << "*0 - ADMM without DD: 1/1 domain for worm/wheel        *" << std::endl;
	std::cout << "*1 - ADMM with DD: 34/18 domains for worm/wheel        *" << std::endl;
	std::cout << "*2 - Dual mortar method with GMG-BiCGSTAB solver       *" << std::endl;
	std::cout << "*3 - ADMM with DD: 34/18 domains for worm/wheel (type  *" << std::endl;
	std::cout << "*    III cross corner)                                 *" << std::endl;
	SHOW_ABBR();
	std::cout << "Please enter a number in 0~3: ";
	long caid;
	std::cin >> caid;
	switch(caid){
		case 0:{
			TANG_PEPA();
			SELE_COSP_1();
			DEHW solv(isnoSelo);
			solv.SOLVE(1,0);
			break;
		}
		case 1:{
			TANG_PEPA();
			DEHW solv(isnoSelo);
			long isnoEige = ISNO_EIGE();
			switch(isnoEige){
				case 0:{
					SELE_COSP_1();
					solv.SOLVE();
					break;
				}
				case 1:{
					solv.SOLVE(-1,1);
					break;
				}
				case 2:{
					solv.SOLVE(0,1);
					break;
				}
			}
			break;
		}
		case 2:{
			DEHW solv(isnoSelo);
			solv.SOLVE(2,0);
			break;
		}
		case 3:{
			TANG_PEPA();
			SELE_COSP_1();
			DEHW_1 solv(isnoSelo);
			solv.SOLVE();
			break;
		}
	}
	return 1;
}

long ISNO_EIGE(){
	std::cout << std::endl;
	std::cout << "***********************Selection************************" << std::endl;
	std::cout << "*0 - Contact analysis                                  *" << std::endl;
	std::cout << "*1 - Eigen analysis of global problem                  *" << std::endl;
	std::cout << "*2 - Eigen analysis of global coarse problem           *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~2: ";
	long caid;
	std::cin >> caid;
	return caid;
}

long TANG_PEPA(){
	std::cout << std::endl;
	std::cout << "*******Coefficient of tangential penalty parameter******" << std::endl;
	std::cout << "*0 - f_t^worm,wheel = 1                                *" << std::endl;
	std::cout << "*1 - f_t^worm,wheel = 10                               *" << std::endl;
	std::cout << "*2 - f_t^worm,wheel = 100                              *" << std::endl;
	std::cout << "*3 - f_t^worm,wheel = 1000                             *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~2: ";
	long caid;
	std::cin >> caid;
	switch(caid){
		case 0:{
			tapeCoef = 0.5;
			break;
		}
		case 1:{
			tapeCoef = 5.0;
			break;
		}
		case 2:{
			tapeCoef = 50.0;
			break;
		}
		case 3:{
			tapeCoef = 500.0;
			break;
		}
	}
	return 1;
}

long ISNO_SELO(){
	std::cout << std::endl;
	std::cout << "******Double enveloping hourglass worm drive example****" << std::endl;
	std::cout << "*0 - Contact analysis with driving worm                *" << std::endl;
	std::cout << "*1 - Slef-locking analysis with driving wheel          *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~1: ";
	long caid;
	std::cin >> caid;
	// if(caid == 0){
		// tapeCoef = 5.0;
	// }
	// else{
		// tapeCoef = 50.0;
	// }
	return 1 - caid;
}

long SHOW_ABBR(){
	std::cout << "********************************************************" << std::endl;
	std::cout << "*ADMM: alternating direction method of multiplier      *" << std::endl;
	std::cout << "*BiCGSTAB: biconjugate gradient stabilized             *" << std::endl;
	std::cout << "*CG: conjugate gradient, DD: domain decomposition      *" << std::endl;
	std::cout << "*GMG: geometric multigrid (preconditioner)             *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	return 1;
}

long SELE_COSP_1(){
	std::cout << std::endl;
	std::cout << "************Type of coarse space correction*************" << std::endl;
	std::cout << "*0 - Macroscopic problem in LATIN                      *" << std::endl;
	std::cout << "*1 - Without coarse space correction                   *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~1: ";
	long caid;
	std::cin >> caid;
	if(caid == 0){
		whadCosp = (1 << 0);
		return (1 << 0);
	}
	else{
		whadCosp = 0;
		return 0;
	}
	return caid;
}
