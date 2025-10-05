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
 * \file   CYLINDER.cpp
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  Implementation.
 */

#include "CYLINDER.h"
#include "CYLINDER_1.h"

long EXAM_CYLI();//example cylinder
long SHOW_ABBR();//show abbreviation
long SELE_COSP_1();

int main(int argc, char **argv){
	omp_set_nested(1);//nested parallelism
	omp_set_dynamic(1);//dynamic thread adjustment
	//
	std::cout << std::endl;
	std::cout << "**************Cylinder Hertz contact example************" << std::endl;
	std::cout << "*0 - ADMM without DD: 4*2 domain per cylinder          *" << std::endl;
	std::cout << "*1 - ADMM with DD: 8*2 domains per cylinder            *" << std::endl;
	std::cout << "*2 - ADMM with DD: 16*2 domains per cylinder           *" << std::endl;
	std::cout << "*3 - Dual mortar method with GMG-BiCGSTAB solver       *" << std::endl;
	SHOW_ABBR();
	std::cout << "Please enter a number in 0~3: ";
	long caid;
	std::cin >> caid;
	switch(caid){
		case 0:{
			CYLINDER cyli;
			cyli.copyNumb = 4;
			cyli.muscSett = SELE_COSP_1();
			cyli.SOLVE();
			break;
		}
		case 1:{
			// for(long ti = 1; ti <= 100; ti ++){
				CYLINDER cyli;
				cyli.copyNumb = 8;
				cyli.muscSett = /*(1 << 0);*/SELE_COSP_1();
				// cyli.charFact = ti;
				cyli.SOLVE();
				//
				// std::ofstream tempOfst(DIRECTORY("resuIterNumb.txt"), std::ios::app);
				// tempOfst << std::setw(10) << ti 
					// << std::setw(10) << cyli.iterNumbReco << std::endl;
				// tempOfst.close();
			// }
			break;
		}
		case 2:{
			CYLINDER cyli;
			cyli.copyNumb = 16;
			cyli.muscSett = SELE_COSP_1();
			cyli.SOLVE();
			break;
		}
		case 3:{
			CYLINDER_1 cyli;
			cyli.copyNumb = 1;
			cyli.SOLVE(2);
			break;
		}
	}
	return 1;
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
		return (1 << 0);
	}
	else{
		return 0;
	}
	return caid;
}
