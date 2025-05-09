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
 * \file   TORSION.cpp
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  Implementation.
 */

#include "TORSION.h"

long EXAM_TORS();//example torsion
long SHOW_ABBR();//show abbreviation
long SELE_COSP();//select multigrid strategy

int main(int argc, char **argv){
	omp_set_nested(1);//nested parallelism
	omp_set_dynamic(1);//dynamic thread adjustment
	//
	std::cout << std::endl;
	std::cout << "***************Torsional cylinder example***************" << std::endl;
	std::cout << "*0 - ADMM with DD: 1*8*4 domains                       *" << std::endl;
	std::cout << "*1 - ADMM with DD: 1*16*4 domains                      *" << std::endl;
	std::cout << "*2 - ADMM with DD: 1*16*8 domains                      *" << std::endl;
	std::cout << "*3 - GMG-CG soler without DD                           *" << std::endl;
	SHOW_ABBR();
	std::cout << "Please enter a number in 0~3: ";
	long caid;
	std::cin >> caid;
	switch(caid){
		case 0:{
			TORSION tors(1);
			tors.muscSett = SELE_COSP();
			tors.domaNumb = {1, 8, 4};//each diviNumb is divisible by each domaNumb
			tors.doleMcsc.assign(tors.domaNumb[0] * tors.domaNumb[1] * tors.domaNumb[2], 2);
			tors.SOLVE();
			break;
		}
		case 1:{
			TORSION tors(1);
			tors.muscSett = SELE_COSP();
			tors.domaNumb = {1, 16, 4};//each diviNumb is divisible by each domaNumb
			tors.doleMcsc.assign(tors.domaNumb[0] * tors.domaNumb[1] * tors.domaNumb[2], 2);
			tors.SOLVE();
			break;
		}
		case 2:{
			TORSION tors(1);
			tors.muscSett = SELE_COSP();
			tors.domaNumb = {1, 16, 8};//each diviNumb is divisible by each domaNumb
			tors.doleMcsc.assign(tors.domaNumb[0] * tors.domaNumb[1] * tors.domaNumb[2], 2);
			tors.SOLVE();
			break;
		}
		case 3:{
			TORSION tors(0);
			tors.SOLVE();
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

long SELE_COSP(){
	std::cout << std::endl;
	std::cout << "************Type of coarse space correction*************" << std::endl;
	std::cout << "*0 - Macroscopic problem in LATIN                      *" << std::endl;
	std::cout << "*1 - Interface-eliminated global problem               *" << std::endl;
	std::cout << "*2 - Both 0 and 1                                      *" << std::endl;
	std::cout << "*3 - Without coarse space correction                   *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~3: ";
	long caid;
	std::cin >> caid;
	if(caid <= 1){
		return (1 << caid);
	}
	else if(caid == 2){
		return ((1 << 0) + (1 << 1));
	}
	else{
		return 0;
	}
	return caid;
}
