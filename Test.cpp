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
 * \file   Test.cpp
 * \author Quancheng Peng <QuanchengPeng@sylu.edu.cn>
 * \brief  Implementation.
 */

#include "PREP.h"

int main(int argc, char **argv){
	std::cout << "************************Examples************************" << std::endl;
	std::cout << "*0 - Torsional cylinder example                        *" << std::endl;
	std::cout << "*1 - Flexible beam example                             *" << std::endl;
	std::cout << "*2 - Block patch test example                          *" << std::endl;
	std::cout << "*3 - Cylinder Hertz contact example                    *" << std::endl;
	std::cout << "*4 - Double enveloping hourglass worm drive example    *" << std::endl;
	std::cout << "********************************************************" << std::endl;
	std::cout << "Please enter a number in 0~4: ";
	long meth;
	std::cin >> meth;
	switch(meth){
		case 0:{
			int tempResu = system("./examples/TORSION");
			break;
		}
		case 1:{
			int tempResu = system("./examples/BEAM");
			break;
		}
		case 2:{
			int tempResu = system("./examples/BLOCK");
			break;
		}
		case 3:{
			int tempResu = system("./examples/CYLINDER");
			break;
		}
		case 4:{
			int tempResu = system("./examples/DEHW");
			break;
		}
	}
	return 1;
}
