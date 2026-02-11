#include "../General/DenseMatrix.hpp"

#include <iostream>

void Test();

int main(/*int argc, char **argv*/) {
    std::cout << "!********************************************************************************!\n";
    Test();
}

void Test(){
    Ddpca::DenseMatrix A(3, 3, std::vector<Ddpca::Real>({1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 3.0, 0.0}));
    Ddpca::DenseMatrix B(3, 3, std::vector<Ddpca::Real>({0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0}));
    Ddpca::DenseMatrix C(3, 3);
    GEMTMT(A, B, C);
    std::cout << "A:\n";
    A.Output();
    std::cout << "B:\n";
    B.Output();
    std::cout << "C:\n";
    C.Output();
}