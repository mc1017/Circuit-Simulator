#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <complex>
#include "library/Eigen/Dense"

using namespace Eigen;

int main(){
    double Ieq, Geq, Vt, Vd, Is, Id;
    Is = 1 * pow(10, -15);
    Vt = 0.025875;
    Vd = 0.7;
    Id = Is * (exp(Vd/Vt) - 1);
    Geq = Is / Vt * exp(Vd/Vt);
    Ieq = Id - Geq * Vd;

    Matrix<double, 2, 2> A;
    Matrix<double, 2, 1> B;
    Matrix<double, 2, 1> X;

    for(int i = 0; i < 10; i++){
        A.setZero();
        B.setZero();
        X.setZero();
        A(0,0) = 0.0101;
        A(0,1) = -0.0001;
        B(0,0) = 0.1;
        B(1,0) = -Ieq;
        A(1,0) = -0.0001;
        A(1,1) = 0.0001+Geq;
        X = A.fullPivLu().solve(B);
        std::cout << A << std::endl;
        std::cout << std::endl;
        std::cout << B << std::endl;
        std::cout << std::endl;
        std::cout << X << std::endl;
        std::cout << std::endl;
        Vd = X(1,0);
        std::cout << Vd << std::endl;
        std::cout << std::endl;
        Id = Is * (exp(Vd/Vt) - 1);
        Geq = Is / Vt * exp(Vd/Vt);
        Ieq = Id - Geq * Vd;

    }

}