#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

int main()
{
    complex<double> c(1,1);

    MatrixXcf matrixA1;//for complex numbers we chose the type of matrix as being Xcf
    matrixA1.setZero(3,3);// can't define as a 3f matrix. Type of matrix is important.
    matrixA1 (0,0) = 30;
    matrixA1(0,1) = c;
    cout<<"\n"<<matrixA1<<endl;

    MatrixXcf mA2;
    mA2.setIdentity(3,3);
    cout<<"\n"<<mA2<<endl;
    MatrixXcf test = mA2*matrixA1;

    //Matrix3f out_res;

    cout<<"\n"<<test<<endl;

}

