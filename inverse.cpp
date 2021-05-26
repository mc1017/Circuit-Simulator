#include <iostream>
#include <complex>
#include "eigen/Eigen/Dense"

using namespace std;
using namespace Eigen;


 
int main(){
    int A_rows, A_columns;
    int b_rows;
    const std::complex<double> If(0.0, 1.0);
    A_rows =2;
    A_columns = 2;
    b_rows = 2;
    MatrixXcd A(A_rows,A_columns);
    MatrixXcd b(b_rows, 1);

    A << 1.0 + 2.0 * If, 2.0 + 1.0 * If, 3.0 - 1.0 * If, 4.0 - 2.0 * If;
    b << 16.0 ,-8.0 + 2.0 * If ;
    MatrixXcd x = A.fullPivLu().solve(b);
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the vector b:\n" << b << endl;
    cout << "The solution is:\n" << x << endl;
}