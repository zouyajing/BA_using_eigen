#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

// Kinematics functions
Matrix3d skew(Vector3d v);
Matrix3d fast_skewexp(Vector3d v);
Vector3d skewcoords(Matrix3d M);
Matrix3d skewlog(Matrix3d M);
MatrixXd kroen_product(MatrixXd A, MatrixXd B);
Matrix3d v_logmap(VectorXd x);
MatrixXd diagonalMatrix(MatrixXd M, unsigned int N);


Matrix4d inverse_se3(Matrix4d T);
Matrix4d expmap_se3(Vector6d x);
Vector6d logmap_se3(Matrix4d T);
Matrix6d adjoint_se3(Matrix4d T);
Matrix6d uncTinv_se3(Matrix4d T, Matrix6d covT );
Matrix6d unccomp_se3(Matrix4d T1, Matrix6d covT1, Matrix6d covTinc );
Vector6d reverse_se3(Vector6d x);


Vector3d logarithm_map_so3(Matrix3d R);
MatrixXd der_logarithm_map(Matrix4d T);
MatrixXd der_logarithm_map_appr(Matrix4d T, double delta);

double robustWeightCauchy(double norm_res);







#endif // HELPER_H
