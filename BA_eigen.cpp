#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "helper.h"
using namespace  std;
using namespace Eigen;
vector<Vector3d> p3d;
vector<Vector2d> p2d;
vector<Matrix4d> poses;
double fx = 520.9, fy = 521.0, cx = 325.1, cy = 249.7;
bool use_point = true;
bool readData(string p3d_file, string p2d_file) {
    ifstream f3d(p3d_file);
    if(!f3d)
    {
        cout<<"No file for 3D points."<<endl;
    }
    else {
        while(!f3d.eof())
        {
            double pt3[3] = {0};
            f3d >> pt3[0] >> pt3[1] >> pt3[2];
            Vector3d Point3d;
            Point3d << pt3[0], pt3[1], pt3[2];
            p3d.push_back(Point3d);
        }
    }

    ifstream f2d(p2d_file);
    if (!f2d) {
        cout << "No file for 2D pixels." << endl;
        return false;
    } else {
        while (!f2d.eof()) {
            double pt2[2] = {0};
            f2d >> pt2[0] >> pt2[1];
            Vector2d Point2d;
            Point2d << pt2[0], pt2[1];
            p2d.push_back(Point2d);
        }
    }
    assert(p3d.size() == p2d.size());
    cout << "The test data contains " << p3d.size() << " 3D-2D pairs." << endl;
    return true;
}

void optimizeFunctions(MatrixXd &J, VectorXd &r, double &e )
{

    // define Jacobians, and residuals
    int n_poses = poses.size();
    int n_points = p3d.size();
    int n_pixels = p2d.size();
    if(use_point)
        J = MatrixXd::Zero(2 * n_pixels, 6 * n_poses + 3 * n_points);// add points
    else
        J = MatrixXd::Zero(2 * n_pixels, 6 * n_poses);
    r = VectorXd::Zero(2 * n_pixels, 1);
    MatrixXd J_pose_i(MatrixXd::Zero(2, 6));
    MatrixXd J_point_j(MatrixXd::Zero(2, 3));
    VectorXd r_pixels_ij(VectorXd::Zero(2, 1));
    MatrixXd J_cam(MatrixXd::Zero(2, 3));
    for(int i = 0; i < n_poses; i++)
    {
        for(int j = 0; j < n_points; j++)
        {
            Matrix4d Tcw = poses[i];
            Vector3d point = p3d[j];
            Vector2d pixel = p2d[i * n_points + j]; // assume all the points are observed by all the poses
            Matrix3d Rcw = Tcw.block(0, 0, 3, 3);
            Vector3d tcw = Tcw.block(0, 3, 3, 1);
            Vector3d p = Rcw*point + tcw;
            double x = p(0), y = p(1), z =p(2), z_2 = z * z;
            double p_u = cx + fx * x / z, p_v = cy + fy * y / z;
            double du = pixel(0) - p_u, dv = pixel(1) - p_v;
            e += (du * du + dv * dv);
            r_pixels_ij<< du, dv;
            J_cam << -fx/z, 0, fx*x/z_2,
                    0, -fy/z, fy*y/z_2;
            J_point_j = J_cam * Rcw;
            J_pose_i.block(0, 3, 2, 3) = - J_cam * skew(p);
            J_pose_i.block(0, 0, 2, 3) = J_cam;
            J.block<2, 6>(2*(i*n_points + j), 6*i) = J_pose_i;
            r.segment(2*(i*n_points + j), 2) = r_pixels_ij;
//            if(use_point)
//                J.block<2,3>(2*(i*n_points + j), 6*n_poses + 3*j) = J_point_j;
        }
    }
    e /= n_pixels;
}

void updatePosesAndPoints(VectorXd& delta_x)
{
    int n_poses = poses.size();
    int n_points = p3d.size();
    int n_pixels = p2d.size();
    for(int i = 0; i < n_poses; i++)
    {
        Vector6d delta_i = delta_x.segment(6*i, 6);
        Matrix4d T = poses[i];
        poses[i] = expmap_se3(delta_i) * T;
    }
    if(!use_point) return;
    for(int j = 0; j < n_points; j++)
    {
        Vector3d delta_j = delta_x.segment(6*n_poses + 3*j, 3);
        Vector3d p = p3d[j];
        p3d[j] = p + delta_j;
    }
}

void gaussNewtonOptimization(int max_iters)
{

    MatrixXd J;
    VectorXd r;
    MatrixXd H;
    MatrixXd g;
    double err = 0, err_prev = 999999999.9;

    // plot initial solution
    int iters;
    for( iters = 0; iters < max_iters; iters++)
    {
        optimizeFunctions( J, r, err );
        // if the difference is very small stop
        if( ( fabs(err-err_prev) < 1e-6 ) || ( err < 1e-6) )// || err > err_prev )
        {
            cout<<"Small reprojection error. Stop optimization."<<endl;
        }
        // update step
        H = J.transpose() * J;
        g = -J.transpose() * r;
        VectorXd delta_x = H.ldlt().solve(g);
        Vector6d delta = delta_x.head(6);
        //cout<<"Update vector: "<<delta_x.transpose()<<endl;
        if(delta.norm() < 1e-6)
        {
            cout<<"Small update vector. Stop optimization"<<endl;
            break;
        }
        // update previous values
        updatePosesAndPoints(delta_x);
        err_prev = err;
    }


}

int main(int argc, char *argv[]) {

    string p3d_file = string(argv[1]);
    string p2d_file = string(argv[2]);
    if (!readData(p3d_file, p2d_file)) {
        cerr << "Read data failed." << endl;
        return -1;
    }
    poses.push_back(Matrix4d::Identity());
    cout<<"BA problem contains "<<poses.size()<<" camera pose and "<< p3d.size()<<" points."<<endl;
    assert(poses.size() * p3d.size() == p2d.size());
    gaussNewtonOptimization(100);
    Eigen::Matrix4d T = poses[0];
    std::cout << "T = \n" << T << std::endl;
    return 0;
}
