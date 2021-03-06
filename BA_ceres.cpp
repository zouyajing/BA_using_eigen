#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <fstream>
#include "helper.h"

using namespace Eigen;
using namespace ceres;
using namespace std;

vector<Vector3d> p3d;
vector<Vector2d> p2d;


class ReprojectionErrorAnalytic : public SizedCostFunction<2,6>
{
public:
    ReprojectionErrorAnalytic(Vector3d point_, Vector2d observed_)
        : point(point_), observed(observed_) {
}
    virtual ~ReprojectionErrorAnalytic(){}
    virtual bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const
    {
        Map<const Vector6d> se3(parameters[0]);
        Matrix4d Tcw = expmap_se3(se3);
        Matrix3d Rcw = Tcw.block(0,0,3,3);
        Vector3d tcw = Tcw.block(0,3,3,1);
        Vector3d pt2 = Rcw*point + tcw;
        const double xp = double(K[0] * (pt2[0] / pt2[2]) + K[2]);
        const double yp = double(K[1] * (pt2[1] / pt2[2]) + K[3]);
        const double u = double(observed(0));
        const double v = double(observed(1));
        residuals[0] = u - xp;
        residuals[1] = v - yp;
        Eigen::Matrix<double, 2, 3, RowMajor> J_cam;
        J_cam << K[0]/pt2[2], 0, -K[0]*pt2[0]/(pt2[2]*pt2[2]),
                0, K[1]/pt2[2], -K[1]*pt2[1]/(pt2[2]*pt2[2]);
        if(jacobians != NULL)
        {
            if(jacobians[0] != NULL)
            {
                Map<Eigen::Matrix<double, 2, 6, RowMajor>> J_se3(jacobians[0]);
                J_se3.setZero();
                J_se3.block(0,3,2,3) = - J_cam * skew(pt2);
                J_se3.block(0,0,2,3) = J_cam;
            }
        }
        return true;
    }
private:
    Vector3d point;
    Vector2d observed;
    // Camera intrinsics
    double K[4] = {520.9, 521.0, 325.1, 249.7}; // fx,fy,cx,cy
};


class ReprojectionError6D {
public:

    ReprojectionError6D(Vector3d point_, Vector2d observed_)
            : point(point_), observed(observed_) {
    }
    bool operator()(const double *const camera, double *residuals) const {
          Vector6d se3 ;
          se3 << camera[0], camera[1], camera[2], camera[3], camera[4], camera[5];
          Matrix4d Tcw = expmap_se3(se3);
          Matrix3d Rcw = Tcw.block(0,0,3,3);
          Vector3d tcw = Tcw.block(0,3,3,1);
          Vector3d pt2 = Rcw*point + tcw;
          const double xp = double(K[0] * (pt2[0] / pt2[2]) + K[2]);
          const double yp = double(K[1] * (pt2[1] / pt2[2]) + K[3]);
          const double u = double(observed(0));
          const double v = double(observed(1));

          residuals[0] = u - xp;
          residuals[1] = v - yp;

        return true;
    }

    static ceres::CostFunction *Create(Vector3d points, Vector2d observed) {
        return (new ceres::NumericDiffCostFunction<ReprojectionError6D,FORWARD,2, 6>(
                new ReprojectionError6D(points, observed)));
    }

private:
    Vector3d point;
    Vector2d observed;
    // Camera intrinsics
    double K[4] = {520.9, 521.0, 325.1, 249.7}; // fx,fy,cx,cy
};






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

int main(int argc, char *argv[]) {

    // Google log
    google::InitGoogleLogging(argv[0]);
    string p3d_file = string(argv[1]);
    string p2d_file = string(argv[2]);

    if (!readData(p3d_file, p2d_file)) {
        cerr << "Read data failed." << endl;
        return -1;
    }

    ceres::Problem problem;
    ceres::LossFunction *lossfunction = NULL;
    double camera[6] = {0,0,0,0,0,0};

    for (uint i = 0; i < p3d.size(); i++) {
        ceres::CostFunction *costfunction = ReprojectionError6D::Create(p3d[i], p2d[i]);
        problem.AddResidualBlock(costfunction, lossfunction, camera);
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
    options.max_num_iterations = 100;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;

    Vector6d se3;
    se3 << camera[0], camera[1], camera[2], camera[3], camera[4], camera[5];
    Eigen::Matrix4d T = expmap_se3(se3);
    std::cout << "T = \n" << T << std::endl;
    return 0;
}
