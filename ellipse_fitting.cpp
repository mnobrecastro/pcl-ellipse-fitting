/** 
 * Miguel Nobre Castro
 * mnobrecastro@gmail.com
 * Created on Apr 8th, 2021
 *
 * Implementation of "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)
 */

#include <iostream>
#include <thread>
#include <array>
#include <vector>
#include <string>
#include <cmath>

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/sample_consensus/sac_model_sphere.h>
#include <Eigen/Eigenvalues>

using namespace std::chrono_literals;

Eigen::VectorXf fit_ellipse(std::array<Eigen::Vector3d, 6>);
std::array<float, 5> conic2parametric(std::array<float, 6>);
std::array<float, 6> parametric2conic(std::array<float, 5>);
std::vector<std::array<float, 2>> ellipse_generator(std::array<float, 5>, std::array<float, 2>, int, float);

int main()//(int argc, char** argv)
{
    // Initialize the PointClouds
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_data(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZ>);

    int N(50);
    std::array<float, 5> ellipse_par_model = { 2.0, 1.0, 1.0, 1.0, M_PI/10.0 };
    std::vector<std::array<float, 2>> ellipse = ellipse_generator(ellipse_par_model, { 0, M_PI/3 }, N, 0.01); //0.01

    // PointCloud with points
    cloud_in->width = N;
    cloud_in->height = 1;
    cloud_in->is_dense = false;
    cloud_in->points.resize(cloud_in->width * cloud_in->height);
    for (std::size_t i = 0; i < cloud_in->points.size(); ++i)
    {
        cloud_in->points[i].x = ellipse[i][0];
        cloud_in->points[i].y = ellipse[i][1];
        cloud_in->points[i].z = 0.0;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "A random number " << rand() % N << "/" << N << std::endl;
    /*
    // Dataset of dummy ellipse points (testing purposes)
    Eigen::Vector3d p0(1.6823, 2.1258, 0.0); 
    Eigen::Vector3d p1(-0.2154, 1.5107, 0.0);
    Eigen::Vector3d p2(-0.8957, 0.3741, 0.0);
    Eigen::Vector3d p3(0.3129, -0.1252, 0.0);
    Eigen::Vector3d p4(2.2218, 0.4890, 0.0);
    Eigen::Vector3d p5(2.9022, 1.6246, 0.0);
    */
    // Dataset of points selected from the ellipse_generator()
    std::array<int, 6> ridx;
    for (size_t i(0); i < ridx.size(); ++i) {
        ridx[i] = int(rand()%N);
    }
    Eigen::Vector3d p0(cloud_in->points[ridx[0]].x, cloud_in->points[ridx[0]].y, cloud_in->points[ridx[0]].z);
    Eigen::Vector3d p1(cloud_in->points[ridx[1]].x, cloud_in->points[ridx[1]].y, cloud_in->points[ridx[1]].z);
    Eigen::Vector3d p2(cloud_in->points[ridx[2]].x, cloud_in->points[ridx[2]].y, cloud_in->points[ridx[2]].z);
    Eigen::Vector3d p3(cloud_in->points[ridx[3]].x, cloud_in->points[ridx[3]].y, cloud_in->points[ridx[3]].z);
    Eigen::Vector3d p4(cloud_in->points[ridx[4]].x, cloud_in->points[ridx[4]].y, cloud_in->points[ridx[4]].z);
    Eigen::Vector3d p5(cloud_in->points[ridx[5]].x, cloud_in->points[ridx[5]].y, cloud_in->points[ridx[5]].z);
    ////Eigen::Vector3d p6(cloud_in->points[600].x, cloud_in->points[600].y, cloud_in->points[600].z);

    // Dataset of six (6) point used for fitting the ellipse
    cloud_data->push_back(pcl::PointXYZ(p0(0), p0(1), p0(2)));
    cloud_data->push_back(pcl::PointXYZ(p1(0), p1(1), p1(2)));
    cloud_data->push_back(pcl::PointXYZ(p2(0), p2(1), p2(2)));
    cloud_data->push_back(pcl::PointXYZ(p3(0), p3(1), p3(2)));
    cloud_data->push_back(pcl::PointXYZ(p4(0), p4(1), p4(2)));
    cloud_data->push_back(pcl::PointXYZ(p5(0), p5(1), p5(2)));

    std::cout << "**** Points (global):\n"
        << p0.transpose() << std::endl
        << p1.transpose() << std::endl
        << p2.transpose() << std::endl
        << p3.transpose() << std::endl
        << p4.transpose() << std::endl
        << p5.transpose() << std::endl;

    //Eigen::Vector3d helper_vec01 = p0 - p1;
    //Eigen::Vector3d helper_vec02 = p0 - p2;
    //Eigen::Vector3d helper_vec10 = p1 - p0;
    //Eigen::Vector3d helper_vec12 = p1 - p2;
    //Eigen::Vector3d helper_vec20 = p2 - p0;
    //Eigen::Vector3d helper_vec21 = p2 - p1;

    //Eigen::Vector3d common_helper_vec = helper_vec10.cross(helper_vec21);
    //Eigen::Vector3d ellipse_normal = common_helper_vec.normalized();

    //// Coordinate transformation to a local reference frame (2D)
    //Eigen::Vector3d x_axis = (p1 - p0).normalized();
    //Eigen::Vector3d z_axis = ellipse_normal;
    //Eigen::Vector3d y_axis = z_axis.cross(x_axis).normalized();

    //// Create the transposed rotation matrix
    //Eigen::Matrix3d Rot;
    //Rot << x_axis(0), x_axis(1), x_axis(2),
    //    y_axis(0), y_axis(1), y_axis(2),
    //    z_axis(0), z_axis(1), z_axis(2);
    //std::cout << "**** Rot^T matrix:\n" << Rot << std::endl;

    //// Convert the points to local coordinates
    //p0 = Rot * (p0 - p6);
    //p1 = Rot * (p1 - p6);
    //p2 = Rot * (p2 - p6);
    //p3 = Rot * (p3 - p6);
    //p4 = Rot * (p4 - p6);
    //p5 = Rot * (p5 - p6);

    //std::cout << "**** Points (local):\n"
    //    << p0.transpose() << std::endl
    //    << p1.transpose() << std::endl
    //    << p2.transpose() << std::endl
    //    << p3.transpose() << std::endl
    //    << p4.transpose() << std::endl
    //    << p5.transpose() << std::endl;


    Eigen::VectorXf neigvec = fit_ellipse({p0,p1,p2,p3,p4,p5});
    std::cout << "The negative eigenvector: " << std::endl << neigvec << std::endl;

    // Convert the CONIC model to a PARAMETRIC model
    std::array<float, 6> con;
    for (size_t i(0); i < neigvec.size(); ++i) { con[i] = neigvec(i); }
    std::array<float, 5> par_model_out(conic2parametric(con));
 
    // Generate the output ellipse
    int N_out(1000);
    std::vector<std::array<float, 2>> ellipse_out = ellipse_generator(par_model_out, { 0, 2 * M_PI }, N_out, 0.0);
    cloud_out->width = N_out;
    cloud_out->height = 1;
    cloud_out->is_dense = false;
    cloud_out->points.resize(cloud_out->width * cloud_out->height);
    for (std::size_t i = 0; i < cloud_out->points.size(); ++i)
    {
        cloud_out->points[i].x = ellipse_out[i][0];
        cloud_out->points[i].y = ellipse_out[i][1];
        cloud_out->points[i].z = 0;
    }


    std::cout << std::endl;
    // Print the inpu Ellipse's Parametric parameters
    std::string parstr_in("");
    for (auto p : ellipse_par_model) { parstr_in += std::to_string(p) + " "; }
    std::cout << "The input PAR model: " << parstr_in << std::endl;
    // Print the output Ellipse's Parametric parameters
    std::string parstr_out("");
    for (auto p : par_model_out) { parstr_out += std::to_string(p) + " "; }
    std::cout << "The output PAR model: " << parstr_out << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////  


    // Initialize the Visualizer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    int vp(0);
    viewer->createViewPort(0.0, 0.0, 1.0, 1.0, vp);
    viewer->setCameraPosition(0.0, 0.0, -0.5, 0.0, -1.0, 0.0, vp);
    viewer->setSize(800, 600);
    float bckgr_gray_level = 0.0;  // Black:=0.0
    viewer->setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, vp);
    viewer->addCoordinateSystem(1.0);

    viewer->removeAllShapes();
    viewer->removeAllPointClouds();

    // Input ellipse cloud
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_in_color_h(255, 255, 255);
    cloud_in_color_h.setInputCloud(cloud_in);
    viewer->addPointCloud(cloud_in, cloud_in_color_h, "cloud");

    // Dataset cloud
    std::vector<pcl::ModelCoefficients::Ptr> sphs;
    size_t k(0);
    for (auto p : cloud_data->points) {
        pcl::ModelCoefficients::Ptr coefs(new pcl::ModelCoefficients);
        coefs->values.push_back(p.x);
        coefs->values.push_back(p.y);
        coefs->values.push_back(p.z);
        coefs->values.push_back(0.05);
        sphs.push_back(coefs);
        ++k;
        viewer->addSphere(*coefs, std::to_string(k));
    }

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_out_color_h(255, 255, 0);
    cloud_out_color_h.setInputCloud(cloud_out);
    viewer->addPointCloud(cloud_out, cloud_out_color_h, "cloud_out");

    while (!viewer->wasStopped())
    {
        viewer->spinOnce(1,true);
        std::this_thread::sleep_for(100ms);
    }
    return 0;
}




Eigen::VectorXf fit_ellipse(std::array<Eigen::Vector3d, 6> pts)
{
    /*
     * 
     */

    // 2D projections only
    Eigen::VectorXf X(6);
    X << pts[0](0), pts[1](0), pts[2](0), pts[3](0), pts[4](0), pts[5](0);
    Eigen::VectorXf Y(6);
    Y << pts[0](1), pts[1](1), pts[2](1), pts[3](1), pts[4](1), pts[5](1);

    // Design matrix D  
    Eigen::MatrixXf D(6, 6);
    D << std::pow(X(0), 2), X(0)* Y(0), std::pow(Y(0), 2), X(0), Y(0), 1.0,
        std::pow(X(1), 2), X(1)* Y(1), std::pow(Y(1), 2), X(1), Y(1), 1.0,
        std::pow(X(2), 2), X(2)* Y(2), std::pow(Y(2), 2), X(2), Y(2), 1.0,
        std::pow(X(3), 2), X(3)* Y(3), std::pow(Y(3), 2), X(3), Y(3), 1.0,
        std::pow(X(4), 2), X(4)* Y(4), std::pow(Y(4), 2), X(4), Y(4), 1.0,
        std::pow(X(5), 2), X(5)* Y(5), std::pow(Y(5), 2), X(5), Y(5), 1.0;
    std::cout << "**** D matrix:\n" << D << std::endl;

    // Scatter matrix S
    Eigen::MatrixXf S = Eigen::MatrixXf::Random(6, 6);
    S = D.transpose() * D;
    std::cout << "**** S matrix:\n" << S << std::endl;

    // Constraint matrix C
    Eigen::MatrixXf C = Eigen::MatrixXf::Zero(6, 6);
    C(0, 2) = -2.0;
    C(1, 1) = 1.0;
    C(2, 0) = -2.0;
    std::cout << "**** C matrix:\n" << C << std::endl;

    /* ---- [ BUG FIX] ----
     * > Compilation error may occur from using 'GeneralizedEigenSolver' due to overloading of functions in PCL and Eigen libs.
     * + Solution: add the 'namespace pcl' to both 'aligned_{malloc,free}' in 'pcl/common/include/pcl/pcl_macros.h' (lines 379-415).
     * This requires rebuilding the PCL lib from the source. More info can be found in the link below:
     * https://github.com/PointCloudLibrary/pcl/issues/4734#issuecomment-830115801
     */
     // Solve the Generalized Eigensystem: S*a = lambda*C*a
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> solver;
    solver.compute(S, C);
    Eigen::VectorXf eigvals = solver.eigenvalues().real();
    std::cout << "Real generalized eigenvalues: " << std::endl << eigvals.transpose() << std::endl;
    std::cout << "Real generalized eigenvectors: " << std::endl << solver.eigenvectors().real() << std::endl;

    // Find the negative eigenvalue 'neigvec' (the largest, if many exist)
    int idx(0); float absmin(0.0);
    for (size_t i(0); i < eigvals.size(); ++i) {
        if (eigvals(i) < absmin && !std::isinf(eigvals(i))) { idx = i; }// break;
    }
    
    Eigen::VectorXf neigvec = solver.eigenvectors().real().col(idx);// .normalized(); 
    /* Matlab normalization (by the largest absolute value)
    //if (std::abs(neigvec.maxCoeff()) > std::abs(neigvec.minCoeff())) {
    //    neigvec /= std::abs(neigvec.maxCoeff());
    //} else {
    //    neigvec /= std::abs(neigvec.minCoeff());
    //}
    */

    return neigvec;
}


std::array<float, 5> conic2parametric(std::array<float,6> con)
{
    /**
     * Converts Conic parameters to Parametric parameters
     */

    // Conic params
    float con_A(0.0), con_B(0.0), con_C(0.0), con_D(0.0), con_E(0.0), con_F(0.0);
    con_A = con[0]; con_B = con[1]; con_C = con[2]; con_D = con[3]; con_E = con[4]; con_F = con[5];
    
    // Build matrix M0
    Eigen::MatrixXf M0(3,3);
    M0 << con_F, con_D / 2.0, con_E / 2.0,
        con_D / 2.0, con_A, con_B / 2.0,
        con_E / 2.0, con_B / 2.0, con_C;
    std::cout << "**** M0 matrix:\n" << M0 << std::endl;
    
    // Build matrix M
    Eigen::MatrixXf M(2,2);
    M << con_A, con_B/2.0,
        con_B/2.0, con_C;
    std::cout << "**** M matrix:\n" << M << std::endl;

    // Calculate the eigenvalues and eigenvectors of matrix M
    Eigen::EigenSolver<Eigen::MatrixXf> solver(M);
    std::cout << "Real eigenvalues of M are:" << endl << solver.eigenvalues().real().transpose() << std::endl;
    //std::cout << "Real matrix of eigenvectors, V, is:" << endl << solver.eigenvectors().real() << std::endl << std::endl;
    Eigen::VectorXf eigvals = solver.eigenvalues().real();
    
    // Order the eigenvalues so that | lambda_0 - con_A| <= |lambda_0 - con_C |
    if (std::abs(eigvals(0) - con_A) > std::abs(eigvals(0) - con_C)) {
        float aux = eigvals(0);
        eigvals(0) = eigvals(1);
        eigvals(1) = aux;
    }
    std::cout << "Verify " << std::abs(eigvals(0) - con_A) << " <= " << std::abs(eigvals(0) - con_C) <<  std::endl;
    std::cout << "So that " << std::abs(eigvals(1) - con_C) << " <= " << std::abs(eigvals(1) - con_A) << std::endl;
    std::cout << "Real eigenvalues: " << eigvals.transpose() << std::endl;
    
    // Parametric eq.params
    float par_a(0.0), par_b(0.0), par_h(0.0), par_k(0.0), par_t(0.0);
    par_a = std::sqrt( -M0.determinant() / (M.determinant() * eigvals(0)));
    par_b = std::sqrt( -M0.determinant() / (M.determinant() * eigvals(1)));
    par_h = (con_B * con_E - 2.0 * con_C * con_D) / (4.0 * con_A * con_C - std::pow(con_B, 2));
    par_k = (con_B * con_D - 2.0 * con_A * con_E) / (4.0 * con_A * con_C - std::pow(con_B, 2));
    par_t = (M_PI/2.0 - std::atan((con_A - con_C) / con_B)) / 2.0; // equivalent to: acot((con_A - con_C) / con_B) / 2.0;

    return { par_a, par_b, par_h, par_k, par_t };
}

std::array<float,6> parametric2conic(std::array<float, 5> par)
{
    /**
     * Converts Parametric parameters to Conic parameters
     */
     
    // Parametric eq.params
    float par_a(0.0), par_b(0.0), par_h(0.0), par_k(0.0), par_t(0.0);
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    // Conic eq.params
    float con_A, con_B, con_C, con_D, con_E, con_F;
    con_A = std::pow(par_b * std::cos(par_t), 2) + std::pow(par_a * std::sin(par_t), 2);
    con_B = -2 * std::cos(par_t) * std::sin(par_t) * (std::pow(par_a, 2) - std::pow(par_b, 2));
    con_C = std::pow(par_b * std::sin(par_t), 2) + std::pow(par_a * std::cos(par_t), 2);
    con_D = -2 * con_A * par_h - par_k * con_B;
    con_E = -2 * con_C * par_k - par_h * con_B;
    con_F = -std::pow(par_a * par_b, 2) + std::pow(con_A * par_h, 2) + con_B * par_h * par_k + con_C * std::pow(par_k, 2);

    return { con_A, con_B, con_C, con_D, con_E, con_F };
}

std::vector<std::array<float, 2>> ellipse_generator(std::array<float, 5> par, std::array<float, 2> arc = { 0,2 * 3.1415 }, int n = 1000, float noise = 0.0)
{
    /**
     * Generates the points of an ellipse using its Parametric 'par'
     * equation parameters, with a given arc 'ths' and noise.
     */

     // Parametric eq.params
    float par_a, par_b, par_h, par_k, par_t;
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    // Generate the angle steps
    std::vector<float> th(n, 0.0);
    for (std::size_t i(0); i < th.size(); ++i) { th[i] = arc[0] + i * (arc[1] - arc[0]) / (n); }

    // Calculate Ellipse data points
    std::vector<std::array<float, 2>> pts(th.size(), { 0.0,0.0 });
    for (std::size_t i(0); i < th.size(); ++i) {
        pts[i][0] = par_h + std::cos(par_t) * par_a * std::cos(th[i]) - std::sin(par_t) * par_b * std::sin(th[i]);
        pts[i][1] = par_k + std::sin(par_t) * par_a * std::cos(th[i]) + std::cos(par_t) * par_b * std::sin(th[i]);

        // Introduce noise in the data
        pts[i][0] += std::rand() % 100 / 100.0 * 2 * noise - noise;
        pts[i][1] += std::rand() % 100 / 100.0 * 2 * noise - noise;
    }

    return pts;
}