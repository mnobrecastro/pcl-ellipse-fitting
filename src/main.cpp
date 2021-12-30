/** 
 * Miguel Nobre Castro
 * mnobrecastro@gmail.com
 * Created on Apr 8th, 2021
 *
 * Implementation of "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)
 * and calculation of the minimal distance from an arbitrary point to the ellipse.
 */

#include <iostream>
#include "ellipse_fitting.h"

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/sample_consensus/sac_model_sphere.h>

using namespace std::chrono_literals;

int main(int argc, char** argv)
{
    // Define the parametric model of an Ellipse
    size_t model(2);
    std::array<float, 5> ellipse_par_model;
    std::vector<std::array<float, 2>> ellipse;
    int N(50); // n points

    switch (model) {
    case 0:
        ellipse_par_model = { 2.0, 1.0, 0.0, 0.0, 0.0 };
        ellipse = ellipse_generator(ellipse_par_model, { 0, 3.0*M_PI /2.0 }, N, 0.01);
        break;
    case 1:
        ellipse_par_model = { 2.0, 1.0, 0.0, 0.0, M_PI/4.0 }; //1/180*M_PI
        ellipse = ellipse_generator(ellipse_par_model, { 0, 3.0*M_PI/2.0 }, N, 0.01);
        break;
    case 2:
        ellipse_par_model = { 2.0, 1.0, 1.0, 1.0, M_PI/10.0 };
        ellipse = ellipse_generator(ellipse_par_model, { 0, M_PI/3.0 }, N, 0.01);
        break;
    default:
        return 1;
    }

    // PointCloud with points
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);
    cloud_in->width = N;
    cloud_in->height = 1;
    cloud_in->is_dense = false;
    cloud_in->points.resize(cloud_in->width * cloud_in->height);
    for (std::size_t i = 0; i < cloud_in->points.size(); ++i) {
        cloud_in->points[i].x = ellipse[i][0];
        cloud_in->points[i].y = ellipse[i][1];
        cloud_in->points[i].z = 0.0;
    }


    //************************************************************************************************
    std::cout << "A random number " << rand() % N << "/" << N << std::endl;

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

    // Dataset of six (6) points used for fitting the ellipse
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_data(new pcl::PointCloud<pcl::PointXYZ>);
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

    Eigen::VectorXf neigvec = fit_ellipse({p0,p1,p2,p3,p4,p5});
    std::cout << "The negative eigenvector: " << std::endl << neigvec << std::endl;

    // Convert the CONIC model to a PARAMETRIC model
    std::array<float, 6> con;
    for (size_t i(0); i < neigvec.size(); ++i) { con[i] = neigvec(i); }
    std::array<float, 5> par_model_out(conic2parametric(con));
 
    // Generate the output ellipse
    int N_out(1000);
    std::vector<std::array<float, 2>> ellipse_out = ellipse_generator(par_model_out, { 0, 2 * M_PI }, N_out, 0.0);
    
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZ>);
    cloud_out->width = N_out;
    cloud_out->height = 1;
    cloud_out->is_dense = false;
    cloud_out->points.resize(cloud_out->width * cloud_out->height);
    for (std::size_t i = 0; i < cloud_out->points.size(); ++i) {
        cloud_out->points[i].x = ellipse_out[i][0];
        cloud_out->points[i].y = ellipse_out[i][1];
        cloud_out->points[i].z = 0;
    }


    std::cout << std::endl;
    // Print the input Ellipse's Parametric parameters
    std::string parstr_in("");
    for (auto p : ellipse_par_model) { parstr_in += std::to_string(p) + " "; }
    std::cout << "The input PAR model: " << parstr_in << std::endl;
    // Print the output Ellipse's Parametric parameters
    std::string parstr_out("");
    for (auto p : par_model_out) { parstr_out += std::to_string(p) + " "; }
    std::cout << "The output PAR model: " << parstr_out << std::endl;


    //************************************************************************************************
    // Dataset of six (6) point used for fitting the ellipse
    std::array<float, 12> vec_th_opt({0.0});
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_dist(new pcl::PointCloud<pcl::PointXYZ>);
    cloud_dist->push_back(pcl::PointXYZ(2.5, 0.0, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(0.0, 2.5, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(-2.5, 0.0, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(0.0, -2.5, 0.0));

    cloud_dist->push_back(pcl::PointXYZ(0.5, 0.5, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(-0.5, 0.5, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(-0.5, -0.5, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(0.5, -0.5, 0.0));

    cloud_dist->push_back(pcl::PointXYZ(1.5, 1.0, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(-1.5, 1.0, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(-1.5, -1.0, 0.0));
    cloud_dist->push_back(pcl::PointXYZ(1.5, -1.0, 0.0));

    std::cout << std::endl;
    for (size_t i = 0; i < cloud_dist->points.size(); ++i) {
        std::cout << "Distancy (" << i << "): " << dist2ellipse(par_model_out, cloud_dist->points[i].x, cloud_dist->points[i].y, vec_th_opt[i]) << std::endl;
    }

    //************************************************************************************************

    // Initialise the Visualizer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
    int vp(0);
    viewer->createViewPort(0.0, 0.0, 1.0, 1.0, vp);
    viewer->setCameraPosition(0.0, 0.0, -10.0, 0.0, -1.0, 0.0, vp);
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
        //viewer->addSphere(*coefs, "set" + std::to_string(k));
        viewer->addSphere(p, 0.05, 1.0, 0.0, 1.0, "set" + std::to_string(k));// , int viewport)
    }

    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud_out_color_h(255, 255, 0);
    cloud_out_color_h.setInputCloud(cloud_out);
    viewer->addPointCloud(cloud_out, cloud_out_color_h, "cloud_out");

    // Dataset points
    std::vector<pcl::ModelCoefficients::Ptr> sphs_dist;
    size_t j(0);
    for (auto p : cloud_dist->points) {
        pcl::ModelCoefficients::Ptr coefs(new pcl::ModelCoefficients);
        coefs->values.push_back(p.x);
        coefs->values.push_back(p.y);
        coefs->values.push_back(p.z);
        coefs->values.push_back(0.05);
        sphs_dist.push_back(coefs);
        viewer->addSphere(*coefs, "pts" + std::to_string(j));

        // Draw the respective "shortest distance to ellipse" line
        float elp_x(0.0), elp_y(0.0), elp_z(0.0);
        get_ellipse_point(par_model_out, vec_th_opt[j], elp_x, elp_y);
        viewer->addLine(pcl::PointXYZ(p.x, p.y, p.z), pcl::PointXYZ(elp_x, elp_y, elp_z), "line" + std::to_string(j));

        ++j;
    }

    while (!viewer->wasStopped()) {
        viewer->spinOnce(1,true);
        std::this_thread::sleep_for(100ms);
    }
    return 0;
}