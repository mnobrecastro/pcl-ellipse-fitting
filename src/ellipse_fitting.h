/** 
 * Miguel Nobre Castro
 * mnobrecastro@gmail.com
 * Created on Apr 8th, 2021
 *
 * Implementation of "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)
 * and calculation of the minimal distance from an arbitrary point to the ellipse.
 */

#include <iostream>
#include <functional>
#include <thread>
#include <array>
#include <vector>
#include <string>
#include <cmath>

#include <Eigen/Eigenvalues>

#define M_PI    3.14159265358979323846   // pi
#define M_PHI   1.61803398874989484820   // golden_ratio

/*
 * Implementation of "Direct Least Square Fitting of Ellipses"
 * (Fitzgibbon et al., 1999)
 */
Eigen::VectorXf fit_ellipse(std::array<Eigen::Vector3d, 6> pts);

/*
 * Converts Conic parameters to Parametric parameters
 */
std::array<float, 5> conic2parametric(std::array<float, 6> con);

/*
 * Converts Parametric parameters to Conic parameters
 */
std::array<float, 6> parametric2conic(std::array<float, 5> par);

/*
 * Generates the points of an ellipse using its Parametric 'par' equation parameters, with a given arc 'ths' and noise.
 */
std::vector<std::array<float, 2>> ellipse_generator(std::array<float, 5> par, std::array<float, 2> arc = { 0,2 * 3.1415 }, int n = 1000, float noise = 0.0);

/*
 * Calculates a point on the ellipse model 'par' using the angle 'th'.
 */
void get_ellipse_point(std::array<float, 5> par, float th, float& x, float& y);

/*
 * Calculates a point (x,y) and its derivative (dx_dth,dy_dth) on the ellipse model 'par' using the angle 'th'.
 */
void get_ellipse_point(std::array<float, 5> par, float th, float& x, float& y, float& dx_dth, float& dy_dth);

/*
 * Minimum distance from point p=(u,v) to the ellipse model 'par'.
 */
float dist2ellipse(std::array<float, 5> par, float u, float v, float& th_opt);

/*
 * Distance between a point (u,v) to a given point in the ellipse model 'par' at an angle 'th'.
 */
float point2ellipse(std::array<float, 5> par, float th, float u, float v);

/*
 * Golden-section search
 */
float golden_section_search(std::array<float, 5>par, float u, float v, std::function<float(std::array<float, 5>, float, float, float)> fobj,
    float th_min, float th_max, float epsilon = 1e-6);