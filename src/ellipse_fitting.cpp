/** 
 * Miguel Nobre Castro
 * mnobrecastro@gmail.com
 * Created on Apr 8th, 2021
 *
 * Implementation of "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)
 * and calculation of the minimal distance from an arbitrary point to the ellipse.
 */

#include "ellipse_fitting.h"

Eigen::VectorXf fit_ellipse(std::array<Eigen::Vector3d, 6> pts)
{
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
    std::cout << "* D matrix:\n" << D << '\n';

    // Scatter matrix S
    Eigen::MatrixXf S = Eigen::MatrixXf::Random(6, 6);
    S = D.transpose() * D;
    std::cout << "* S matrix:\n" << S << '\n';

    // Constraint matrix C
    Eigen::MatrixXf C = Eigen::MatrixXf::Zero(6, 6);
    C(0, 2) = -2.0;
    C(1, 1) = 1.0;
    C(2, 0) = -2.0;
    std::cout << "* C matrix:\n" << C << '\n';

    // Solve the Generalized Eigensystem: S*a = lambda*C*a
    /* ---- [ BUG FIX] ----
     * > Compilation error may occur from using 'GeneralizedEigenSolver' due to overloading of functions in PCL and Eigen libs.
     * + Solution: add the 'namespace pcl' to both 'aligned_{malloc,free}' in 'pcl/common/include/pcl/pcl_macros.h' (lines 379-415).
     * Please see: https://github.com/PointCloudLibrary/pcl/issues/4734#issuecomment-830115801
     */
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> solver;
    solver.compute(S, C);
    Eigen::VectorXf eigvals = solver.eigenvalues().real();
    std::cout << "Real generalized eigenvalues:\n" << eigvals.transpose() << '\n';
    std::cout << "Real generalized eigenvectors:\n" << solver.eigenvectors().real() << '\n';

    // Find the negative eigenvalue 'neigvec' (the largest, if many exist)
    int idx(0); float absmin(0.0);
    for (size_t i(0); i < eigvals.size(); ++i) {
        if (eigvals(i) < absmin && !std::isinf(eigvals(i)))
            idx = i;
    }
    
    Eigen::VectorXf neigvec = solver.eigenvectors().real().col(idx);

    return neigvec;
}


std::array<float, 5> conic2parametric(std::array<float,6> con)
{
    // Conic equation vars
    float con_A(0.0), con_B(0.0), con_C(0.0), con_D(0.0), con_E(0.0), con_F(0.0);
    con_A = con[0]; con_B = con[1]; con_C = con[2]; con_D = con[3]; con_E = con[4]; con_F = con[5];
    
    // Build matrix M0
    Eigen::MatrixXf M0(3,3);
    M0 << con_F, con_D / 2.0, con_E / 2.0,
        con_D / 2.0, con_A, con_B / 2.0,
        con_E / 2.0, con_B / 2.0, con_C;
    std::cout << "* M0 matrix:\n" << M0 << '\n';
    
    // Build matrix M
    Eigen::MatrixXf M(2,2);
    M << con_A, con_B/2.0,
        con_B/2.0, con_C;
    std::cout << "* M matrix:\n" << M << '\n';

    // Calculate the eigenvalues and eigenvectors of matrix M
    Eigen::EigenSolver<Eigen::MatrixXf> solver(M);
    std::cout << "Real eigenvalues of M:\n" << solver.eigenvalues().real().transpose() << '\n';
    //std::cout << "Real matrix of eigenvectors:\n" << solver.eigenvectors().real() << '\n' << '\n';
    Eigen::VectorXf eigvals = solver.eigenvalues().real();
    
    // Order the eigenvalues so that | lambda_0 - con_A| <= |lambda_0 - con_C |
    if (std::abs(eigvals(0) - con_A) > std::abs(eigvals(0) - con_C)) {
        float aux = eigvals(0);
        eigvals(0) = eigvals(1);
        eigvals(1) = aux;
    }
    std::cout << "Verify " << std::abs(eigvals(0) - con_A) << " <= " << std::abs(eigvals(0) - con_C) <<  '\n';
    std::cout << "So that " << std::abs(eigvals(1) - con_C) << " <= " << std::abs(eigvals(1) - con_A) << '\n';
    std::cout << "Real eigenvalues: " << eigvals.transpose() << '\n';
    
    // Parametric equation vars
    float par_a(0.0), par_b(0.0), par_h(0.0), par_k(0.0), par_t(0.0);
    par_a = std::sqrt( -M0.determinant() / (M.determinant() * eigvals(0)));
    par_b = std::sqrt( -M0.determinant() / (M.determinant() * eigvals(1)));
    par_h = (con_B * con_E - 2.0 * con_C * con_D) / (4.0 * con_A * con_C - std::pow(con_B, 2));
    par_k = (con_B * con_D - 2.0 * con_A * con_E) / (4.0 * con_A * con_C - std::pow(con_B, 2));
    par_t = (M_PI/2.0 - std::atan((con_A - con_C) / con_B)) / 2.0; // equivalent to acot((con_A - con_C) / con_B) / 2.0;

    return { par_a, par_b, par_h, par_k, par_t };
}


std::array<float,6> parametric2conic(std::array<float, 5> par)
{
    // Parametric equation vars
    float par_a(0.0), par_b(0.0), par_h(0.0), par_k(0.0), par_t(0.0);
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    // Conic equation vars
    float con_A, con_B, con_C, con_D, con_E, con_F;
    con_A = std::pow(par_b * std::cos(par_t), 2) + std::pow(par_a * std::sin(par_t), 2);
    con_B = -2 * std::cos(par_t) * std::sin(par_t) * (std::pow(par_a, 2) - std::pow(par_b, 2));
    con_C = std::pow(par_b * std::sin(par_t), 2) + std::pow(par_a * std::cos(par_t), 2);
    con_D = -2 * con_A * par_h - par_k * con_B;
    con_E = -2 * con_C * par_k - par_h * con_B;
    con_F = -std::pow(par_a * par_b, 2) + std::pow(con_A * par_h, 2) + con_B * par_h * par_k + con_C * std::pow(par_k, 2);

    return { con_A, con_B, con_C, con_D, con_E, con_F };
}


std::vector<std::array<float, 2>> ellipse_generator(std::array<float, 5> par, std::array<float, 2> arc, int n, float noise)
{
    // Parametric equation vars
    float par_a, par_b, par_h, par_k, par_t;
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    // Generate the angle steps
    std::vector<float> th(n, 0.0);
    for (std::size_t i(0); i < th.size(); ++i) { th[i] = arc[0] + i * (arc[1] - arc[0]) / (n); }

    // Calculate the Ellipse data points
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


void get_ellipse_point(std::array<float, 5> par, float th, float& x, float& y)
{
    // Parametric equation vars
    float par_a, par_b, par_h, par_k, par_t;
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    x = par_h + std::cos(par_t) * par_a * std::cos(th) - std::sin(par_t) * par_b * std::sin(th);
    y = par_k + std::sin(par_t) * par_a * std::cos(th) + std::cos(par_t) * par_b * std::sin(th);

    return;
}


void get_ellipse_point(std::array<float, 5> par, float th, float& x, float& y, float& dx_dth, float& dy_dth)
{
    // Parametric equation vars
    float par_a, par_b, par_h, par_k, par_t;
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    get_ellipse_point(par, th, x, y);

    dx_dth = -std::cos(par_t) * par_a * std::sin(th) - std::sin(par_t) * par_b * std::cos(th);
    dy_dth = -std::sin(par_t) * par_a * std::sin(th) + std::cos(par_t) * par_b * std::cos(th);

    return;
}


float dist2ellipse(std::array<float, 5> par, float u, float v, float& th_opt)
{
    // Parametric equation vars
    float par_a, par_b, par_h, par_k, par_t;
    par_a = par[0]; par_b = par[1]; par_h = par[2]; par_k = par[3]; par_t = par[4];

    Eigen::Vector2f center(par_h, par_k);
    Eigen::Vector2f p(u, v);
    p -= center;

    // Local x-axis of the ellipse
    Eigen::Vector2f x_axis;
    get_ellipse_point(par, 0.0, x_axis(0), x_axis(1));
    x_axis -= center;

    // Local y-axis of the ellipse
    Eigen::Vector2f y_axis;
    get_ellipse_point(par, M_PI/2.0, y_axis(0), y_axis(1));
    y_axis -= center;

    // Convert the point p=(u,v) to local ellipse coordinates
    float x_proj = p.dot(x_axis) / x_axis.norm();
    float y_proj = p.dot(y_axis) / y_axis.norm();

    // Find the ellipse quandrant to where the point p=(u,v) belongs,
    // and limit the search interval to 'th_min' and 'th_max'.
    float th_min(0.0), th_max(0.0);
    float th = std::atan2(y_proj, x_proj);
    std::cout << "Point (" << u << "," << v << "): ";

    if (-M_PI <= th && th < -M_PI/2.0) {
        th_min = -M_PI;
        th_max = -M_PI/2.0;
        // std::cout << "Q3 " << th * 180 / M_PI << '\n';
    }
    if (-M_PI/2.0 <= th && th < 0.0) {
        th_min = -M_PI/2.0;
        th_max = 0.0;
        // std::cout << "Q4 " << th * 180 / M_PI << '\n';
    }
    if (0.0 <= th && th < M_PI/2.0) {
        th_min = 0.0;
        th_max = M_PI/2.0;
        // std::cout << "Q1 " << th * 180 / M_PI << '\n';
    }
    if (M_PI/2.0 <= th && th <= M_PI) {
        th_min = M_PI/2.0;
        th_max = M_PI;
        // std::cout << "Q2 " << th * 180 / M_PI << '\n';
    }

    // Use an unconstrained line search optimizer
    th_opt = golden_section_search(par, u, v, &point2ellipse, th_min, th_max, 1.0e-3);
    float d = point2ellipse(par, th_opt, u, v);
    return d;
}


float point2ellipse(std::array<float, 5> par, float th, float u, float v)
{
    float x(0.0), y(0.0);
    get_ellipse_point(par, th, x, y);
    float d = std::sqrt(std::pow(u - x, 2) + std::pow(v - y, 2));
    return d;
}


float golden_section_search(std::array<float, 5>par, float u, float v, std::function<float(std::array<float, 5>, float, float, float)> fobj,
    float th_min, float th_max, float epsilon)
{
    float tl(th_min), tu(th_max), ta(0.0), tb(0.0);
    ta = tl + (tu - tl) * (1 - 1 / M_PHI);
    tb = tl + (tu - tl) * 1 / M_PHI;

    std::cout << "---- Golden-section search (" << u << "," << v << ") ----\n";
    std::cout << '\t' << (tl + tu) / 2.0 << '\t' << fobj(par, (tl + tu) / 2.0, u, v) << '\n';
    while ((tu - tl) > epsilon) {
        if (fobj(par, ta, u, v) < fobj(par, tb, u, v)) {
            tu = tb;
            tb = ta;
            ta = tl + (tu - tl) * (1.0 - 1.0 / M_PHI);
        }
        else if (fobj(par, ta, u, v) > fobj(par, tb, u, v)) {
            tl = ta;
            ta = tb;
            tb = tl + (tu - tl) * 1.0 / M_PHI;
        }
        else {
            tl = ta;
            tu = tb;
            ta = tl + (tu - tl) * (1.0 - 1.0 / M_PHI);
            tb = tl + (tu - tl) * 1.0 / M_PHI;
        }
        std::cout << '\t' << (tl + tu) / 2.0 << '\t' << fobj(par, (tl + tu) / 2.0, u, v) << '\n';
    }
    return (tl + tu) / 2.0;
}