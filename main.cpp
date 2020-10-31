#include "NIP.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

#define co 2.99792458e8 // Speed of ligth in vacuum [m/s]

template <typename T>
T cubic(T x) { return x * x * x - x * x - x - 1; }

template <typename T>
T cubic_p(T x) { return 3 * x * x - 2 * x - 1; }

template <typename T>
T cubic_pp(T x) { return 6 * x - 2; }

template <typename T>
void print_res(RootResults<T>& res, int prec)
{
    std::cout << std::setprecision(prec) << res.root << " " << res.iterations << " " << res.function_calls
              << " " << res.flag << std::endl;
}

template <typename T>
void test_NIP(T x0, int maxiter, int prec)
{
    T tol = pow(10, -prec);
    NewtonIterativeProcedure<T>* NIP;
    RootResults<T> res;

    NIP = new NewtonIterativeProcedure<T>(cubic<T>, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);

    NIP = new NewtonIterativeProcedure<T>(cubic<T>, cubic_p<T>, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);

    NIP = new NewtonIterativeProcedure<T>(cubic<T>, cubic_p<T>, cubic_pp<T>, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);
}

// 3D dispersion (as function of wavevector magnitude) and derivatives

double disp3_k(double k, double w, double* khat, double* dxyzt)
{
    double ret = pow(sin(w * dxyzt[0] / 2) / dxyzt[0] / co, 2);
    for (int i = 0; i < 3; ++i) {
        ret -= pow(sin(k * khat[i] * dxyzt[i + 1] / 2) / dxyzt[i + 1], 2);
    }
    return ret;
}

double disp3_k_p(double k, double* khat, double* dxyzt)
{
    double ret = 0;
    for (int i = 0; i < 3; ++i) {
        ret -= khat[i] * sin(k * khat[i] * dxyzt[i + 1]) / dxyzt[i];
    }
    return ret / 2;
}

double disp3_k_pp(double k, double* khat, double* dxyzt)
{
    double ret = 0;
    for (int i = 0; i < 3; ++i) {
        ret -= khat[i] * khat[i] * cos(k * khat[i] * dxyzt[i + 1]);
    }
    return ret / 2;
}

// 1D dispersion (as function of grid spacing) and derivatives

double disp1_d(double d, double dt, double w, double k)
{
    return pow(sin(w * dt / 2) / dt / co, 2) - pow(sin(k * d / 2) / d, 2);
}

double disp1_d_p(double d, double k)
{
    return (2 * pow(sin(k * d / 2), 2) - 0.5 * d * k * sin(d * k)) / (d * d * d);
}

double disp1_d_pp(double d, double k)
{
    return ((6 - d * d * k * k) * cos(k * d) + 4 * d * k * sin(d * k) - 6) / (2 * d * d * d * d);
}

int main(int argc, char* argv[])
{
    std::cout << "float" << std::endl;
    test_NIP<float>(1.5, 20, 8);
    std::cout << "double" << std::endl;
    test_NIP<double>(1.5, 40, 16);
    std::cout << "long double" << std::endl;
    test_NIP<long double>(1.5, 80, 32);

    double th = M_PI;
    double ph = M_PI / 2;
    double dx = 60e-9;

    double w = 2 * M_PI * co / 4000e-9;
    double khat[3] = { sin(th) * cos(ph), sin(th) * sin(ph), cos(ph) };
    double dxyzt[4] = { 0.99 * sqrt(3) * dx / co, dx, dx, dx };

    using namespace std::placeholders;
    int prec = 16;
    double tol = pow(10, -prec);

    auto d3 = std::bind(disp3_k, _1, w, khat, dxyzt);
    auto d3_p = std::bind(disp3_k_p, _1, khat, dxyzt);
    auto d3_pp = std::bind(disp3_k_pp, _1, khat, dxyzt);
    NewtonIterativeProcedure<double> k_solver(d3, d3_p, d3_pp, w / co);
    // NewtonIterativeProcedure<double> k_solver(d3, d3_p, w / co);
    // NewtonIterativeProcedure<double> k_solver(d3, w / co);
    auto k_res = k_solver.solve(1000, tol);
    print_res(k_res, prec);
    if (k_res.flag != NIP_SUCCESS) {
        std::cerr << "ERROR: k failed to converge" << std::endl;
        abort();
    }
    double k = k_res.root;
    std::cout << "vp/c = " << k / w / co << std::endl;

    auto d1 = std::bind(disp1_d, _1, dxyzt[0], w, k);
    auto d1_p = std::bind(disp1_d_p, _1, k);
    auto d1_pp = std::bind(disp1_d_pp, _1, k);
    // NewtonIterativeProcedure<double> d_solver(d1, d1_p, d1_pp, dx);
    NewtonIterativeProcedure<double> d_solver(d1, dx);
    auto d_res = d_solver.solve(1000, tol);
    print_res(d_res, prec);
    std::cout << "d/dx = " << d_res.root / dx << std::endl;
}
