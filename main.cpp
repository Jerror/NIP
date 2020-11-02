#include "NIP.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

// #define co 2.99792458e8 // Speed of ligth in vacuum [m/s]
#define co 1.0 // Naturalized

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
              << " " << res.flag << std::setprecision(6) << std::endl;
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
    double ret = pow(sin(w * dxyzt[3] / 2) / dxyzt[3] / co, 2);
    for (int i = 0; i < 3; ++i) {
        ret -= pow(sin(k * khat[i] * dxyzt[i] / 2) / dxyzt[i], 2);
    }
    return ret;
}

double disp3_k_p(double k, double* khat, double* dxyzt)
{
    double ret = 0;
    for (int i = 0; i < 3; ++i) {
        ret -= khat[i] * sin(k * khat[i] * dxyzt[i]) / dxyzt[i];
    }
    return ret / 2;
}

double disp3_k_pp(double k, double* khat, double* dxyzt)
{
    double ret = 0;
    for (int i = 0; i < 3; ++i) {
        ret -= khat[i] * khat[i] * cos(k * khat[i] * dxyzt[i]);
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

double findroot(std::function<double(double)> f,
    std::function<double(double)> fp,
    std::function<double(double)> fpp, double x0, int maxiter, int prec)
{
    double tol = pow(10, -prec);
    RootResults<double> res;
    NewtonIterativeProcedure<double>* solver;

    solver = new NewtonIterativeProcedure<double>(f, fp, fpp, x0);
    res = solver->solve(maxiter, tol);
    print_res(res, prec);
    if (res.flag != NIP_SUCCESS) {
        std::cout << "Halley's method failed to converge, trying something else" << std::endl;
        solver = new NewtonIterativeProcedure<double>(f, fp, x0);
        res = solver->solve(maxiter, tol);
        print_res(res, prec);
        if (res.flag != NIP_SUCCESS) {
            std::cout << "Newton's method failed to converge, trying something else" << std::endl;
            solver = new NewtonIterativeProcedure<double>(f, x0);
            res = solver->solve(maxiter, tol);
            print_res(res, prec);
            if (res.flag != NIP_SUCCESS) {
                std::cerr << "Secant method failed to converge, giving up" << std::endl;
                abort();
            }
        }
    }
    return res.root;
}

int main(int argc, char* argv[])
{
    std::cout << "float" << std::endl;
    test_NIP<float>(1.5, 20, 8);
    std::cout << "double" << std::endl;
    test_NIP<double>(1.5, 40, 16);
    std::cout << "long double" << std::endl;
    test_NIP<long double>(1.5, 80, 32);

    std::cout << std::endl
              << "Dispersion matching" << std::endl
              << std::endl;

    double th, ph;
    double nlambda;
    double dx = 1;
    int maxiter, prec;

    int next_arg = 1;

    th = M_PI / atof(argv[next_arg++]);
    ph = M_PI / atof(argv[next_arg++]);
    nlambda = atof(argv[next_arg++]);
    if (argc > 6)
        dx = atof(argv[next_arg++]);
    maxiter = atoi(argv[next_arg++]);
    prec = atoi(argv[next_arg++]);

    double w = 2 * M_PI * co / nlambda / dx;
    double khat[3] = { sin(th) * cos(ph), sin(th) * sin(ph), cos(th) };
    double dxyzt[4] = { dx, dx, dx, 0.99 / sqrt(3) * dx / co };

    using namespace std::placeholders;
    auto d3 = std::bind(disp3_k, _1, w, khat, dxyzt);
    auto d3_p = std::bind(disp3_k_p, _1, khat, dxyzt);
    auto d3_pp = std::bind(disp3_k_pp, _1, khat, dxyzt);
    double k = findroot(d3, d3_p, d3_pp, w / co, maxiter, prec);
    std::cout << std::setprecision(prec) << "vp/c = " << w / k / co << std::endl
              << std::endl
              << std::setprecision(6);

    auto d1 = std::bind(disp1_d, _1, dxyzt[3], w, k);
    auto d1_p = std::bind(disp1_d_p, _1, k);
    auto d1_pp = std::bind(disp1_d_pp, _1, k);
    double d = findroot(d1, d1_p, d1_pp, dx, maxiter, prec);
    std::cout << std::setprecision(prec) << "d/dx = " << d / dx << std::endl
              << std::endl
              << std::setprecision(6);
}
