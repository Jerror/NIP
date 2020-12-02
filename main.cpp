#include "NIP.hpp"
#include "TFSF_MND.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <array>

// #define co 2.99792458e8 // Speed of ligth in vacuum [m/s]
#define co 1.0 // Naturalized
#define LOG10_2 0.301029995

template <typename T>
T cubic(T x) { return x * x * x - x * x - x - 1; }

template <typename T>
T cubic_p(T x) { return 3 * x * x - 2 * x - 1; }

template <typename T>
T cubic_pp(T x) { return 6 * x - 2; }

template <typename T>
void print_res(RootResults<T>& res, int prec)
{
    std::cout << std::setprecision(int(LOG10_2 * prec + 1)) << res.root << " " << res.iterations << " " << res.function_calls
              << " " << res.flag << std::setprecision(6) << std::endl;
}

template <typename T>
void test_NIP(std::function<T(T)> f,
    std::function<T(T)> fp,
    std::function<T(T)> fpp, T x0, int maxiter, int prec)
{
    T tol = pow(2, -prec);
    NewtonIterativeProcedure<T>* NIP;
    RootResults<T> res;

    NIP = new NewtonIterativeProcedure<T>(f, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);

    NIP = new NewtonIterativeProcedure<T>(f, fp, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);

    NIP = new NewtonIterativeProcedure<T>(f, fp, fpp, x0);
    res = NIP->solve(maxiter, tol);
    print_res(res, prec);
}

double findroot(std::function<double(double)> f,
    std::function<double(double)> fp,
    std::function<double(double)> fpp, double x0, int maxiter, int prec)
{
    double tol = pow(2, -prec);
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

void test_cubic() {
  std::cout << "Test: solve x^3 - x^2 - x - 1" << std::endl << std::endl;

  std::cout << "float" << std::endl;
  test_NIP<float>(cubic<float>, cubic_p<float>, cubic_pp<float>, 1.5, 20, 20);
  std::cout << "double" << std::endl;
  test_NIP<double>(cubic<double>, cubic_p<double>, cubic_pp<double>, 1.5, 40,
                   50);
  std::cout << "long double" << std::endl;
  test_NIP<long double>(cubic<long double>, cubic_p<long double>,
                        cubic_pp<long double>, 1.5, 80, 112);
}

void test_dispersion() {
    std::cout << std::endl
              << "Dispersion matching" << std::endl
              << std::endl;

    int maxiter = 50;
    double dx = 1;

    int n_N = 1000;
    int m_N = 1000;

    int failures = 0;

    for(int i = 3; i < 100; ++i)
    {
        double nlambda = i;
        for(int n = 0; n < n_N; ++n)
        {
            double th = M_PI*(n+0.5)/n_N;
            for(int m = 0; m < m_N; ++m)
            {
                double ph = M_PI*(m+0.5)/m_N;

                double w = 2 * M_PI * co / nlambda / dx;
                std::array<double, 3> khat{ sin(th) * cos(ph), sin(th) * sin(ph), cos(th) };
                std::array<double, 4> dxyzt{dx, dx, dx,
                    0.99 / sqrt(3) * dx / co};

                RootResults<double> kres, dres;
                int kprec, dprec;
                int x, y;

                DispersionWRTk3D d3(w, khat, dxyzt, co);
                for(kprec = 52; kprec > 0; --kprec)
                {
                  for (x = 0; x < 20; ++x) {
                    NewtonIterativeProcedure<double> ksolver(
                        d3, w / co * (20 - x) / 30);
                    // NewtonIterativeProcedure<double> ksolver(d3, d3.D(),
                    // d3.D2(), w / co * (20 - x)/30);
                    kres = ksolver.solve(maxiter, pow(2, -kprec));
                    if (kres.flag == NIP_SUCCESS)
                      goto k_done;
                  }
                }
            k_done:

                DispersionWRTd1D d1 = d3.Dispersion1D_to_match(kres.root);
                for(dprec = 52; dprec > 0; --dprec)
                {
                  for (y = 0; y < 20; ++y) {
                    NewtonIterativeProcedure<double> dsolver(d1, dx * (20 - y) /
                                                                     30);
                    // NewtonIterativeProcedure<double> dsolver(d1, d1.D(),
                    // d1.D2(), dx * (20 - y)/30);
                    dres = dsolver.solve(maxiter, pow(2, -dprec));
                    if (dres.flag == NIP_SUCCESS)
                      goto d_done;
                  }
                }
            d_done:

                if ((kprec < 52) || (dprec < 52))
                    std::cout<<kprec<<", "<<dprec<<std::endl;

                if ((kprec == 0) || (dprec == 0))
                    std::cerr<<++failures<<std::endl;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    // test_cubic();

    test_dispersion();

    // std::cout << std::endl
    //           << "Dispersion matching" << std::endl
    //           << std::endl;

    // double th, ph;
    // double nlambda;
    // double dx = 1;
    // int maxiter, prec;

    // int next_arg = 1;

    // th = M_PI / atof(argv[next_arg++]);
    // ph = M_PI / atof(argv[next_arg++]);
    // nlambda = atof(argv[next_arg++]);
    // if (argc > 6)
    //     dx = atof(argv[next_arg++]);
    // maxiter = atoi(argv[next_arg++]);
    // prec = atoi(argv[next_arg++]);

    // double w = 2 * M_PI * co / nlambda / dx;
    // std::array<double, 3> khat{ sin(th) * cos(ph), sin(th) * sin(ph), cos(th) };
    // std::array<double, 4> dxyzt{ dx, dx, dx, 0.99 / sqrt(3) * dx / co };

    // DispersionWRTk3D d3(w, khat, dxyzt, co);
    // double k = findroot(d3, d3.D(), d3.D2(), w / co, maxiter, prec);
    // std::cout << std::setprecision(LOG10_2 * prec + 1) << "vp/c = " << w / k / co << std::endl
    //           << std::endl
    //           << std::setprecision(6);

    // DispersionWRTd1D d1 = d3.Dispersion1D_to_match(k);
    // double d = findroot(d1, d1.D(), d1.D2(), dx, maxiter, prec);
    // std::cout << std::setprecision(LOG10_2 * prec + 1) << "d/dx = " << d / dx << std::endl
    //           << std::setprecision(6);

    // std::cout << std::endl
    //           << "Compare methods for dispersion matching" << std::endl
    //           << std::endl;
    // std::cout << "Solve 3D for k" << std::endl;
    // test_NIP<double>(d3, d3.D(), d3.D2(), w / co, maxiter, prec);
    // std::cout << std::endl
    //           << "Solve 1D for d" << std::endl;
    // test_NIP<double>(d1, d1.D(), d1.D2(), dx, maxiter, prec);
}
