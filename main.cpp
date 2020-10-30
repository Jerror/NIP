#include "NIP.hpp"
#include <iostream>

template <typename T>
T cubic(T x) { return x * x * x + x * x + x + 1; }

template <typename T>
T cubic_p(T x) { return 3 * x * x + 2 * x + 1; }

template <typename T>
T cubic_pp(T x) { return 6 * x + 2; }

template <typename T> void test_NIP(int maxiter, T tol) {
  NewtonIterativeProcedure<T> *NIP;
  RootResults<T> res;

  NIP = new NewtonIterativeProcedure<T>(cubic<T>, 0.5);
  res = NIP->solve(maxiter, tol);
  std::cout << res.root << " " << res.iterations << " " << res.function_calls
            << " " << res.flag << std::endl;

  NIP = new NewtonIterativeProcedure<T>(cubic<T>, cubic_p<T>, 0.5);
  res = NIP->solve(maxiter, tol);
  std::cout << res.root << " " << res.iterations << " " << res.function_calls
            << " " << res.flag << std::endl;

  NIP = new NewtonIterativeProcedure<T>(cubic<T>, cubic_p<T>, cubic_pp<T>, 0.5);
  res = NIP->solve(maxiter, tol);
  std::cout << res.root << " " << res.iterations << " " << res.function_calls
            << " " << res.flag << std::endl;
}

int main(int argc, char *argv[]) {
    std::cout<<"float"<<std::endl;
    test_NIP<float>(20, 1e-8);
    std::cout<<"double"<<std::endl;
    test_NIP<double>(40, 1e-16);
    std::cout<<"long double"<<std::endl;
    test_NIP<long double>(80, 1e-32);
}
