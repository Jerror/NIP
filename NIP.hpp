#ifndef NIP_H
#define NIP_H

#include <functional>

enum NIP_ROOT_FLAG { NIP_SUCCESS = 0,
    NIP_CONVERR = -1 };

template <typename T>
struct RootResults {
    T root;
    int iterations;
    int function_calls;
    NIP_ROOT_FLAG flag;
};

template <typename T>
struct IterResults {
    T x;
    T xprev;
    T est_err;
    int f_calls;
};

template <typename T> struct NIPIterGen {
    IterResults<T> res;
    NIPIterGen(T xprev, T x) : res{xprev, x, 1, 0} {}
    T operator()(){return 0;};
};

template <typename T>
class NewtonIterativeProcedure {

    NIPIterGen<T> method_iter_gen;
    int iterations = 0;

public:

    // Custom method constructor
    NewtonIterativeProcedure(NIPIterGen<T> my_method_iter_gen) : method_iter_gen{my_method_iter_gen} {}

    // Secant method constructors
    NewtonIterativeProcedure(std::function<T(T)> f, T x0, T x1);
    NewtonIterativeProcedure(std::function<T(T)> f, T x0);
    // Newton-Raphson method constructor
    NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime, T x0);
    // Halley's method constructor
    NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime, std::function<T(T)> fprime2, T x0);

    RootResults<T> solve(int maxiter = 50, T tol = 1e-8)
    {
        NIP_ROOT_FLAG flag = NIP_CONVERR;
        for (int i = 0; i < maxiter; ++i, ++iterations) {
            T est_err = method_iter_gen();
            if (est_err <= tol) {
                flag = NIP_SUCCESS;
                break;
            }
        }
        return RootResults<T> { method_iter_gen.res.x, iterations,
            method_iter_gen.res.f_calls, flag };
    }
};

#endif // NIP_H
