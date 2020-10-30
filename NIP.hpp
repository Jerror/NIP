#ifndef NIP_H
#define NIP_H

#include <cmath> // fabs
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

template <typename T>
struct SecantIterGen {
    std::function<T(T)> f;

    SecantIterGen(std::function<T(T)> f_)
        : f { f_ }
    {
    }

    T operator()(IterResults<T>& state)
    {
        T tmp = state.x - f(state.x) * (state.x - state.xprev) / (f(state.x) - f(state.xprev));
        state.xprev = state.x;
        state.x = tmp;
        state.est_err = fabs((state.x - state.xprev) / state.x);
        state.f_calls += 3;
        return state.est_err;
    }
};

template <typename T>
struct NewtonRaphsonIterGen {
    std::function<T(T)> f;
    std::function<T(T)> fprime;

    NewtonRaphsonIterGen(std::function<T(T)> f_, std::function<T(T)> fprime_)
        : f { f_ }
        , fprime { fprime_ }
    {
    }

    T operator()(IterResults<T>& state)
    {
        state.xprev = state.x;
        state.x = state.xprev - f(state.xprev) / fprime(state.xprev);
        state.est_err = fabs((state.x - state.xprev) / state.x);
        state.f_calls += 2;
        return state.est_err;
    }
};

template <typename T>
struct HalleyIterGen {
    std::function<T(T)> f;
    std::function<T(T)> fprime;
    std::function<T(T)> fprime2;

    HalleyIterGen(std::function<T(T)> f_, std::function<T(T)> fprime_,
        std::function<T(T)> fprime2_)
        : f { f_ }
        , fprime { fprime_ }
        , fprime2 { fprime2_ }
    {
    }

    T operator()(IterResults<T>& state)
    {
        state.xprev = state.x;
        T f_fp = f(state.xprev) / fprime(state.xprev);
        state.x = state.xprev - f_fp / (1 - f_fp * fprime2(state.xprev) / fprime(state.xprev) / 2);
        state.est_err = fabs((state.x - state.xprev) / state.x);
        state.f_calls += 4;
        return state.est_err;
    }
};

template <typename T>
class NewtonIterativeProcedure {

    IterResults<T> state;
    std::function<void(IterResults<T>&)> method_iter_gen;
    int iterations = 0;

public:
    // Custom method constructor
    NewtonIterativeProcedure(IterResults<T> state_, std::function<void(IterResults<T>&)> method_iter_gen_)
        : state { state_ }
        , method_iter_gen { method_iter_gen_ }
    {
    }

    // Secant method constructors
    NewtonIterativeProcedure(std::function<T(T)> f, T x0, T x1)
        : state { x0, x1, 1, 0 }
        , method_iter_gen { SecantIterGen<T>(f) }
    {
    }
    NewtonIterativeProcedure(std::function<T(T)> f, T x0)
        : NewtonIterativeProcedure<T>(f, x0,
            x0 + (1 + x0) * (x0 >= 0 ? 1e-4 : -1e-4))
    {
    }

    // Newton-Raphson method constructor
    NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime,
        T x0)
        : state { x0, x0, 1, 0 }
        , method_iter_gen { NewtonRaphsonIterGen<T>(f, fprime) }
    {
    }

    // Halley's method constructor
    NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime, std::function<T(T)> fprime2, T x0)
        : state { x0, x0, 1, 0 }
        , method_iter_gen { HalleyIterGen<T>(f, fprime, fprime2) }
    {
    }

    RootResults<T> solve(int maxiter = 50, T tol = 1e-8)
    {
        NIP_ROOT_FLAG flag = NIP_CONVERR;
        for (int i = 0; i < maxiter; ++i, ++iterations) {
            method_iter_gen(state);
            if (state.est_err <= tol) {
                flag = NIP_SUCCESS;
                ++iterations;
                break;
            }
        }
        return RootResults<T> { state.x, iterations,
            state.f_calls, flag };
    }
};

#endif // NIP_H
