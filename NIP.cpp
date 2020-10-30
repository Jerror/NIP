#include "NIP.hpp"

template <typename T>
struct SecantIterGen : NIPIterGen<T>{
    using NIPIterGen<T>::res;
    std::function<T(T)> f;

    SecantIterGen(T xprev, T x, std::function<T(T)> f_)
        : NIPIterGen<T> { xprev, x, 1, 0 }
        , f { f_ }
    {
    }

    T operator()()
    {
        T tmp = res.x - f(res.x) * (res.x - res.xprev) / (f(res.x) - f(res.xprev));
        res.xprev = res.x;
        res.x = tmp;
        res.est_err = abs((res.x - res.xprev) / res.x);
        res.f_calls += 3;
        return res.est_err;
    }
};

template <typename T>
struct NewtonRaphsonIterGen : NIPIterGen<T> {
    using NIPIterGen<T>::res;
    std::function<T(T)> f;
    std::function<T(T)> fprime;

    NewtonRaphsonIterGen(T x, std::function<T(T)> f_, std::function<T(T)> fprime_)
        : NIPIterGen<T> { x, x, 1, 0 }
        , f { f_ }
        , fprime { fprime_ }
    {
    }

    T operator()()
    {
        res.xprev = res.x;
        res.x = res.xprev - f(res.xprev) / fprime(res.xprev);
        res.est_err = abs((res.x - res.xprev) / res.x);
        res.f_calls += 2;
        return res.est_err;
    }
};

template <typename T>
struct HalleyIterGen : NIPIterGen<T> {
    using NIPIterGen<T>::res;
    std::function<T(T)> f;
    std::function<T(T)> fprime;
    std::function<T(T)> fprime2;

    HalleyIterGen(T x, std::function<T(T)> f_, std::function<T(T)> fprime_,
        std::function<T(T)> fprime2_)
        : NIPIterGen<T> { x, x, 1, 0 }
        , f { f_ }
        , fprime { fprime_ }
        , fprime2 { fprime2_ }
    {
    }

    T operator()()
    {
        res.xprev = res.x;
        T f_fp = f(res.xprev) / fprime(res.xprev);
        res.x = res.xprev - f_fp / (1 - f_fp * fprime2(res.xprev) / fprime(res.xprev) / 2);
        res.est_err = abs((res.x - res.xprev) / res.x);
        res.f_calls += 4;
        return res.est_err;
    }
};

template <typename T>
NewtonIterativeProcedure<T>::NewtonIterativeProcedure(std::function<T(T)> f, T x0, T x1)
    : method_iter_gen { SecantIterGen<T>(x0, x1, f) }
{
}

template <typename T>
NewtonIterativeProcedure<T>::NewtonIterativeProcedure(std::function<T(T)> f, T x0)
    : NewtonIterativeProcedure<T>(f, x0,
        x0 + (1 + x0) * (x0 >= 0 ? 1e-4 : -1e-4))
{
}

template <typename T>
NewtonIterativeProcedure<T>::NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime,
    T x0)
    : method_iter_gen { NewtonRaphsonIterGen<T>(x0, f, fprime) }
{
}

template <typename T>
NewtonIterativeProcedure<T>::NewtonIterativeProcedure(std::function<T(T)> f, std::function<T(T)> fprime,
    std::function<T(T)> fprime2, T x0)
    : method_iter_gen { HalleyIterGen<T>(x0, f, fprime, fprime2) }
{
}
