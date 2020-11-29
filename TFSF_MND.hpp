#ifndef TFSF_MND_H
#define TFSF_MND_H

#include "DebugUtil.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <iostream>

// 1D dispersion (as function of grid spacing) and derivatives
class DispersionWRTd1D {
    double dt; // temporal lattice constant
    double w; // wavelength
    double k; // numerical wavevector magnitude

    double co; // speed of light

    double f(double d)
    {
        return pow(sin(w * dt / 2) / dt / co, 2) - pow(sin(k * d / 2) / d, 2);
    }

    double df_dd(double d)
    {
        return (2 * pow(sin(k * d / 2), 2) - 0.5 * d * k * sin(d * k)) / (d * d * d);
    }

    double d2f_dd2(double d)
    {
        return ((6 - d * d * k * k) * cos(k * d) + 4 * d * k * sin(d * k) - 6) / (2 * d * d * d * d);
    }

public:
    DispersionWRTd1D(double dt_, double w_, double k_, double co_ = 1)
        // Default to natural units for numerical stability
        : dt { dt_ }
        , w { w_ }
        , k { k_ }
        , co { co_ }
    {
    }

    // Evaluate dispersion relation
    double operator()(double d)
    {
        return f(d);
    }

    // Return dispersion relation first d-derivative
    std::function<double(double)> D()
    {
        return std::bind(&DispersionWRTd1D::df_dd, this, std::placeholders::_1);
    }

    // Return dispersion relation second d-derivative
    std::function<double(double)> D2()
    {
        return std::bind(&DispersionWRTd1D::d2f_dd2, this, std::placeholders::_1);
    }
};

// 3D dispersion (as function of wavevector magnitude) and derivatives
class DispersionWRTk3D {
    double w; // wavelength
    std::array<double, 3> khat; // unit wavevector
    std::array<double, 4> dxyzt; // spacetime lattice constants

    double co; // speed of light

    double f(double k)
    {
        double ret = pow(sin(w * dxyzt[3] / 2) / dxyzt[3] / co, 2);
        for (int i = 0; i < 3; ++i) {
            ret -= pow(sin(k * khat[i] * dxyzt[i] / 2) / dxyzt[i], 2);
        }
        return ret;
    }

    double df_dk(double k)
    {
        double ret = 0;
        for (int i = 0; i < 3; ++i) {
            ret -= khat[i] * sin(k * khat[i] * dxyzt[i]) / dxyzt[i];
        }
        return ret / 2;
    }

    double d2f_dk2(double k)
    {
        double ret = 0;
        for (int i = 0; i < 3; ++i) {
            ret -= khat[i] * khat[i] * cos(k * khat[i] * dxyzt[i]);
        }
        return ret / 2;
    }

public:
    DispersionWRTk3D(double w_, std::array<double, 3> khat_,
        std::array<double, 4> dxyzt_, double co_ = 1.0)
        // Default to natural units for numerical stability
        : w { w_ }
        , khat { khat_ }
        , dxyzt { dxyzt_ }
        , co { co_ }
    {
        double knorm = khat[0]*khat[0] + khat[1]*khat[1] + khat[2]*khat[2];
        if (fabs(knorm - 1) > 1e-4) {
            WARNING("Normalized value of khat will be used in dispersion relation evaluation.")
        }
        for (int i = 0; i < 3; ++i) {
            khat[i] /= knorm;
        }
    }

    // Return 1D dispersion relation object which completes the MND system
    DispersionWRTd1D Dispersion1D_to_match(double k)
    {
        return DispersionWRTd1D(dxyzt[3], w, k, co);
    }

    // Evaluate dispersion relation
    double operator()(double k)
    {
        return f(k);
    }

    // Return dispersion relation first k-derivative
    std::function<double(double)> D()
    {
        return std::bind(&DispersionWRTk3D::df_dk, this, std::placeholders::_1);
    }

    // Return dispersion relation second k-derivative
    std::function<double(double)> D2()
    {
        return std::bind(&DispersionWRTk3D::d2f_dk2, this, std::placeholders::_1);
    }
};

#endif // TFSF_MND_H
