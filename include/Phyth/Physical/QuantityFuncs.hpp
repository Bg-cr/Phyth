#ifndef PHYTH_QUANTITY_FUNCS_H
#define PHYTH_QUANTITY_FUNCS_H

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Quantities.hpp"
#include <cmath>

namespace Phyth::Utils {


    inline double sin(Scalar q) {
        return std::sin(q.value);
    }

    inline double cos(Scalar q) {
        return std::cos(q.value);
    }

    inline double tan(Scalar q) {
        return std::tan(q.value);
    }


    inline Scalar asin(double x) {
        return Scalar(std::asin(x));
    }

    inline Scalar acos(double x) {
        return Scalar(std::acos(x));
    }

    inline Scalar atan(double x) {
        return Scalar(std::atan(x));
    }

    inline Scalar atan2(double y, double x) {
        return Scalar(std::atan2(y, x));
    }


    inline double sinh(double x) { return std::sinh(x); }
    inline double cosh(double x) { return std::cosh(x); }
    inline double tanh(double x) { return std::tanh(x); }
    inline double asinh(double x) { return std::asinh(x); }
    inline double acosh(double x) { return std::acosh(x); }
    inline double atanh(double x) { return std::atanh(x); }


    inline double log(double x) { return std::log(x); }
    inline double log10(double x) { return std::log10(x); }
    inline double log2(double x) { return std::log2(x); }
    inline double log1p(double x) { return std::log1p(x); }
    inline double exp(double x) { return std::exp(x); }
    inline double expm1(double x) { return std::expm1(x); }
    inline double exp2(double x) { return std::exp2(x); }

    template<int Power, typename UnitT>
    constexpr auto pow(Quantity<UnitT> q) {
        if constexpr (Power == 0) {
            return Scalar(1.0);
        } else if constexpr (Power > 0) {
            using DimT = typename UnitT::Dimension;
            using ResultDim = DimPowerTypeT<Power, DimT>;
            using ResultUnit = Unit<ResultDim>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        } else {
            using DimT = typename UnitT::Dimension;
            using ResultDim = DimNegT<DimPowerTypeT<Power, DimT>>;
            using ResultUnit = Unit<ResultDim>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        }
    }

    template<typename UnitT>
    constexpr auto square(Quantity<UnitT> q) {
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimMulT<DimT, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value);
    }

    template<typename UnitT>
    constexpr auto cube(Quantity<UnitT> q) {
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimMulT<DimMulT<DimT, DimT>, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value * q.value);
    }


    template<int N, typename UnitT>
    constexpr auto root(Quantity<UnitT> q) {
        static_assert(N > 0, "Root index must be positive");
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimRootT<N, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(std::pow(q.value, 1.0 / N));
    }

    template<typename UnitT>
    constexpr auto sqrt(Quantity<UnitT> q) {
        return root<2>(q);
    }

    template<typename UnitT>
    constexpr auto cbrt(Quantity<UnitT> q) {
        return root<3>(q);
    }


    template<typename UnitT>
    constexpr auto abs(Quantity<UnitT> q) {
        return Quantity<UnitT>(std::abs(q.value));
    }

}

#endif  //PHYTH_QUANTITY_FUNCS_H