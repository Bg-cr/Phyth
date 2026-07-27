#ifndef PHYTH_QUANTITY_FUNCS_H
#define PHYTH_QUANTITY_FUNCS_H

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Quantities.hpp"
#include <cmath>

namespace Phyth::Utils {
    inline Scalar sin(Scalar q) {
        return Scalar(std::sin(q.value));
    }

    inline Scalar cos(Scalar q) {
        return Scalar(std::cos(q.value));
    }

    inline Scalar tan(Scalar q) {
        return Scalar(std::tan(q.value));
    }


    inline Scalar asin(Scalar q) {
        return Scalar(std::asin(q.value));
    }

    inline Scalar acos(Scalar q) {
        return Scalar(std::acos(q.value));
    }

    inline Scalar atan(Scalar q) {
        return Scalar(std::atan(q.value));
    }

    inline Scalar atan2(Scalar y, Scalar x) {
        return Scalar(std::atan2(y.value, x.value));
    }

    inline Scalar sinh(Scalar q) {
        return Scalar(std::sinh(q.value));
    }

    inline Scalar cosh(Scalar q) {
        return Scalar(std::cosh(q.value));
    }

    inline Scalar tanh(Scalar q) {
        return Scalar(std::tanh(q.value));
    }

    inline Scalar asinh(Scalar q) {
        return Scalar(std::asinh(q.value));
    }

    inline Scalar acosh(Scalar q) {
        return Scalar(std::acosh(q.value));
    }

    inline Scalar atanh(Scalar q) {
        return Scalar(std::atanh(q.value));
    }


    inline Scalar log(Scalar q) {
        return Scalar(std::log(q.value));
    }

    inline Scalar log10(Scalar q) {
        return Scalar(std::log10(q.value));
    }

    inline Scalar log2(Scalar q) {
        return Scalar(std::log2(q.value));
    }

    inline Scalar log1p(Scalar q) {
        return Scalar(std::log1p(q.value));
    }

    inline Scalar exp(Scalar q) {
        return Scalar(std::exp(q.value));
    }

    inline Scalar expm1(Scalar q) {
        return Scalar(std::expm1(q.value));
    }

    inline Scalar exp2(Scalar q) {
        return Scalar(std::exp2(q.value));
    }

    template<int Power, typename UnitT>
    constexpr auto pow(Quantity<UnitT> q) {
        if constexpr (Power == 0) {
            return Scalar(1.0);
        } else if constexpr (Power > 0) {
            using DimT = typename UnitT::Dimension;
            using ResultDim = DimPowerTypeT<Power, DimT>;
            using ResultUnit = Unit<ResultDim, std::ratio<1>>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        } else {
            using DimT = typename UnitT::Dimension;
            using ResultDim = DimNegT<DimPowerTypeT<-Power, DimT>>;
            using ResultUnit = Unit<ResultDim, std::ratio<1>>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        }
    }

    template<typename UnitT>
    constexpr auto square(Quantity<UnitT> q) {
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimMulT<DimT, DimT>;
        using ResultUnit = Unit<ResultDim, std::ratio<1>>;
        return Quantity<ResultUnit>(q.value * q.value);
    }

    template<typename UnitT>
    constexpr auto cube(Quantity<UnitT> q) {
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimMulT<DimMulT<DimT, DimT>, DimT>;
        using ResultUnit = Unit<ResultDim, std::ratio<1>>;
        return Quantity<ResultUnit>(q.value * q.value * q.value);
    }

    template<int N, typename UnitT>
    constexpr auto root(Quantity<UnitT> q) {
        static_assert(N > 0, "Root index must be positive");
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimRootT<N, DimT>;
        using ResultUnit = Unit<ResultDim, std::ratio<1>>;
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

    template<typename UnitT>
    constexpr Quantity<UnitT> fmod(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::fmod(x.value, y.value));
    }

    constexpr Scalar fmod(Scalar x, Scalar y) {
        return {std::fmod(x.value, y.value)};
    }

}

#endif //PHYTH_QUANTITY_FUNCS_H