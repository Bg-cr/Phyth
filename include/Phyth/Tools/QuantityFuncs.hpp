#ifndef PHYTH_QUANTITY_FUNCS_HPP
#define PHYTH_QUANTITY_FUNCS_HPP

#include <algorithm>

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Quantities.hpp"
#include <cmath>

#include "Vector3.hpp"

namespace Phyth::Utils {
    inline Scalar sin(const Scalar q) {
        return {std::sin(q.value)};
    }

    inline Scalar cos(const Scalar q) {
        return {std::cos(q.value)};
    }

    inline Scalar tan(const Scalar q) {
        return {std::tan(q.value)};
    }

    inline Scalar asin(const Scalar q) {
        return {std::asin(q.value)};
    }

    inline Scalar acos(const Scalar q) {
        return {std::acos(q.value)};
    }

    inline Scalar atan(const Scalar q) {
        return {std::atan(q.value)};
    }

    inline Scalar atan2(const Scalar y, Scalar x) {
        return {std::atan2(y.value, x.value)};
    }

    inline Scalar sinh(const Scalar q) {
        return {std::sinh(q.value)};
    }

    inline Scalar cosh(const Scalar q) {
        return {std::cosh(q.value)};
    }

    inline Scalar tanh(const Scalar q) {
        return {std::tanh(q.value)};
    }

    inline Scalar asinh(const Scalar q) {
        return {std::asinh(q.value)};
    }

    inline Scalar acosh(const Scalar q) {
        return {std::acosh(q.value)};
    }

    inline Scalar atanh(const Scalar q) {
        return {std::atanh(q.value)};
    }


    inline Scalar log(const Scalar q) {
        return {std::log(q.value)};
    }

    inline Scalar log10(const Scalar q) {
        return {std::log10(q.value)};
    }

    inline Scalar log2(const Scalar q) {
        return {std::log2(q.value)};
    }

    inline Scalar log1p(const Scalar q) {
        return {std::log1p(q.value)};
    }

    inline Scalar exp(const Scalar q) {
        return {std::exp(q.value)};
    }

    inline Scalar expm1(const Scalar q) {
        return {std::expm1(q.value)};
    }

    inline Scalar exp2(const Scalar q) {
        return {std::exp2(q.value)};
    }

    template<int Power, typename UnitT>
    constexpr auto pow(Quantity<UnitT> q) {
        if constexpr (Power == 0) {
            return Scalar(1.0);
        } else if constexpr (Power > 0) {
            using DimT = typename UnitT::DimensionT;
            using ResultDim = DimPowerTypeT<Power, DimT>;
            using ResultUnit = Unit<ResultDim>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        } else {
            using DimT = typename UnitT::DimensionT;
            using ResultDim = DimRcpT<DimPowerTypeT<-Power, DimT>>;
            using ResultUnit = Unit<ResultDim>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        }
    }

    template<typename UnitT>
    constexpr auto square(Quantity<UnitT> q) {
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimMulT<DimT, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value);
    }

    template<typename UnitT>
    constexpr auto cube(Quantity<UnitT> q) {
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimMulT<DimMulT<DimT, DimT>, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value * q.value);
    }

    template<int N, typename UnitT>
    constexpr auto root(Quantity<UnitT> q) {
        static_assert(N > 0, "Root index must be positive");
        using DimT = typename UnitT::DimensionT;
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

    template<typename UnitT>
    constexpr Quantity<UnitT> fmod(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::fmod(x.value, y.value));
    }

    constexpr Scalar fmod(const Scalar x, const Scalar y) {
        return {std::fmod(x.value, y.value)};
    }

    template<typename UnitT>
    constexpr Quantity<UnitT> max(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::max(x.value, y.value));
    }

    template <typename UnitT>
    constexpr Quantity<UnitT> min(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::min(x.value, y.value));
    }

    template <typename UnitT>
    constexpr Quantity<UnitT> clamp(Quantity<UnitT> x, Quantity<UnitT> min, Quantity<UnitT> max) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::clamp(x.value, min.value, max.value));
    }
}

#endif //PHYTH_QUANTITY_FUNCS_HPP