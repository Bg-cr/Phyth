#ifndef PHYTH_QUANTITY_FUNCS_HPP
#define PHYTH_QUANTITY_FUNCS_HPP

#include <algorithm>
#include <cmath>

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Quantities.hpp"

namespace Phyth::Utils {
    /**
     * @brief Sine of a dimensionless quantity
     *
     * @param q Angle in radians (dimensionless)
     * @return sin(q) as a dimensionless scalar
     */
    inline Scalar sin(const Scalar q) {
        return {std::sin(q.value)};
    }

    /**
     * @brief Cosine of a dimensionless quantity
     *
     * @param q Angle in radians (dimensionless)
     * @return cos(q) as a dimensionless scalar
     */
    inline Scalar cos(const Scalar q) {
        return {std::cos(q.value)};
    }

    /**
     * @brief Tangent of a dimensionless quantity
     *
     * @param q Angle in radians (dimensionless)
     * @return tan(q) as a dimensionless scalar
     */
    inline Scalar tan(const Scalar q) {
        return {std::tan(q.value)};
    }

    /**
     * @brief Arc sine (inverse sine)
     *
     * @param q Value in [-1, 1] (dimensionless)
     * @return Angle in radians (dimensionless)
     */
    inline Scalar asin(const Scalar q) {
        return {std::asin(q.value)};
    }

    /**
     * @brief Arc cosine (inverse cosine)
     *
     * @param q Value in [-1, 1] (dimensionless)
     * @return Angle in radians (dimensionless)
     */
    inline Scalar acos(const Scalar q) {
        return {std::acos(q.value)};
    }

    /**
     * @brief Arc tangent (inverse tangent)
     *
     * @param q Value (dimensionless)
     * @return Angle in radians (dimensionless)
     */
    inline Scalar atan(const Scalar q) {
        return {std::atan(q.value)};
    }

    /**
     * @brief Two-argument arc tangent (atan2)
     *
     * @param y Y component (dimensionless)
     * @param x X component (dimensionless)
     * @return Angle in radians (dimensionless), in [-pi, pi]
     */
    inline Scalar atan2(const Scalar y, const Scalar x) {
        return {std::atan2(y.value, x.value)};
    }

    /**
     * @brief Hyperbolic sine
     *
     * @param q Value (dimensionless)
     * @return sinh(q) as a dimensionless scalar
     */
    inline Scalar sinh(const Scalar q) {
        return {std::sinh(q.value)};
    }

    /**
     * @brief Hyperbolic cosine
     *
     * @param q Value (dimensionless)
     * @return cosh(q) as a dimensionless scalar
     */
    inline Scalar cosh(const Scalar q) {
        return {std::cosh(q.value)};
    }

    /**
     * @brief Hyperbolic tangent
     *
     * @param q Value (dimensionless)
     * @return tanh(q) as a dimensionless scalar
     */
    inline Scalar tanh(const Scalar q) {
        return {std::tanh(q.value)};
    }

    /**
     * @brief Inverse hyperbolic sine
     *
     * @param q Value (dimensionless)
     * @return asinh(q) as a dimensionless scalar
     */
    inline Scalar asinh(const Scalar q) {
        return {std::asinh(q.value)};
    }

    /**
     * @brief Inverse hyperbolic cosine
     *
     * @param q Value (dimensionless)
     * @return acosh(q) as a dimensionless scalar
     */
    inline Scalar acosh(const Scalar q) {
        return {std::acosh(q.value)};
    }

    /**
     * @brief Inverse hyperbolic tangent
     *
     * @param q Value in (-1, 1) (dimensionless)
     * @return atanh(q) as a dimensionless scalar
     */
    inline Scalar atanh(const Scalar q) {
        return {std::atanh(q.value)};
    }

    /**
     * @brief Natural logarithm (base e)
     *
     * @param q Positive value (dimensionless)
     * @return ln(q) as a dimensionless scalar
     */
    inline Scalar log(const Scalar q) {
        return {std::log(q.value)};
    }

    /**
     * @brief Base-10 logarithm
     *
     * @param q Positive value (dimensionless)
     * @return log10(q) as a dimensionless scalar
     */
    inline Scalar log10(const Scalar q) {
        return {std::log10(q.value)};
    }

    /**
     * @brief Base-2 logarithm
     *
     * @param q Positive value (dimensionless)
     * @return log2(q) as a dimensionless scalar
     */
    inline Scalar log2(const Scalar q) {
        return {std::log2(q.value)};
    }

    /**
     * @brief log(1 + q) with higher precision for small q
     *
     * @param q Value (dimensionless)
     * @return log(1+q) as a dimensionless scalar
     */
    inline Scalar log1p(const Scalar q) {
        return {std::log1p(q.value)};
    }

    /**
     * @brief Exponential (e^q)
     *
     * @param q Value (dimensionless)
     * @return e^q as a dimensionless scalar
     */
    inline Scalar exp(const Scalar q) {
        return {std::exp(q.value)};
    }

    /**
     * @brief exp(q) - 1 with higher precision for small q
     *
     * @param q Value (dimensionless)
     * @return e^q - 1 as a dimensionless scalar
     */
    inline Scalar expm1(const Scalar q) {
        return {std::expm1(q.value)};
    }

    /**
     * @brief Base-2 exponential (2^q)
     *
     * @param q Value (dimensionless)
     * @return 2^q as a dimensionless scalar
     */
    inline Scalar exp2(const Scalar q) {
        return {std::exp2(q.value)};
    }

    /**
     * @brief Raise a quantity to an integer power
     *
     * The resulting dimension is the original dimension raised to the power.
     *
     * @tparam Power Integer exponent (can be negative)
     * @tparam UnitT Unit type of the input quantity
     * @param q The quantity to raise
     * @return Quantity with dimension DimT^Power
     *
     * Example:
     *   auto t = 2.0_s;
     *   auto t2 = Utils::pow<2>(t);   // 4 s^2
     *   auto inv_t = Utils::pow<-1>(t); // 0.5 s^-1
     */
    template<int Power, typename UnitT>
    auto pow(Quantity<UnitT> q) {
        if (Power == 0) {
            return Scalar(1.0);
        }
        if (Power > 0) {
            using DimT = typename UnitT::DimensionT;
            using ResultDim = DimPowerTypeT<Power, DimT>;
            using ResultUnit = Unit<ResultDim>;
            return Quantity<ResultUnit>(std::pow(q.value, Power));
        }
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimRcpT<DimPowerTypeT<-Power, DimT> >;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(std::pow(q.value, Power));
    }

    /**
     * @brief Square of a quantity
     *
     * @tparam UnitT Unit type of the input quantity
     * @param q The quantity to square
     * @return Quantity with dimension DimT^2
     *
     * Example:
     *   auto t = 2.0_s;
     *   auto t2 = Utils::square(t);  // 4 s^2
     */
    template<typename UnitT>
    auto square(Quantity<UnitT> q) {
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimMulT<DimT, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value);
    }

    /**
     * @brief Cube of a quantity
     *
     * @tparam UnitT Unit type of the input quantity
     * @param q The quantity to cube
     * @return Quantity with dimension DimT^3
     */
    template<typename UnitT>
    auto cube(Quantity<UnitT> q) {
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimMulT<DimMulT<DimT, DimT>, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(q.value * q.value * q.value);
    }

    /**
     * @brief N-th root of a quantity
     *
     * @tparam N Root index (positive integer)
     * @tparam UnitT Unit type of the input quantity
     * @param q The quantity to take the root of
     * @return Quantity with dimension DimT^(1/N)
     *
     * Example:
     *   auto area = 9.0_m2;
     *   auto length = Utils::root<2>(area);  // 3.0_m
     */
    template<int N, typename UnitT>
    auto root(Quantity<UnitT> q) {
        static_assert(N > 0, "Root index must be positive");
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimRootT<N, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(std::pow(q.value, 1.0 / N));
    }

    /**
     * @brief Square root of a quantity
     *
     * @param q The quantity to take the square root of
     * @return Quantity with dimension DimT^(1/2)
     *
     * Example:
     *   auto area = 9.0_m2;
     *   auto length = Utils::sqrt(area);  // 3.0_m
     */
    template<typename UnitT>
    auto sqrt(Quantity<UnitT> q) {
        return root<2>(q);
    }

    /**
     * @brief Cube root of a quantity
     *
     * @param q The quantity to take the cube root of
     * @return Quantity with dimension DimT^(1/3)
     */
    template<typename UnitT>
    auto cbrt(Quantity<UnitT> q) {
        return root<3>(q);
    }

    /**
     * @brief Absolute value of a quantity
     *
     * @param q The quantity
     * @return Quantity with the same unit, value = |q.value|
     */
    template<typename UnitT>
    auto abs(Quantity<UnitT> q) {
        return Quantity<UnitT>(std::abs(q.value));
    }

    /**
     * @brief Floating-point remainder of x / y
     *
     * @param x Dividend
     * @param y Divisor
     * @return x - n*y for some integer n, with the same unit as x and y
     */
    template<typename UnitT>
    Quantity<UnitT> fmod(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::fmod(x.value, y.value));
    }

    /**
     * @brief Decompose into integer and fractional parts
     *
     * @param x Input quantity
     * @param y Output: integer part (same unit as x)
     * @return Fractional part (same unit as x)
     *
     * Example:
     *   auto x = 3.14_m;
     *   auto frac = Utils::modf(x, int_part);  // int_part = 3_m, frac = 0.14_m
     */
    template<typename UnitT>
    Quantity<UnitT> modf(Quantity<UnitT> x, Quantity<UnitT> &y) {
        return Quantity<UnitT>(std::modf(x.value, &y.value));
    }

    /**
     * @brief Maximum of two quantities
     *
     * @param x First quantity
     * @param y Second quantity
     * @return The larger quantity
     */
    template<typename UnitT>
    Quantity<UnitT> max(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::max(x.value, y.value));
    }

    /**
     * @brief Minimum of two quantities
     *
     * @param x First quantity
     * @param y Second quantity
     * @return The smaller quantity
     */
    template<typename UnitT>
    Quantity<UnitT> min(Quantity<UnitT> x, Quantity<UnitT> y) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::min(x.value, y.value));
    }

    /**
     * @brief Clamp a quantity between min and max bounds
     *
     * @param x Value to clamp
     * @param min Lower bound (inclusive)
     * @param max Upper bound (inclusive)
     * @return If x < min, returns min; if x > max, returns max; otherwise returns x
     */
    template<typename UnitT>
    Quantity<UnitT> clamp(Quantity<UnitT> x, Quantity<UnitT> min, Quantity<UnitT> max) {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");
        return Quantity<UnitT>(std::clamp(x.value, min.value, max.value));
    }

    /**
     * @brief Floor of a quantity (largest integer <= x)
     *
     * @param x Input quantity
     * @return Quantity with the same unit, value = floor(x.value)
     */
    template<typename UnitT>
    Quantity<UnitT> floor(Quantity<UnitT> x) {
        return Quantity<UnitT>(std::floor(x.value));
    }

    /**
     * @brief Determine whether a quantity is nan
     *
     * @param x Input quantity
     * @return Return value = isnan(x.value)
     */
    template<typename UnitT>
    bool isnan(Quantity<UnitT> x) {
        return std::isnan(x.value);
    }

    /**
     * @brief Determine whether a quantity is finite
     *
     * @param x Input quantity
     * @return Return value = isfintie(x.value)
     */
    template<typename UnitT>
    bool isfinite(Quantity<UnitT> x) {
        return std::isfinite(x.value);
    }
}

#endif // PHYTH_QUANTITY_FUNCS_HPP
