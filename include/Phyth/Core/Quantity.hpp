#ifndef PHYTH_QUANTITY_HPP
#define PHYTH_QUANTITY_HPP

#include "Dimension.hpp"
#include "Unit.hpp"
#include "Phyth/System/UnitsRegistry.hpp"
#include "Phyth/System/Config.hpp"
#include <sstream>
#include <type_traits>

namespace Phyth {
    /**
     * @brief Physical quantity with a value and a compile-time unit
     *
     * Quantity is the primary user-facing type. It stores a numeric value and
     * carries a unit type that encodes both dimension (SI exponents) and scale
     * (e.g., meter vs. kilometer).
     *
     * @tparam UnitT A unit type (e.g., Unit<Length, std::ratio<1>>)
     *
     * Example:
     *   Quantity<Meter> length = 5.0_m;
     *   Quantity<Second> time = 3.0_s;
     *   auto speed = length / time; // Quantity<Velocity> with value ~1.666...
     */
    template<typename UnitT>
    struct Quantity {
        static_assert(is_unit_v<UnitT>, "UnitT must be a unit type");

        /** The numeric value of this quantity, stored in the unit's base scale */
        double value;

        /** @brief Construct a quantity with a given value (default: 0) */
        constexpr Quantity(const double v = 0.0) : value(v) {
        }

        /**
         * @brief Convert from another unit of the same dimension
         *
         * The stored value is converted to this unit's scale.
         *
         * @tparam OtherUnit The source unit type
         * @param other The quantity to convert from
         */
        template<typename OtherUnit>
        constexpr Quantity(Quantity<OtherUnit> other)
            : value(other.value * OtherUnit::scale / UnitT::scale) {
            static_assert(std::is_same_v<typename UnitT::DimensionT, typename OtherUnit::DimensionT>,
                          "Cannot convert between different dimensions");
        }

        /**
         * @brief Convert this quantity to a different unit of the same dimension
         *
         * @tparam TargetUnit The target unit type
         * @return Quantity<TargetUnit> with the same physical value
         *
         * Example:
         *   auto length = 1.0_km;
         *   auto in_meters = length.as<Meter>(); // 1000 m
         */
        template<typename TargetUnit>
        constexpr auto as() const {
            static_assert(std::is_same_v<typename UnitT::DimensionT,
                              typename TargetUnit::DimensionT>,
                          "Cannot convert between different dimensions");
            return Quantity<TargetUnit>(value * UnitT::scale / TargetUnit::scale);
        }

        /**
         * @brief Extract the numeric value as a different type
         *
         * @tparam TargetType The numeric type to convert to (e.g., float, int)
         * @return The value cast to TargetType
         */
        template<typename TargetType>
        constexpr auto to() const {
            return static_cast<TargetType>(value);
        }

        /**
         * @brief Add another quantity of the same dimension
         *
         * @tparam OtherUnit The other quantity's unit type
         * @param other The quantity to add
         * @return Reference to this quantity
         */
        template<typename OtherUnit>
        constexpr Quantity &operator+=(Quantity<OtherUnit> other) {
            static_assert(std::is_same_v<typename UnitT::DimensionT, typename OtherUnit::DimensionT>,
                          "Cannot add different dimensions");
            value += other.value * OtherUnit::scale / UnitT::scale;
            return *this;
        }

        /**
         * @brief Subtract another quantity of the same dimension
         *
         * @tparam OtherUnit The other quantity's unit type
         * @param other The quantity to subtract
         * @return Reference to this quantity
         */
        template<typename OtherUnit>
        constexpr Quantity &operator-=(Quantity<OtherUnit> other) {
            static_assert(std::is_same_v<typename UnitT::DimensionT, typename OtherUnit::DimensionT>,
                          "Cannot subtract different dimensions");
            value -= other.value * OtherUnit::scale / UnitT::scale;
            return *this;
        }

        /** @brief Multiply the value by a scalar */
        constexpr Quantity &operator*=(const double scalar) {
            value *= scalar;
            return *this;
        }

        /** @brief Divide the value by a scalar */
        constexpr Quantity &operator/=(const double scalar) {
            value /= scalar;
            return *this;
        }

        /** @brief Pre-increment: increase value by 1 */
        constexpr Quantity &operator++() {
            ++value;
            return *this;
        }

        /** @brief Pre-decrement: decrease value by 1 */
        constexpr Quantity &operator--() {
            --value;
            return *this;
        }

        /** @brief Post-increment: increase value by 1, return old value */
        constexpr Quantity operator++(int) {
            Quantity temp = *this;
            ++value;
            return temp;
        }

        /** @brief Post-decrement: decrease value by 1, return old value */
        constexpr Quantity operator--(int) {
            Quantity temp = *this;
            --value;
            return temp;
        }
    };

    /**
     * @brief Trait: check if a type is a Quantity instantiation
     */
    template<typename>
    struct is_quantity : std::false_type {
    };

    template<typename UnitT>
    struct is_quantity<Quantity<UnitT> > : std::true_type {
    };

    /**
     * @brief Convenience variable template for is_quantity
     */
    template<typename T>
    inline constexpr bool is_quantity_v = is_quantity<T>::value;


    /**
     * @brief Equality comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if both quantities represent the same physical value
     *
     * Example:
     *   1.0_km == 1000.0_m  // true
     *   1.0_km == 999.0_m   // false
     */
    template<typename U1, typename U2>
    constexpr bool operator==(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::DimensionT, typename U2::DimensionT>,
                      "Cannot compare different dimensions");
        return a.value == b.value * (U2::scale / U1::scale);
    }

    /**
     * @brief Inequality comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if the quantities represent different physical values
     *
     * Example:
     *   1.0_km != 1000.0_m  // false
     *   1.0_km != 999.0_m   // true
     */
    template<typename U1, typename U2>
    constexpr bool operator!=(Quantity<U1> a, Quantity<U2> b) {
        return !(a == b);
    }

    /**
     * @brief Less-than comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if a is physically smaller than b
     *
     * Example:
     *   1.0_m < 100.0_cm   // false (equal)
     *   1.0_m < 150.0_cm   // true  (1m < 1.5m)
     */
    template<typename U1, typename U2>
    constexpr bool operator<(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::DimensionT, typename U2::DimensionT>,
                      "Cannot compare different dimensions");
        return a.value < b.value * (U2::scale / U1::scale);
    }

    /**
     * @brief Greater-than comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if a is physically larger than b
     *
     * Example:
     *   1.0_m > 100.0_cm   // false (equal)
     *   1.0_m > 80.0_cm    // true  (1m > 0.8m)
     */
    template<typename U1, typename U2>
    constexpr bool operator>(Quantity<U1> a, Quantity<U2> b) {
        return b < a;
    }

    /**
     * @brief Less-than-or-equal comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if a is physically smaller than or equal to b
     *
     * Example:
     *   1.0_m <= 100.0_cm   // true (equal)
     *   1.0_m <= 150.0_cm   // true
     *   1.0_m <= 80.0_cm    // false
     */
    template<typename U1, typename U2>
    constexpr bool operator<=(Quantity<U1> a, Quantity<U2> b) {
        return !(a > b);
    }

    /**
     * @brief Greater-than-or-equal comparison between two quantities
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return true if a is physically larger than or equal to b
     *
     * Example:
     *   1.0_m >= 100.0_cm   // true (equal)
     *   1.0_m >= 80.0_cm    // true
     *   1.0_m >= 150.0_cm   // false
     */
    template<typename U1, typename U2>
    constexpr bool operator>=(Quantity<U1> a, Quantity<U2> b) {
        return !(a < b);
    }

    /**
     * @brief Equality comparison: dimensionless quantity == scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value == scalar
     *
     * Example:
     *   auto angle = 1.0_rad;
     *   if (angle == 1.0) { ... }  // true
     */
    template<typename UnitT>
    constexpr bool operator==(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value == scalar;
    }

    /**
     * @brief Equality comparison: scalar == dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar == a.value
     *
     * Example:
     *   auto angle = 1.0_rad;
     *   if (1.0 == angle) { ... }  // true
     */
    template<typename UnitT>
    constexpr bool operator==(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar == a.value;
    }

    /**
     * @brief Inequality comparison: dimensionless quantity != scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value != scalar
     */
    template<typename UnitT>
    constexpr bool operator!=(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value != scalar;
    }

    /**
     * @brief Inequality comparison: scalar != dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar != a.value
     */
    template<typename UnitT>
    constexpr bool operator!=(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar != a.value;
    }

    /**
     * @brief Less-than comparison: dimensionless quantity < scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value < scalar
     */
    template<typename UnitT>
    constexpr bool operator<(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value < scalar;
    }

    /**
     * @brief Less-than comparison: scalar < dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar < a.value
     */
    template<typename UnitT>
    constexpr bool operator<(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar < a.value;
    }

    /**
     * @brief Greater-than comparison: dimensionless quantity > scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value > scalar
     */
    template<typename UnitT>
    constexpr bool operator>(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value > scalar;
    }

    /**
     * @brief Greater-than comparison: scalar > dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar > a.value
     */
    template<typename UnitT>
    constexpr bool operator>(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar > a.value;
    }

    /**
     * @brief Less-than-or-equal comparison: dimensionless quantity <= scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value <= scalar
     */
    template<typename UnitT>
    constexpr bool operator<=(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value <= scalar;
    }

    /**
     * @brief Less-than-or-equal comparison: scalar <= dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar <= a.value
     */
    template<typename UnitT>
    constexpr bool operator<=(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar <= a.value;
    }

    /**
     * @brief Greater-than-or-equal comparison: dimensionless quantity >= scalar
     *
     * @tparam UnitT Must be dimensionless
     * @param a The dimensionless quantity
     * @param scalar The scalar value to compare against
     * @return true if a.value >= scalar
     */
    template<typename UnitT>
    constexpr bool operator>=(Quantity<UnitT> a, double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return a.value >= scalar;
    }

    /**
     * @brief Greater-than-or-equal comparison: scalar >= dimensionless quantity
     *
     * @tparam UnitT Must be dimensionless
     * @param scalar The scalar value
     * @param a The dimensionless quantity
     * @return true if scalar >= a.value
     */
    template<typename UnitT>
    constexpr bool operator>=(double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare a dimensionful quantity with a scalar");
        return scalar >= a.value;
    }

    /**
     * @brief Add two quantities of the same dimension
     *
     * @tparam U1 Unit type of the left operand (also the result unit)
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return Quantity<U1> with value a + b converted to U1's scale
     *
     * Example:
     *   1.0_m + 50.0_cm  // Quantity<Meter>(1.5)
     */
    template<typename U1, typename U2>
    constexpr auto operator+(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::DimensionT, typename U2::DimensionT>,
                      "Cannot add different dimensions");
        return Quantity<U1>(a.value + b.value * (U2::scale / U1::scale));
    }

    /**
     * @brief Subtract two quantities of the same dimension
     *
     * @tparam U1 Unit type of the left operand (also the result unit)
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return Quantity<U1> with value a - b converted to U1's scale
     *
     * Example:
     *   1.0_m - 50.0_cm  // Quantity<Meter>(0.5)
     */
    template<typename U1, typename U2>
    constexpr auto operator-(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::DimensionT, typename U2::DimensionT>,
                      "Cannot subtract different dimensions");
        return Quantity<U1>(a.value - b.value * (U2::scale / U1::scale));
    }

    /**
     * @brief Multiply two quantities
     *
     * The resulting dimension is the product of the two input dimensions.
     * The resulting unit is dimensionless (scale = 1) with the combined dimension.
     *
     * @tparam U1 Unit type of the left operand
     * @tparam U2 Unit type of the right operand
     * @param a Left quantity
     * @param b Right quantity
     * @return Quantity<Unit<DimMulT<Dim1, Dim2>>>
     *
     * Example:
     *   5.0_m * 3.0_s  // Quantity<Unit<DimMulT<Length, Time>>>(15)
     */
    template<typename U1, typename U2>
    constexpr auto operator*(Quantity<U1> a, Quantity<U2> b) {
        using Dim1 = typename U1::DimensionT;
        using Dim2 = typename U2::DimensionT;
        using ResultDim = DimMulT<Dim1, Dim2>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(a.value * b.value);
    }

    /**
     * @brief Divide one quantity by another
     *
     * The resulting dimension is the quotient of the two input dimensions.
     * The resulting unit is dimensionless (scale = 1) with the combined dimension.
     *
     * @tparam U1 Unit type of the numerator
     * @tparam U2 Unit type of the denominator
     * @param a Numerator quantity
     * @param b Denominator quantity
     * @return Quantity<Unit<DimDivT<Dim1, Dim2>>>
     *
     * Example:
     *   10.0_m / 2.0_s  // Quantity<Unit<Velocity>>(5)
     */
    template<typename U1, typename U2>
    constexpr auto operator/(Quantity<U1> a, Quantity<U2> b) {
        using Dim1 = typename U1::DimensionT;
        using Dim2 = typename U2::DimensionT;
        using ResultDim = DimDivT<Dim1, Dim2>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(a.value / b.value);
    }

    /**
     * @brief Multiply a quantity by a scalar
     *
     * @param a The quantity
     * @param scalar The scalar multiplier
     * @return Quantity<UnitT> with value a.value * scalar
     */
    template<typename UnitT>
    constexpr auto operator*(Quantity<UnitT> a, double scalar) {
        return Quantity<UnitT>(a.value * scalar);
    }

    /**
     * @brief Multiply a scalar by a quantity (commutative)
     *
     * @param scalar The scalar multiplier
     * @param a The quantity
     * @return Quantity<UnitT> with value scalar * a.value
     */
    template<typename UnitT>
    constexpr auto operator*(double scalar, Quantity<UnitT> a) {
        return Quantity<UnitT>(scalar * a.value);
    }

    /**
     * @brief Divide a quantity by a scalar
     *
     * @param a The quantity
     * @param scalar The scalar divisor
     * @return Quantity<UnitT> with value a.value / scalar
     */
    template<typename UnitT>
    constexpr auto operator/(Quantity<UnitT> a, double scalar) {
        return Quantity<UnitT>(a.value / scalar);
    }

    /**
     * @brief Divide a scalar by a quantity, producing an inverse unit
     *
     * This is the only way to create inverse dimensions from a scalar.
     *
     * @tparam UnitT The unit type of the denominator
     * @param scalar The scalar numerator
     * @param a The quantity denominator
     * @return Quantity<Unit<DimDivT<Dimensionless, DimT>>>
     *
     * Example:
     *   auto freq = 1.0 / 2.0_s;  // Quantity<Unit<DimRcpT<Time>>>(0.5)
     *   // The dimension is s^-1, which is frequency.
     */
    template<typename UnitT>
    constexpr auto operator/(double scalar, Quantity<UnitT> a) {
        using DimT = typename UnitT::DimensionT;
        using ResultDim = DimDivT<Dimensionless, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(scalar / a.value);
    }

    /**
     * @brief Unary negation
     *
     * @param a The quantity to negate
     * @return Quantity<UnitT> with value -a.value
     */
    template<typename UnitT>
    constexpr auto operator-(Quantity<UnitT> a) {
        return Quantity<UnitT>(-a.value);
    }

    /**
     * @brief Unary plus
     *
     * @param a The quantity
     * @return Quantity<UnitT> with value +a.value (a copy)
     */
    template<typename UnitT>
    constexpr auto operator+(Quantity<UnitT> a) {
        return Quantity<UnitT>(+a.value);
    }

    /**
     * @brief Get the symbol string for a unit type (e.g., "m", "kg", "N")
     *
     * @tparam UnitT The unit type to look up
     * @return The unit symbol, or empty string if not registered
     */
    template<typename UnitT>
    std::string UnitSymbol() {
        return UnitRegistry::GetInstance().GetSymbol<UnitT>();
    }

    /**
     * @brief Get the full name string for a unit type (e.g., "meter", "kilogram")
     *
     * @tparam UnitT The unit type to look up
     * @return The unit name, or empty string if not registered
     */
    template<typename UnitT>
    std::string UnitName() {
        return UnitRegistry::GetInstance().GetName<UnitT>();
    }

    /**
     * @brief Format a raw numeric value for string output
     *
     * Uses std::ostringstream. Override this function for custom formatting.
     *
     * @param value The value to format
     * @return String representation of the value
     */
    inline std::string FormatValue(const double value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
    }

    /**
     * @brief Convert a Quantity to a human-readable string
     *
     * Output format depends on Config::output_mode:
     *   - Raw: "<value> <dimension-string>" (e.g., "5 m*s^-1")
     *   - Auto: "<value> <transformed dimension-string>" (e.g., "7 J")
     *
     * If the unit has no registered symbol or name, falls back to dimension string.
     *
     * @tparam UnitT The unit type of the quantity
     * @param q The quantity to convert
     * @return String representation of the quantity
     */
    template<typename UnitT>
    std::string QuantityToString(const Quantity<UnitT> &q) {
        using DimT = typename UnitT::DimensionT;
        if (Config::output_mode == OutputMode::Raw) {
            return FormatValue(q.value) + " " + DimToString<DimT>::value();
        }
        const double display_value = q.value / UnitT::scale;
        std::string unit = UnitSymbol<UnitT>();
        if (unit.empty()) {
            unit = UnitName<UnitT>();
        }
        if (unit.empty()) {
            return FormatValue(display_value) + " " + DimToString<DimT>::value();
        }
        return FormatValue(display_value) + " " + unit;
    }

    /**
     * @brief Stream insertion operator for Quantity
     *
     * @tparam UnitT The unit type of the quantity
     * @param os The output stream
     * @param q The quantity to output
     * @return Reference to the output stream
     */
    template<typename UnitT>
    std::ostream &operator<<(std::ostream &os, const Quantity<UnitT> &q) {
        os << QuantityToString(q);
        return os;
    }
}

/**
 * @brief std::common_type specialization for Quantity
 *
 * Allows mixing Quantity types in templates that use common_type.
 * Prefers the left operand's unit as the common type.
 *
 * @tparam U1 Unit type of the first Quantity
 * @tparam U2 Unit type of the second Quantity
 *
 * Example:
 *   std::common_type_t<Quantity<Meter>, Quantity<Centimeter>>
 *   // -> Quantity<Meter>
 */
template<typename U1, typename U2, typename = void>
struct common_type_impl;

template<typename U1, typename U2>
struct common_type_impl<U1, U2, std::enable_if_t<Phyth::is_unit_v<U1> && Phyth::is_unit_v<U2>> > {
    using type = Phyth::Quantity<U1>;
};

template<typename U1, typename U2>
struct common_type_impl<U1, U2, void> {
};

template<typename U1, typename U2>
struct std::common_type<Phyth::Quantity<U1>, Phyth::Quantity<U2> >
        : common_type_impl<U1, U2> {
};

#endif // PHYTH_QUANTITY_HPP
