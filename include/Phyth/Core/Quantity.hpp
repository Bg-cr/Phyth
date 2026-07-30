#ifndef PHYTH_QUANTITY_HPP
#define PHYTH_QUANTITY_HPP

#include "Dimension.hpp"
#include "Unit.hpp"
#include "Phyth/System/UnitsRegistry.hpp"
#include "Phyth/System/Config.hpp"
#include <sstream>
#include <type_traits>

namespace Phyth {
    template<typename UnitT>
    struct Quantity {
        static_assert(is_unit_v<UnitT>, "UnitT is not a unit!");

        long double value;
        constexpr Quantity(const long double v = 0.0) : value(v) {}

        template<typename OtherUnit>
        constexpr Quantity(Quantity<OtherUnit> other)
            : value(other.value * OtherUnit::scale / UnitT::scale) {
            static_assert(std::is_same_v<typename UnitT::Dimension, typename OtherUnit::Dimension>,
                          "Cannot convert between different dimensions!");
        }

        template<typename TargetUnit>
        constexpr auto as() const {
            static_assert(std::is_same_v<typename UnitT::Dimension,
                                         typename TargetUnit::Dimension>,
                          "Cannot convert between different dimensions!");
            return Quantity<TargetUnit>(value * UnitT::scale / TargetUnit::scale);
        }

        template<typename OtherUnit>
        constexpr Quantity& operator+=(Quantity<OtherUnit> other) {
            static_assert(std::is_same_v<typename UnitT::Dimension, typename OtherUnit::Dimension>,
                          "Cannot add different dimensions!");
            value += other.value * OtherUnit::scale / UnitT::scale;
            return *this;
        }

        template<typename OtherUnit>
        constexpr Quantity& operator-=(Quantity<OtherUnit> other) {
            static_assert(std::is_same_v<typename UnitT::Dimension, typename OtherUnit::Dimension>,
                          "Cannot subtract different dimensions!");
            value -= other.value * OtherUnit::scale / UnitT::scale;
            return *this;
        }

        constexpr Quantity& operator*=(const long double scalar) {
            value *= scalar;
            return *this;
        }

        constexpr Quantity& operator/=(const long double scalar) {
            value /= scalar;
            return *this;
        }

        constexpr Quantity& operator++() {
            ++value;
            return *this;
        }

        constexpr Quantity& operator--() {
            --value;
            return *this;
        }

        constexpr Quantity operator++(int) {
            Quantity temp = *this;
            ++value;
            return temp;
        }

        constexpr Quantity operator--(int) {
            Quantity temp = *this;
            --value;
            return temp;
        }
    };

    template<typename>
    struct is_quantity : std::false_type {};

    template<typename UnitT>
    struct is_quantity<Quantity<UnitT>> : std::true_type {};

    template<typename T>
    inline constexpr bool is_quantity_v = is_quantity<T>::value;

    template<typename U1, typename U2>
    constexpr bool operator==(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::Dimension, typename U2::Dimension>,
                      "Cannot compare different dimensions!");
        return a.value == b.value * (U2::scale / U1::scale);
    }

    template<typename U1, typename U2>
    constexpr bool operator!=(Quantity<U1> a, Quantity<U2> b) {
        return !(a == b);
    }

    template<typename U1, typename U2>
    constexpr bool operator<(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::Dimension, typename U2::Dimension>,
                      "Cannot compare different dimensions!");
        return a.value < b.value * (U2::scale / U1::scale);
    }

    template<typename U1, typename U2>
    constexpr bool operator>(Quantity<U1> a, Quantity<U2> b) {
        return b < a;
    }

    template<typename U1, typename U2>
    constexpr bool operator<=(Quantity<U1> a, Quantity<U2> b) {
        return !(a > b);
    }

    template<typename U1, typename U2>
    constexpr bool operator>=(Quantity<U1> a, Quantity<U2> b) {
        return !(a < b);
    }

    template<typename UnitT>
    constexpr bool operator==(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value == scalar;
    }

    template<typename UnitT>
    constexpr bool operator==(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar == a.value;
    }

    template<typename UnitT>
    constexpr bool operator!=(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value != scalar;
    }

    template<typename UnitT>
    constexpr bool operator!=(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar != a.value;
    }

    template<typename UnitT>
    constexpr bool operator<(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value < scalar;
    }

    template<typename UnitT>
    constexpr bool operator<(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar < a.value;
    }

    template<typename UnitT>
    constexpr bool operator>(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value > scalar;
    }

    template<typename UnitT>
    constexpr bool operator>(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar > a.value;
    }

    template<typename UnitT>
    constexpr bool operator<=(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value <= scalar;
    }

    template<typename UnitT>
    constexpr bool operator<=(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar <= a.value;
    }

    template<typename UnitT>
    constexpr bool operator>=(Quantity<UnitT> a, long double scalar) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return a.value >= scalar;
    }

    template<typename UnitT>
    constexpr bool operator>=(long double scalar, Quantity<UnitT> a) {
        static_assert(is_dimensionless_v<UnitT>,
                      "Cannot compare dimensionful quantity with scalar!");
        return scalar >= a.value;
    }

    template<typename U1, typename U2>
    constexpr auto operator+(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::Dimension, typename U2::Dimension>,
                      "Cannot add different dimensions!");
        return Quantity<U1>(a.value + b.value * (U2::scale / U1::scale));
    }

    template<typename U1, typename U2>
    constexpr auto operator-(Quantity<U1> a, Quantity<U2> b) {
        static_assert(std::is_same_v<typename U1::Dimension, typename U2::Dimension>,
                      "Cannot subtract different dimensions!");
        return Quantity<U1>(a.value - b.value * (U2::scale / U1::scale));
    }

    template<typename U1, typename U2>
    constexpr auto operator*(Quantity<U1> a, Quantity<U2> b) {
        using Dim1 = typename U1::Dimension;
        using Dim2 = typename U2::Dimension;
        using ResultDim = DimMulT<Dim1, Dim2>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(a.value * b.value);
    }

    template<typename U1, typename U2>
    constexpr auto operator/(Quantity<U1> a, Quantity<U2> b) {
        using Dim1 = typename U1::Dimension;
        using Dim2 = typename U2::Dimension;
        using ResultDim = DimDivT<Dim1, Dim2>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(a.value / b.value);
    }

    template<typename UnitT>
    constexpr auto operator*(Quantity<UnitT> a, long double scalar) {
        return Quantity<UnitT>(a.value * scalar);
    }

    template<typename UnitT>
    constexpr auto operator*(long double scalar, Quantity<UnitT> a) {
        return Quantity<UnitT>(scalar * a.value);
    }

    template<typename UnitT>
    constexpr auto operator/(Quantity<UnitT> a, long double scalar) {
        return Quantity<UnitT>(a.value / scalar);
    }

    template<typename UnitT>
    constexpr auto operator/(long double scalar, Quantity<UnitT> a) {
        using DimT = typename UnitT::Dimension;
        using ResultDim = DimDivT<Dimensionless, DimT>;
        using ResultUnit = Unit<ResultDim>;
        return Quantity<ResultUnit>(scalar / a.value);
    }

    template<typename UnitT>
    constexpr auto operator-(Quantity<UnitT> a) {
        return Quantity<UnitT>(-a.value);
    }

    template<typename UnitT>
    constexpr auto operator+(Quantity<UnitT> a) {
        return Quantity<UnitT>(+a.value);
    }

    template<typename UnitT>
    std::string UnitSymbol() {
        return UnitRegistry::GetInstance().Symbol<UnitT>();
    }

    template<typename UnitT>
    std::string UnitName() {
        return UnitRegistry::GetInstance().Name<UnitT>();
    }

    template<typename>
    std::string FormatValue(const long double value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
    }

    template<typename UnitT>
    std::string QuantityToString(const Quantity<UnitT>& q) {
        if (Config::output_mode == OutputMode::Raw) {
            using DimT = typename UnitT::Dimension;
            return FormatValue<UnitT>(q.value) + " " + DimToString<DimT>::value();
        }
        const long double display_value = q.value / UnitT::scale;
        std::string unit = UnitSymbol<UnitT>();
        if (unit.empty()) {
            unit = UnitName<UnitT>();
        }
        if (unit.empty()) {
            return FormatValue<UnitT>(display_value);
        }
        return FormatValue<UnitT>(display_value) + " " + unit;
    }

    template<typename UnitT>
    std::ostream& operator<<(std::ostream& os, const Quantity<UnitT>& q) {
        os << QuantityToString(q);
        return os;
    }
}

template<typename U1, typename U2, typename = void>
struct common_type_impl;

template<typename U1, typename U2>
struct common_type_impl<U1, U2, std::enable_if_t<Phyth::is_unit_v<U1> && Phyth::is_unit_v<U2>>> {
    using type = Phyth::Quantity<U1>;
};

template<typename U1, typename U2>
struct common_type_impl<U1, U2, void> {};

template<typename U1, typename U2>
struct std::common_type<Phyth::Quantity<U1>, Phyth::Quantity<U2>>
    : common_type_impl<U1, U2> {};

#endif //PHYTH_QUANTITY_HPP