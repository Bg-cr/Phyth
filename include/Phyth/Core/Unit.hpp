#ifndef PHYTH_UNIT_HPP
#define PHYTH_UNIT_HPP

#include <ratio>

namespace Phyth {
    /**
     * @brief Physical unit with dimension and scale factor relative to SI base
     *
     * Unit is the compile-time representation of a unit. It combines:
     *   - A dimension (e.g., Length, Mass, or derived like Velocity)
     *   - A scale factor relative to the SI base unit of that dimension
     *   - An optional tag for same-dimension differentiation
     *
     * @tparam DimT The dimension type (e.g., Length, Mass, Velocity)
     * @tparam Scale std::ratio scale relative to base unit (default: 1)
     * @tparam TagT An optional tag type to distinguish units with same dimension
     *              and scale (e.g., radians vs. degrees)
     *
     * Example:
     *   using Meter       = Unit<Length, std::ratio<1>>;
     *   using Kilometer   = Unit<Length, std::ratio<1000>>;
     *   using Second      = Unit<Time, std::ratio<1>>;
     *   using MilliSecond = Unit<Time, std::ratio<1, 1000>>;
     *
     *   // Same dimension, different scale, no tag
     *   static_assert(std::is_same_v<Meter::DimensionT, Length>);
     *   static_assert(Meter::scale == 1.0);
     *   static_assert(Kilometer::scale == 1000.0);
     *
     *   // Differentiating units with the same dimension and scale via TagT:
     *   struct RadianTag {};
     *   struct DegreeTag {};
     *   using Radian = Unit<Dimensionless, std::ratio<1>, RadianTag>;
     *   using Degree = Unit<Dimensionless, std::ratio<1>, DegreeTag>;
     *   // Radian and Degree are distinct types despite having the same
     *   // DimensionT and scale, because their TagT differ.
     *   static_assert(!std::is_same_v<Radian, Degree>);
     */
    template<typename DimT, typename Scale = std::ratio<1>, typename TagT = void>
    struct Unit {
        /** The dimension type of this unit */
        using DimensionT = DimT;
        using ScaleT = Scale;

        /**
         * @brief Scale factor relative to the SI base unit
         *
         * Value is computed at compile-time as Scale::num / Scale::den.
         * Examples:
         *   Meter::scale = 1.0
         *   Kilometer::scale = 1000.0
         *   MilliSecond::scale = 0.001
         */
        static constexpr double scale = static_cast<double>(Scale::num) / Scale::den;

        /** @brief Default constructor */
        Unit() = default;
    };

    /**
     * @brief Trait: check if a type is a Unit instantiation
     *
     * @tparam T The type to check
     *
     * Example:
     *   static_assert(is_unit_v<Unit<Length>>);    // true
     *   static_assert(is_unit_v<int>);             // false
     */
    template<typename T>
    struct is_unit : std::false_type {
    };

    template<typename DimT, typename Scale, typename TagT>
    struct is_unit<Unit<DimT, Scale, TagT> > : std::true_type {
    };

    /**
     * @brief Convenience variable template for is_unit
     *
     * Example:
     *   if constexpr (is_unit_v<T>) { ... }
     */
    template<typename T>
    inline constexpr bool is_unit_v = is_unit<T>::value;

    template<typename U1, typename U2>
    struct UnitMul {
        using type = Unit<DimMulT<typename U1::DimensionT, typename U1::DimensionT>, std::ratio_multiply<typename
            U1::ScaleT, typename U2::ScaleT> >;
    };

    template<typename U1, typename U2>
    using UnitMulT = typename UnitMul<U1, U2>::type;

    template<typename U1, typename U2>
    struct UnitDiv {
        using type = Unit<DimDivT<typename U1::DimensionT, typename U1::DimensionT>, std::ratio_divide<typename
            U1::ScaleT, typename U2::ScaleT> >;
    };

    template<typename U1, typename U2>
    using UnitDivT = typename UnitDiv<U1, U2>::type;

    template<typename U, typename Ratio = std::ratio<1> >
    struct UnitRec {
        using type = Unit<DimDivT<Dimensionless, typename U::DimensionT>, std::ratio_divide<Ratio, typename U::ScaleT> >
        ;
    };

    template<typename U1, typename Ratio = std::ratio<1> >
    using UnitRecT = typename UnitRec<U1, Ratio>::type;
}

#endif // PHYTH_UNIT_HPP
