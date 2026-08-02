#ifndef PHYTH_UNIT_HPP
#define PHYTH_UNIT_HPP

#include <ratio>

namespace Phyth {

    /**
     * @brief Physical unit with dimension and scale factor
     * @tparam DimT Dimension type (e.g., Length, Mass)
     * @tparam Scale std::ratio scale relative to base unit
     * @tparam TagT Tag for same-dimension differentiation
     */
    template<typename DimT, typename Scale = std::ratio<1>, typename TagT = void>
    struct Unit {
        using DimensionT = DimT;
        static constexpr double scale = static_cast<double>(Scale::num) / Scale::den;

        Unit() = default;
    };

    template <typename>
    struct is_unit : std::false_type {};

    template <typename DimT, typename Scale, typename TagT>
    struct is_unit<Unit<DimT, Scale, TagT>> : std::true_type {};

    template <typename T>
    inline constexpr bool is_unit_v = is_unit<T>::value;

}

#endif //PHYTH_UNIT_HPP