#ifndef PHYTH_UNIT_HPP
#define PHYTH_UNIT_HPP

#include <ratio>

namespace Phyth {
    template<typename DimT, typename Scale = std::ratio<1>, typename = void>
    struct Unit {
        using Dimension = DimT;
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
