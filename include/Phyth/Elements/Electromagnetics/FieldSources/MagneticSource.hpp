#ifndef PHYTH_MAGNETIC_SOURCE_HPP
#define PHYTH_MAGNETIC_SOURCE_HPP
#include "Phyth/Tools/Vector3.hpp"


namespace Phyth::Electromagnetics {
    class MagneticSource {
    public:
        virtual ~MagneticSource() = default;

        [[nodiscard]] virtual Vector3<Quantity<Tesla>>
        GetMagneticFieldAt(const Vector3<Quantity<Meter>>& point) const = 0;
    };

    template <typename T>
    inline constexpr bool is_magnetic_source_v = std::is_base_of_v<MagneticSource, T>;
}

#endif //PHYTH_MAGNETIC_SOURCE_HPP
