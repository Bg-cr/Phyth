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
}

#endif //PHYTH_MAGNETIC_SOURCE_HPP
