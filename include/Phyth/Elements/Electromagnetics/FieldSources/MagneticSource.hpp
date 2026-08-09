#ifndef PHYTH_MAGNETIC_SOURCE_HPP
#define PHYTH_MAGNETIC_SOURCE_HPP

#include "Phyth/Tools/Vector3.hpp"

namespace Phyth::Electromagnetics {

    /**
     * @brief Abstract base class for all magnetic field sources
     *
     * MagneticSource defines the interface for any object that produces
     * a magnetic field. Subclasses represent specific sources such as:
     *   - Current-carrying wire segments
     *   - Magnetic dipoles
     *   - Solenoids
     *   - Permanent magnets
     *
     * Example:
     *   class MagneticDipole : public MagneticSource {
     *       [[nodiscard]] Vector3<Quantity<Tesla>>
     *       GetMagneticFieldAt(const Vector3<Quantity<Meter>>& point) const override {
     *           // Compute dipole field at point
     *           return B;
     *       }
     *   };
     */
    class MagneticSource {
    public:
        virtual ~MagneticSource() = default;

        /**
         * @brief Compute the magnetic field at a point in space
         *
         * @param point Position in 3D space (meters)
         * @return Magnetic flux density B (Tesla) as a Vector3
         *
         * Subclasses must implement this method according to their
         * specific geometry and current distribution.
         *
         * @note: Some sources have singularities at certain locations
         *       (e.g., on the wire axis). Subclasses should document
         *       these edge cases.
         */
        [[nodiscard]] virtual Vector3<Quantity<Tesla>>
        GetMagneticFieldAt(const Vector3<Quantity<Meter>>& point) const = 0;
    };

    /**
     * @brief Type trait to check if a type is a MagneticSource
     *
     * Example:
     *   static_assert(is_magnetic_source_v<CurrentSegment>);  // true
     *   static_assert(is_magnetic_source_v<MagneticDipole>);   // true
     *   static_assert(is_magnetic_source_v<int>);              // false
     */
    template <typename T>
    inline constexpr bool is_magnetic_source_v = std::is_base_of_v<MagneticSource, T>;

}

#endif // PHYTH_MAGNETIC_SOURCE_HPP