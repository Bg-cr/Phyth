#ifndef PHYTH_FORCE_SOURCE_HPP
#define PHYTH_FORCE_SOURCE_HPP

#include "Phyth/Tools/Vector3.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Matrix3x3.hpp"

namespace Phyth::Mechanics {
    /**
     * @brief Abstract base class for all force/traction distributions
     *
     * ForceSource defines the interface for any mechanical loading that
     * produces a stress field. Subclasses represent specific loading
     * geometries such as point forces, surface tractions, body forces,
     * or distributed loads.
     *
     * The interface provides methods for computing:
     *   - Stress tensor σ at a point (primary for coupling)
     *   - Force vector at a point (simplified view)
     *   - Total force (integrated over the source)
     *   - Body force density at a point (for distributed sources)
     *   - Minimum distance from a point to the source (for regularization)
     *
     * Example:
     *   class PointForce : public ForceSource {
     *       // Implements GetStressAt, GetForceAt, GetTotalForce,
     *       // GetBodyForceDensityAt, GetMinimumDistanceTo
     *   };
     */
    class ForceSource {
    public:
        virtual ~ForceSource() = default;

        /**
         * @brief Compute the stress tensor at a point
         *
         * @param point Position in 3D space (meters)
         * @return Stress tensor σ (Pa) as a 3x3 matrix
         *
         * This is the primary method for mechanics-EM coupling.
         * For point forces, this uses Kelvin's solution.
         * For distributed loads, this integrates over the source geometry.
         */
        [[nodiscard]] virtual Matrix3x3<Quantity<Pascal>> GetStressAt(
            const Vector3<Quantity<Meter>>& point
        ) const = 0;

        /**
         * @brief Compute the force vector at a point
         *
         * @param point Position in 3D space (meters)
         * @return Force vector (Newtons)
         *
         * This provides a simplified view of the loading.
         * For point forces, this returns the force at the point.
         * For distributed loads, this may be zero (use GetBodyForceDensityAt instead).
         */
        [[nodiscard]] virtual Vector3<Quantity<Newton>> GetForceAt(
            const Vector3<Quantity<Meter>>& point
        ) const = 0;

        /**
         * @brief Get the total force of this source
         *
         * @return Total force (Newtons) as a vector
         *
         * For point forces, this is the force value itself.
         * For distributed loads, this is the integral of body force density.
         */
        [[nodiscard]] virtual Vector3<Quantity<Newton>> GetTotalForce() const = 0;

        /**
         * @brief Get the body force density at a point
         *
         * @param point Position in 3D space (meters)
         * @return Body force density (N/m^3)
         *
         * The default implementation returns zero for non-distributed sources
         * (point forces, surface tractions). Subclasses for distributed sources
         * (volume forces) should override this.
         */
        [[nodiscard]] virtual Vector3<Quantity<NewtonPerMeterCubed>>
        GetBodyForceDensityAt(const Vector3<Quantity<Meter>>& point) const {
            return Vector3<Quantity<NewtonPerMeterCubed>>();
        }

        /**
         * @brief Get the minimum distance from a point to this force source
         *
         * @param point Position in 3D space (meters)
         * @return Minimum distance (meters)
         *
         * This is used for regularization in field computations to avoid
         * singularities. For point forces, this returns the distance to the
         * point force location. For extended sources, this returns the distance
         * to the nearest point on the source geometry.
         */
        [[nodiscard]] virtual Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const = 0;
    };

    /**
     * @brief Type trait to check if a type is a ForceSource
     *
     * Example:
     *   static_assert(is_force_source_v<PointForce>);  // true
     *   static_assert(is_force_source_v<int>);         // false
     */
    template<typename T>
    inline constexpr bool is_force_source_v = std::is_base_of_v<ForceSource, T>;
}

#endif // PHYTH_FORCE_SOURCE_HPP