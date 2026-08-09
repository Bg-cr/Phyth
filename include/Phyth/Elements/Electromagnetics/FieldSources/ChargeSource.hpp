#ifndef PHYTH_CHARGE_SOURCE_HPP
#define PHYTH_CHARGE_SOURCE_HPP

#include "Phyth/Tools/Vector3.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Abstract base class for all charge distributions
     *
     * ChargeSource defines the interface for any charge distribution that
     * produces an electric field and potential. Subclasses represent specific
     * geometries such as point charges, line charges, surface charges, volume
     * charges, or dipoles.
     *
     * The interface provides methods for computing:
     *   - Electric field E at a point
     *   - Electric potential V at a point
     *   - Total charge (integrated over the source)
     *   - Charge density at a point (for distributed sources)
     *   - Minimum distance from a point to the source (for regularization)
     *
     * Example:
     *   class PointCharge : public ChargeSource {
     *       // Implements GetElectricFieldAt, GetElectricPotentialAt,
     *       // GetTotalChargeValue, GetMinimumDistanceTo
     *   };
     */
    class ChargeSource {
    public:
        virtual ~ChargeSource() = default;

        /**
         * @brief Compute the electric field at a point
         *
         * @param point Position in 3D space (meters)
         * @return Electric field intensity E (N/C) as a Vector3
         *
         * For point charges, this uses Coulomb's law.
         * For distributed charges, this integrates over the source geometry.
         */
        [[nodiscard]] virtual Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const = 0;

        /**
         * @brief Compute the electric potential at a point
         *
         * @param point Position in 3D space (meters)
         * @return Electric potential V (Volts)
         *
         * For point charges, this uses V = k * q / r.
         * For distributed charges, this integrates over the source geometry.
         */
        [[nodiscard]] virtual Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const = 0;

        /**
         * @brief Get the total charge of this source
         *
         * @return Total charge (Coulombs)
         *
         * For point charges, this is the charge value itself.
         * For distributed charges, this is the integral of charge density.
         */
        [[nodiscard]] virtual Quantity<Coulomb> GetTotalChargeValue() const = 0;

        /**
         * @brief Get the charge density at a point
         *
         * @param point Position in 3D space (meters)
         * @return Charge density (C/m^3)
         *
         * The default implementation returns zero for non-distributed sources
         * (point charges, dipoles). Subclasses for distributed sources
         * (line, surface, volume charges) should override this.
         */
        [[nodiscard]] virtual Quantity<CoulombPerCubicMeter>
        GetChargeDensityAt(const Vector3<Quantity<Meter> > &point) const {
            return 0_Cpm3;
        }

        /**
         * @brief Get the minimum distance from a point to this charge source
         *
         * @param point Position in 3D space (meters)
         * @return Minimum distance (meters)
         *
         * This is used for regularization in field computations to avoid
         * singularities. For point charges, this returns the distance to the
         * point charge location. For extended sources, this returns the distance
         * to the nearest point on the source geometry.
         */
        [[nodiscard]] virtual Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const = 0;
    };

    /**
     * @brief Type trait to check if a type is a ChargeSource
     *
     * Example:
     *   static_assert(is_charge_source_v<PointCharge>);  // true
     *   static_assert(is_charge_source_v<int>);          // false
     */
    template<typename T>
    inline constexpr bool is_charge_source_v = std::is_base_of_v<ChargeSource, T>;
}

#endif // PHYTH_CHARGE_SOURCE_HPP
