#ifndef PHYTH_DIPOLE_CHARGE_HPP
#define PHYTH_DIPOLE_CHARGE_HPP

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Ideal electric dipole with zero net charge
     *
     * DipoleCharge represents an ideal point electric dipole. It has zero
     * net charge (total charge = 0) but produces a non-zero electric field
     * and potential due to the separation of positive and negative charges.
     *
     * The dipole is defined by its position and dipole moment vector p,
     * which points from the negative to the positive charge.
     *
     * Example:
     *   auto pos = Vector3<Quantity<Meter>>(0_m, 0_m, 0_m);
     *   auto p = Vector3<Quantity<CoulombMeter>>(1.0_C*m, 0_C*m, 0_C*m);
     *   DipoleCharge dipole(pos, p);
     *
     *   auto point = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     *   auto E = dipole.GetElectricFieldAt(point);
     *   auto V = dipole.GetElectricPotentialAt(point);
     *
     * @note: The field diverges at the dipole location. Querying the field
     *       at the dipole position throws an exception.
     */
    class DipoleCharge : public ChargeSource {
    public:
        /**
         * @brief Construct an electric dipole
         *
         * @param position Location of the dipole (meters)
         * @param dipole_moment Dipole moment vector p (C*m)
         */
        DipoleCharge(const Vector3<Quantity<Meter> > &position,
                     const Vector3<Quantity<CoulombMeter> > &dipole_moment)
            : position_(position), p_(dipole_moment) {
        }

        /**
         * @brief Compute the electric field of the dipole
         *
         * The field at point r (relative to the dipole) is:
         *   E = k_E * (3 * (p dot hat_r) * hat_r - p) / r^3
         *
         * @param point Position to evaluate the field (meters)
         * @return Electric field intensity E (N/C)
         *
         * @throws std::runtime_error If point is at the dipole location
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("DipoleCharge: field diverges at dipole location");
            }

            const auto hat_r = r.Normalized();
            return Consts::k_E * (3 * p_.Dot(hat_r) * hat_r - p_) / (dist * dist * dist);
        }

        /**
         * @brief Compute the electric potential of the dipole
         *
         * The potential at point r (relative to the dipole) is:
         *   V = k_E * (p dot hat_r) / r^2
         *
         * @param point Position to evaluate the potential (meters)
         * @return Electric potential V (Volts)
         *
         * @throws std::runtime_error If point is at the dipole location
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("DipoleCharge: potential diverges at dipole location");
            }

            const auto hat_r = r.Normalized();
            return Consts::k_E * p_.Dot(hat_r) / (dist * dist);
        }

        /**
         * @brief Get the total charge of the dipole
         *
         * @return 0_C (an ideal dipole has zero net charge)
         */
        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return 0_C;
        }

        /**
         * @brief Get the minimum distance from a point to the dipole location
         *
         * @param point Position to compute distance from
         * @return Distance from the point to the dipole location (meters)
         */
        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            return (point - position_).Length();
        }

        /** @return The dipole moment vector p (C*m) */
        [[nodiscard]] Vector3<Quantity<CoulombMeter> > GetDipoleMoment() const {
            return p_;
        }

    private:
        Vector3<Quantity<Meter> > position_;
        Vector3<Quantity<CoulombMeter> > p_;
    };
} // namespace Phyth::Electromagnetics

#endif // PHYTH_DIPOLE_CHARGE_HPP
