#ifndef PHYTH_DIPOLE_CHARGE_HPP
#define PHYTH_DIPOLE_CHARGE_HPP

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {

    /** @brief Ideal electric dipole with zero net charge */
    class DipoleCharge : public ChargeSource {
    public:
        DipoleCharge(const Vector3<Quantity<Meter>>& position,
                     const Vector3<Quantity<CoulombMeter>>& dipole_moment)
            : position_(position), p_(dipole_moment) {}

        /** @return k_E * (3 * (p dot hat_r) * hat_r - p) / r^3 */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at dipole location");
            }

            const auto hat_r = r.Normalized();
            return Consts::k_E * (3 * p_.Dot(hat_r) * hat_r - p_) / (dist * dist * dist);
        }

        /** @return k_E * (p dot hat_r) / r^2 */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at dipole location");
            }

            const auto hat_r = r.Normalized();
            return Consts::k_E * p_.Dot(hat_r) / (dist * dist);
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return 0_C;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const override {
            return (point - position_).Length();
        }

        [[nodiscard]] Vector3<Quantity<CoulombMeter>> GetDipoleMoment() const {
            return p_;
        }

    private:
        Vector3<Quantity<Meter>> position_;
        Vector3<Quantity<CoulombMeter>> p_;
    };

}

#endif //PHYTH_DIPOLE_CHARGE_HPP