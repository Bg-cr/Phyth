#ifndef PHYTH_DIPOLE_CHARGE_HPP
#define PHYTH_DIPOLE_CHARGE_HPP

#include "ChargeSource.hpp"

namespace Phyth::Electromagnetics {
    class DipoleCharge : public ChargeSource {
    public:
        DipoleCharge(const Vector3<Quantity<Meter>> &position,
                     Vector3<Quantity<CoulombMeter>> dipole_moment)
            : position_(position), p_(dipole_moment) {}

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at dipole location");
            }

            const auto r_hat = r.Normalized();
            const auto p_dot_r = p_.Dot(r_hat);
            const auto E = Consts::k_E * (3 * p_dot_r * r_hat - p_) / (dist * dist * dist);
            return E;
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at dipole location");
            }

            const auto r_hat = r.Normalized();
            return Consts::k_E * p_.Dot(r_hat) / (dist * dist);
        }

        [[nodiscard]] static Quantity<Coulomb> GetTotalChargeValue() override {
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
