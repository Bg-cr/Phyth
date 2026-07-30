#ifndef PHYTH_CHARGE_SOURCES_HPP
#define PHYTH_CHARGE_SOURCES_HPP

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    class PointCharge : public ChargeSource {
    public:
        PointCharge(const Quantity<Kilogram> mass, const Vector3<Quantity<Meter>> &position, const Quantity<Coulomb> charge)
            : ChargeSource(mass, position, charge) {
        }

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - GetPosition();
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at charge location");
            }
            return Consts::k_E * GetChargeValue() / (dist * dist) * r.Normalized();
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - GetPosition();
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at charge location");
            }
            return Consts::k_E * GetChargeValue() / dist;
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return GetChargeValue();
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const override {
            return (point - GetPosition()).Length();
        }
    };
}

#endif //PHYTH_CHARGE_SOURCES_HPP
