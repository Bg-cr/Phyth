#ifndef PHYTH_POINT_CHARGE_HPP
#define PHYTH_POINT_CHARGE_HPP

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {

    /**
     * @brief Point charge with Coulomb field and potential
     *
     * Field and potential diverge at the charge location.
     */
    class PointCharge : public ChargeSource {
    public:
        explicit PointCharge(const Vector3<Quantity<Meter>>& position, const Quantity<Coulomb> charge_value)
            : position_(position), charge_value_(charge_value) {}

        /**
         * @throw: std::runtime_error when evaluated at distance < Config::epsilon
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at charge location");
            }
            return Consts::k_E * charge_value_ / (dist * dist) * r.Normalized();
        }

        /**
         * @throw: std::runtime_error when evaluated at distance < Config::epsilon
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at charge location");
            }
            return Consts::k_E * charge_value_ / dist;
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return charge_value_;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const override {
            return (point - position_).Length();
        }

        [[nodiscard]] Vector3<Quantity<Meter>> GetPosition() const { return position_; }
        void SetPosition(const Vector3<Quantity<Meter>>& position) { position_ = position; }

    private:
        Vector3<Quantity<Meter>> position_;
        Quantity<Coulomb> charge_value_;
    };

}

#endif //PHYTH_POINT_CHARGE_HPP