#ifndef PHYTH_POINT_CHARGE_HPP
#define PHYTH_POINT_CHARGE_HPP

#include "ChargeSource.hpp"
#include "MagneticSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Point charge with Coulomb field and potential
     *
     * Field and potential diverge at the charge location.
     */
    class PointCharge : public Mechanics::Particle,
                        public ChargeSource,
                        public MagneticSource {
    public:
        explicit PointCharge(const Quantity<Kilogram> mass, const Vector3<Quantity<Meter>> &position, const Quantity<Coulomb> charge_value)
            : Particle(mass, position), charge_value_(charge_value) {
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return k_E * Q / r^2 * hat_r
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at charge location");
            }
            return Consts::k_E * charge_value_ / (dist * dist) * r.Normalized();
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return mu_0 / 4pi * q * (v times r) / r^3
         */
        [[nodiscard]] Vector3<Quantity<Tesla>>
        GetMagneticFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();
            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at charge location");
            }
            return Consts::mu_0 / Consts::varpi * charge_value_ * velocity_.Cross(r) / Utils::cube(dist);
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return k_E * Q / r
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
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
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            return (point - position_).Length();
        }

    private:
        Quantity<Coulomb> charge_value_;
    };
}

#endif //PHYTH_POINT_CHARGE_HPP
