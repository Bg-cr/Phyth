#ifndef PHYTH_CHARGE_PARTICLE_HPP
#define PHYTH_CHARGE_PARTICLE_HPP

#include "Phyth/Elements/Mechanics/Particle.hpp"

namespace Phyth::Electromagnetics {
    class ChargeParticle : public Mechanics::Particle {
    public:
        ChargeParticle(const Quantity<Kilogram> mass, const Vector3<Quantity<Meter>> &position, const Quantity<Coulomb> charge_value)
            : Particle(mass, position), charge_value_(charge_value) {}

        [[nodiscard]] Quantity<Coulomb> GetChargeValue() const { return charge_value_; }
        void SetChargeValue(const Quantity<Coulomb> charge_value) { charge_value_ = charge_value; }

    private:
        Quantity<Coulomb> charge_value_;
    };
}

#endif //PHYTH_CHARGE_PARTICLE_HPP
