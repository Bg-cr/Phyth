#ifndef PHYTH_CHARGE_PARTICLE_HPP
#define PHYTH_CHARGE_PARTICLE_HPP

#include <memory>

#include "Phyth/Elements/Mechanics/Particle.hpp"
#include "ChargeSources/PointCharge.hpp"

namespace Phyth::Electromagnetics {

    /**
     * @brief Mechanical particle carrying a point charge
     *
     * Combines a Mechanics::Particle with a PointCharge source.
     * Position is synchronized between the particle and its charge.
     */
    class ChargeParticle : public Mechanics::Particle {
    public:
        ChargeParticle(const Quantity<Kilogram> mass, std::shared_ptr<PointCharge> charge)
            : Particle(mass, charge->GetPosition()), charge_(std::move(charge)) {}

        /**
         * @brief Sets position and synchronizes with the charge
         */
        void SetPosition(const Vector3<Quantity<Meter>>& position) {
            position_ = position;
            charge_->SetPosition(position_);
        }

        void SetCharge(std::shared_ptr<PointCharge> charge) { charge_ = std::move(charge); }
        [[nodiscard]] std::shared_ptr<const PointCharge> GetCharge() const { return charge_; }

    private:
        std::shared_ptr<PointCharge> charge_;
    };

}

#endif //PHYTH_CHARGE_PARTICLE_HPP