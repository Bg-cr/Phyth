#ifndef PHYTH_PARTICLE_H
#define PHYTH_PARTICLE_H

#include <functional>

#include "Phyth/Physical/PhysicalConsts.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Vector3.hpp"

namespace Phyth::Mechanics {

    /**
     * @brief Physical particle with mass, position, velocity, and force accumulation
     */
    class Particle {
    public:
        Particle(const Quantity<Kilogram>& mass, const Vector3<Quantity<Meter>>& position) noexcept
            : mass_(mass), position_(position) {}

        virtual ~Particle() = default;

        /**
         * @brief Accumulates external force (cleared each integration step)
         */
        void ApplyForce(const Vector3<Quantity<Newton>>& force) {
            if (fixed_) return;
            external_force_ += force;
        }

        /**
         * @brief Integrates equations of motion using velocity Verlet.
         */
        virtual void Integrate(const Quantity<Second> dt) noexcept {
            if (fixed_) {
                external_force_ = Vector3<Quantity<Newton>>();
                return;
            }
            compute_forces_func_(this);

            const auto acc_old = external_force_ / mass_;
            position_ = position_ + velocity_ * dt + 0.5 * acc_old * dt * dt;

            external_force_ = Vector3<Quantity<Newton>>();
            compute_forces_func_(this);

            const auto acc_new = external_force_ / mass_;
            velocity_ = velocity_ + 0.5 * (acc_old + acc_new) * dt;

            external_force_ = Vector3<Quantity<Newton>>();
        }

        [[nodiscard]] constexpr const auto& GetPosition() const noexcept { return position_; }
        void SetPosition(const Vector3<Quantity<Meter>>& position) noexcept { position_ = position; }

        [[nodiscard]] constexpr auto GetMass() const noexcept { return mass_; }
        void SetMass(const Quantity<Kilogram>& mass) noexcept { mass_ = mass; }

        [[nodiscard]] constexpr const auto& GetVelocity() const noexcept { return velocity_; }
        void SetVelocity(const Vector3<Quantity<MeterPerSecond>>& velocity) noexcept { velocity_ = velocity; }

        [[nodiscard]] constexpr bool IsFixed() const noexcept { return fixed_; }
        void SetFixed(const bool fixed) noexcept { fixed_ = fixed; }

        [[nodiscard]] constexpr const auto& GetExternalForce() const noexcept { return external_force_; }
        void SetExternalForce(const Vector3<Quantity<Newton>>& external_force) noexcept { external_force_ = external_force; }

        void SetComputeForcesFunction(const std::function<void(Particle*)>& func) noexcept { compute_forces_func_ = func; }

        [[nodiscard]] Quantity<Joule> GetPotentialEnergy(const Quantity<Meter> y_ref = 0) const {
            return mass_ * Consts::g * (position_.y - y_ref);
        }

        [[nodiscard]] Quantity<Joule> GetKineticEnergy() const noexcept {
            return 0.5 * mass_ * velocity_.LengthSquared();
        }

        [[nodiscard]] Quantity<Joule> GetTotalEnergy(const Quantity<Meter> y_ref = 0) const noexcept {
            return GetPotentialEnergy(y_ref) + GetKineticEnergy();
        }

        friend std::ostream& operator<<(std::ostream& os, const Particle& p) {
            os << "Particle{pos=" << p.position_ << ", mass=" << p.mass_ << "}";
            return os;
        }

    protected:
        Quantity<Kilogram> mass_;
        Vector3<Quantity<MeterPerSecond>> velocity_{};
        bool fixed_{false};

        Vector3<Quantity<Newton>> external_force_{};
        Vector3<Quantity<Meter>> position_;

        std::function<void(Particle*)> compute_forces_func_ = [](Particle*) {};
    };
}

#endif //PHYTH_PARTICLE_H