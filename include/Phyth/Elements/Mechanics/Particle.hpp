#ifndef PHYTH_PARTICLE_H
#define PHYTH_PARTICLE_H

#include "Phyth/Physical/PhysicalConsts.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Vector3.hpp"

namespace Phyth::Mechanics {
    class Particle {
    public:
        Particle(const Quantity<Kilogram> &mass, const Vector3<Quantity<Meter>> &position) noexcept
            : mass_(mass), position_(position) {}

        void ApplyForce(const Vector3<Quantity<Newton>> &force) {
            if (fixed_) return;
            external_force_ += force;
        }

        void Integrate(const Quantity<Second>& dt) noexcept {
            if (fixed_) {
                external_force_ = Vector3<Quantity<Newton>>();
                return;
            }

            const auto acc = external_force_ / mass_;
            position_ = position_ + velocity_ * dt + 0.5 * acc * dt * dt;
            velocity_ = velocity_ + acc * dt;

            external_force_ = Vector3<Quantity<Newton>>();
        }

        [[nodiscard]] constexpr const auto &GetPosition() const noexcept { return position_; }
        void SetPosition(const Vector3<Quantity<Meter>> &position) noexcept { position_ = position; }

        [[nodiscard]] constexpr auto GetMass() const noexcept { return mass_; }
        void SetMass(const Quantity<Kilogram> &mass) noexcept { mass_ = mass; }

        [[nodiscard]] constexpr const auto &GetVelocity() const noexcept { return velocity_; }
        void SetVelocity(const Vector3<Quantity<MeterPerSecond>> &velocity) noexcept { velocity_ = velocity; }

        [[nodiscard]] constexpr bool IsFixed() const noexcept { return fixed_; }
        void SetFixed(const bool fixed) noexcept { fixed_ = fixed; }

        [[nodiscard]] constexpr const auto &GetExternalForce() const noexcept { return external_force_; }
        void SetExternalForce(const Vector3<Quantity<Newton>> &external_force) noexcept { external_force_ = external_force; }

        [[nodiscard]] Quantity<Joule> GetPotentialEnergy(const Quantity<Meter> y_ref = 0) const {
            return mass_ * Consts::g * (position_.y - y_ref);
        }

        [[nodiscard]] Quantity<Joule> GetKineticEnergy() const noexcept {
            return 0.5 * mass_ * velocity_.LengthSquared();
        }

        [[nodiscard]] Quantity<Joule> GetTotalEnergy(const Quantity<Meter> y_ref = 0) const noexcept {
            return GetPotentialEnergy(y_ref) + GetKineticEnergy();
        }

        friend std::ostream &operator<<(std::ostream &os, const Particle &p) {
            os << "Particle{pos=" << p.position_ << ", mass=" << p.mass_ << "}";
            return os;
        }

    protected:
        Quantity<Kilogram> mass_;
        Vector3<Quantity<MeterPerSecond>> velocity_ {};
        bool fixed_ { false };

        Vector3<Quantity<Newton>> external_force_ {};
        Vector3<Quantity<Meter>> position_;
    };
}

#endif //PHYTH_PARTICLE_H