#ifndef PHYTH_PARTICLE_H
#define PHYTH_PARTICLE_H

#include <functional>

#include "Phyth/Physical/PhysicalConsts.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Physical/PhysicalConfig.hpp"
#include "Phyth/Tools/Vector3.hpp"

namespace Phyth::Mechanics {
    /**
     * @brief Physical particle with mass, position, velocity, and force accumulation
     *
     * Particle is the fundamental building block for mechanical simulations.
     * It stores the state (position, velocity, mass) and accumulates forces
     * applied by constraints, springs, or user-defined force functions.
     *
     * Integration is performed using the velocity Verlet algorithm, which
     * is symplectic and energy-conserving for Hamiltonian systems.
     *
     * Particles can be marked as fixed (immovable), which is useful for
     * representing anchors or boundary conditions.
     *
     * Example:
     *   auto p = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
     *   p->SetComputeForcesFunction([](Particle* p) {
     *       // Apply gravity
     *       p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
     *   });
     *
     *   // In simulation loop:
     *   p->Integrate();
     *   auto total_energy = p->GetTotalEnergy();
     */
    class Particle {
    public:
        /**
         * @brief Construct a particle with mass and initial position
         *
         * @param mass Mass of the particle (kg)
         * @param position Initial position (m)
         *
         * Initial velocity is zero and the particle is not fixed.
         */
        Particle(const Quantity<Kilogram> &mass, const Vector3<Quantity<Meter> > &position) noexcept
            : mass_(mass), position_(position) {
        }

        virtual ~Particle() = default;

        /**
         * @brief Accumulate an external force on the particle
         *
         * Forces are accumulated during the force computation phase and
         * cleared after each integration step. This allows multiple force
         * sources (gravity, springs, collisions) to contribute.
         *
         * @param force Force to apply (N)
         *
         * Note: If the particle is fixed, this function does nothing.
         */
        void ApplyForce(const Vector3<Quantity<Newton> > &force) {
            if (fixed_) return;
            external_force_ += force;
        }

        /**
         * @brief Integrate equations of motion using velocity Verlet
         *
         * The velocity Verlet algorithm consists of three steps:
         *   1. Update position: x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
         *   2. Compute new acceleration: a(t+dt) = F(t+dt) / m
         *   3. Update velocity: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
         *
         * The user-defined force function is called twice per integration step:
         *   - Once before position update (using old forces)
         *   - Once after position update (using new position)
         *
         * External forces are automatically cleared at the start and end.
         * If the particle is fixed, this function does nothing.
         */
        virtual void Integrate() noexcept {
            if (fixed_) {
                external_force_ = Vector3<Quantity<Newton> >();
                return;
            }
            compute_forces_func_(this);

            const auto acc_old = external_force_ / mass_;
            position_ = position_ + velocity_ * Config::dt + 0.5 * acc_old * Config::dt * Config::dt;

            external_force_ = Vector3<Quantity<Newton> >();
            compute_forces_func_(this);

            const auto acc_new = external_force_ / mass_;
            velocity_ = velocity_ + 0.5 * (acc_old + acc_new) * Config::dt;

            external_force_ = Vector3<Quantity<Newton> >();
        }

        /** @return Current position (m) */
        [[nodiscard]] constexpr const auto &GetPosition() const noexcept { return position_; }

        /** @brief Set the particle position (m) */
        void SetPosition(const Vector3<Quantity<Meter> > &position) noexcept { position_ = position; }

        /** @return Mass (kg) */
        [[nodiscard]] constexpr auto GetMass() const noexcept { return mass_; }

        /** @brief Set the particle mass (kg) */
        void SetMass(const Quantity<Kilogram> &mass) noexcept { mass_ = mass; }

        /** @return Current velocity (m/s) */
        [[nodiscard]] constexpr const auto &GetVelocity() const noexcept { return velocity_; }

        /** @brief Set the particle velocity (m/s) */
        void SetVelocity(const Vector3<Quantity<MeterPerSecond> > &velocity) noexcept { velocity_ = velocity; }

        /** @return Whether the particle is fixed (immovable) */
        [[nodiscard]] constexpr bool IsFixed() const noexcept { return fixed_; }

        /** @brief Set whether the particle is fixed (immovable) */
        void SetFixed(const bool fixed) noexcept { fixed_ = fixed; }

        /** @return Accumulated external force (N) */
        [[nodiscard]] constexpr const auto &GetExternalForce() const noexcept { return external_force_; }

        /** @brief Manually set the external force (N) */
        void SetExternalForce(const Vector3<Quantity<Newton> > &external_force) noexcept {
            external_force_ = external_force;
        }

        /**
         * @brief Set the user-defined force computation function
         *
         * @param func Function that applies forces to the particle via ApplyForce() (Only Apply Force)
         *
         * The function is called automatically during Integrate().
         * The default function does nothing.
         *
         * Example:
         *   p->SetComputeForcesFunction([](Particle* p) {
         *       p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});  // gravity
         *   });
         */
        void SetComputeForcesFunction(const std::function<void(Particle *)> &func) noexcept {
            compute_forces_func_ = func;
        }

        /**
         * @brief Compute gravitational potential energy relative to a reference height
         *
         * @param y_ref Reference height (m), default 0
         * @return m * g * (position.y - y_ref) (J)
         *
         * The default reference height is y=0. Positive energy means the
         * particle is above the reference.
         */
        [[nodiscard]] Quantity<Joule> GetPotentialEnergy(const Quantity<Meter> y_ref = 0) const {
            return mass_ * Consts::g * (position_.y - y_ref);
        }

        /**
         * @brief Compute kinetic energy
         *
         * @return 0.5 * m * |v|^2 (J)
         */
        [[nodiscard]] Quantity<Joule> GetKineticEnergy() const noexcept {
            return 0.5 * mass_ * velocity_.LengthSquared();
        }

        /**
         * @brief Compute total mechanical energy
         *
         * @param y_ref Reference height for potential energy, default 0
         * @return PotentialEnergy + KineticEnergy (J)
         *
         * For conservative systems, this quantity should remain constant.
         * Energy drift can be used as a measure of numerical accuracy.
         */
        [[nodiscard]] Quantity<Joule> GetTotalEnergy(const Quantity<Meter> y_ref = 0) const noexcept {
            return GetPotentialEnergy(y_ref) + GetKineticEnergy();
        }

        /**
         * @brief Stream insertion operator for Particle
         *
         * Output format: "Particle{pos=(x, y, z), mass=value}"
         */
        friend std::ostream &operator<<(std::ostream &os, const Particle &p) {
            os << "Particle{pos=" << p.position_ << ", mass=" << p.mass_ << "}";
            return os;
        }

    protected:
        Quantity<Kilogram> mass_;
        Vector3<Quantity<MeterPerSecond> > velocity_{};
        bool fixed_{false};

        Vector3<Quantity<Newton> > external_force_{};
        Vector3<Quantity<Meter> > position_;

        std::function<void(Particle *)> compute_forces_func_ = [](Particle *) {
        };
    };
}

#endif // PHYTH_PARTICLE_H
