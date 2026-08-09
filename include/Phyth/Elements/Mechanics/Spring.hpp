#ifndef PHYTH_SPRING_H
#define PHYTH_SPRING_H

#include <stdexcept>
#include <cmath>
#include <memory>
#include <utility>

#include "Particle.hpp"

namespace Phyth::Mechanics {
    /**
     * @brief Spring constraint between two particles
     *
     * Spring applies a Hooke's law force between two particles:
     *   F = -k * (current_length - rest_length) * direction
     *
     * Additionally, a damping force proportional to the relative velocity
     * along the spring axis can be added:
     *   F_damp = -c * (v1 - v2) · direction
     *
     * The spring is bi-directional: equal and opposite forces are applied
     * to both particles. If one particle is fixed, the other still moves
     * correctly.
     *
     * Example:
     *   auto p1 = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
     *   auto p2 = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
     *
     *   Spring spring(p1, p2, 10.0_N/m, 0.5_kgps);
     *
     *   // In simulation loop:
     *   spring.Constrain();
     *
     *   auto energy = spring.GetElasticEnergy();  // 1/2 * k * elongation^2
     */
    class Spring {
    public:
        /**
         * @brief Construct a spring between two particles
         *
         * @param p1 First particle (force applied in direction p2->p1)
         * @param p2 Second particle (force applied in direction p1->p2)
         * @param stiffness Spring constant k (N/m)
         * @param damping Damping coefficient c (kg/s), default 0
         *
         * The rest length is automatically set to the initial distance
         * between the two particles.
         *
         * @throws std::invalid_argument If rest length <= 0, stiffness < 0,
         *         or damping < 0
         */
        Spring(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2,
               const Quantity<NewtonPerMeter> &stiffness,
               const Quantity<KilogramPerSecond> &damping = Quantity<KilogramPerSecond>{0})
            : p1_(std::move(p1)), p2_(std::move(p2)),
              stiffness_(stiffness),
              damping_(damping),
              rest_length_((p1_->GetPosition() - p2_->GetPosition()).Length()) {
            validateRestLength(rest_length_);
            validateStiffness(stiffness_);
            validateDamping(damping_);
        }

        /**
         * @brief Apply spring and damping forces to both particles (Force Constrain)
         *
         * The spring force is:
         *   F = k * (current_length - rest_length) * direction
         *
         * If damping > 0, an additional damping force is applied:
         *   F_damp = c * (v1 - v2) · direction
         *
         * Both forces are applied in opposite directions to p1 and p2.
         */
        void Constrain() const {
            const Vector3<Quantity<Meter> > displacement = p2_->GetPosition() - p1_->GetPosition();
            const Quantity<Meter> current_length = displacement.Length();
            const Quantity<Meter> elongation = current_length - rest_length_;

            const Quantity<Newton> force_magnitude = stiffness_ * elongation;
            const Vector3<Scalar> direction = displacement.Normalized();

            const Vector3<Quantity<Newton> > spring_force1 = direction * force_magnitude;
            const Vector3<Quantity<Newton> > spring_force2 = -direction * force_magnitude;

            if (damping_ > 0.0_kgps) {
                const Vector3<Quantity<MeterPerSecond> > relative_velocity = p2_->GetVelocity() - p1_->GetVelocity();
                const Quantity<Newton> damping_force_magnitude = damping_ * relative_velocity.Dot(direction);
                const Vector3<Quantity<Newton> > damping_force1 = direction * damping_force_magnitude;
                const Vector3<Quantity<Newton> > damping_force2 = -direction * damping_force_magnitude;
                p1_->ApplyForce(damping_force1);
                p2_->ApplyForce(damping_force2);
            }

            p1_->ApplyForce(spring_force1);
            p2_->ApplyForce(spring_force2);
        }

        /** @return Shared pointer to the first particle */
        [[nodiscard]] const std::shared_ptr<Particle> &GetParticle1() const { return p1_; }

        /** @return Shared pointer to the second particle */
        [[nodiscard]] const std::shared_ptr<Particle> &GetParticle2() const { return p2_; }

        /** @return Current rest length (constant unless modified) */
        [[nodiscard]] const Quantity<Meter> &GetRestLength() const { return rest_length_; }

        /** @return Spring stiffness constant k (N/m) */
        [[nodiscard]] const Quantity<NewtonPerMeter> &GetStiffness() const { return stiffness_; }

        /** @return Damping coefficient c (kg/s) */
        [[nodiscard]] const Quantity<KilogramPerSecond> &GetDamping() const { return damping_; }

        /** @return Current distance between the two particles */
        [[nodiscard]] Quantity<Meter> GetCurrentLength() const {
            return (p1_->GetPosition() - p2_->GetPosition()).Length();
        }

        /** @return Current elongation (current_length - rest_length) */
        [[nodiscard]] Quantity<Meter> GetElongation() const { return GetCurrentLength() - rest_length_; }

        /**
         * @brief Set a new rest length
         *
         * @param rest_length New rest length (must be positive finite)
         * @throws std::invalid_argument If rest_length <= 0
         */
        void SetRestLength(const Quantity<Meter> &rest_length) {
            validateRestLength(rest_length);
            rest_length_ = rest_length;
        }

        /**
         * @brief Set a new stiffness
         *
         * @param stiffness New stiffness (must be non-negative)
         * @throws std::invalid_argument If stiffness < 0
         */
        void SetStiffness(const Quantity<NewtonPerMeter> &stiffness) {
            validateStiffness(stiffness);
            stiffness_ = stiffness;
        }

        /**
         * @brief Set a new damping coefficient
         *
         * @param damping New damping coefficient (must be non-negative)
         * @throws std::invalid_argument If damping < 0
         */
        void SetDamping(const Quantity<KilogramPerSecond> &damping) {
            validateDamping(damping);
            damping_ = damping;
        }

        /**
         * @brief Compute the elastic potential energy stored in the spring
         *
         * @return 1/2 * k * elongation^2 (Joules)
         */
        [[nodiscard]] Quantity<Joule> GetElasticEnergy() const {
            const auto elongation = GetElongation();
            return 0.5 * stiffness_ * elongation * elongation;
        }

        /**
         * @brief Compute energy per unit length
         *
         * @return ElasticEnergy / CurrentLength (Newtons)
         *
         * Note: If current_length is zero, this will produce a division by zero.
         *       Caller should ensure the particles are separated.
         */
        [[nodiscard]] Quantity<Newton> GetEnergyDensity() const {
            return GetElasticEnergy() / GetCurrentLength();
        }

    private:
        std::shared_ptr<Particle> p1_;
        std::shared_ptr<Particle> p2_;

        Quantity<NewtonPerMeter> stiffness_;
        Quantity<KilogramPerSecond> damping_;
        Quantity<Meter> rest_length_;

        static void validateRestLength(const Quantity<Meter> &length) {
            if (const auto len = length.value; !std::isfinite(len) || len <= 0.0) {
                throw std::invalid_argument("Spring::validateRestLength: rest length must be a positive finite number");
            }
        }

        static void validateStiffness(const Quantity<NewtonPerMeter> &stiffness) {
            if (stiffness.value < 0.0) {
                throw std::invalid_argument("Spring::validateStiffness: stiffness must be non-negative");
            }
        }

        static void validateDamping(const Quantity<KilogramPerSecond> &damping) {
            if (damping.value < 0.0) {
                throw std::invalid_argument("Spring::validateDamping: damping must be non-negative");
            }
        }
    };
}

#endif // PHYTH_SPRING_H
