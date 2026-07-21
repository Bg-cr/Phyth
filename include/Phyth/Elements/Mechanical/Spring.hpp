#ifndef PHYTH_SPRING_H
#define PHYTH_SPRING_H

#include <stdexcept>
#include <cmath>
#include <memory>
#include <utility>

#include "Particle.hpp"

namespace Phyth::Mechanical {
    class Spring {
    public:
        Spring(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2,
               const Quantity<NewtonPerMeter> &stiffness,
               const Quantity<KilogramPerSecond> &damping = Quantity<KilogramPerSecond>{0})
            : p1_(std::move(p1)), p2_(std::move(p2)),
              stiffness_(stiffness),
              damping_(damping),
              rest_length_((p1_->GetPosition() - p2_->GetPosition()).Length()) {
            validate_rest_length(rest_length_);
            validate_stiffness(stiffness_);
            validate_damping(damping_);
        }

        void Constrain() const {
            const Vector3<Quantity<Meter> > pos1 = p1_->GetPosition();
            const Vector3<Quantity<Meter> > pos2 = p2_->GetPosition();
            const Vector3<Quantity<MeterPerSecond> > vel1 = p1_->GetVelocity();
            const Vector3<Quantity<MeterPerSecond> > vel2 = p2_->GetVelocity();

            const Vector3<Quantity<Meter> > displacement = pos2 - pos1;
            const Quantity<Meter> current_length = displacement.Length();
            const Quantity<Meter> elongation = current_length - rest_length_;

            const Quantity<Newton> force_magnitude = stiffness_ * elongation;
            const Vector3<Scalar> direction = displacement.Normalized();

            const Vector3<Quantity<Newton> > spring_force1 = direction * force_magnitude;
            const Vector3<Quantity<Newton> > spring_force2 = -direction * force_magnitude;

            if (damping_ > 0.0_kgps) {
                const Vector3<Quantity<MeterPerSecond> > relative_velocity = vel2 - vel1;
                const Quantity<Newton> damping_force_magnitude = damping_ * relative_velocity.Dot(direction);
                const Vector3<Quantity<Newton> > damping_force1 = direction * damping_force_magnitude;
                const Vector3<Quantity<Newton> > damping_force2 = -direction * damping_force_magnitude;
                p1_->ApplyForce(damping_force1);
                p2_->ApplyForce(damping_force2);
            }

            p1_->ApplyForce(spring_force1);
            p2_->ApplyForce(spring_force2);
        }

        [[nodiscard]] const std::shared_ptr<Particle> &GetParticle1() const { return p1_; }
        [[nodiscard]] const std::shared_ptr<Particle> &GetParticle2() const { return p2_; }

        [[nodiscard]] const Quantity<Meter> &GetRestLength() const { return rest_length_; }
        [[nodiscard]] const Quantity<NewtonPerMeter> &GetStiffness() const { return stiffness_; }
        [[nodiscard]] const Quantity<KilogramPerSecond> &GetDamping() const { return damping_; }

        [[nodiscard]] Quantity<Meter> GetCurrentLength() const { return (p1_->GetPosition() - p2_->GetPosition()).Length(); }
        [[nodiscard]] Quantity<Meter> GetElongation() const { return GetCurrentLength() - rest_length_; }

        void SetRestLength(const Quantity<Meter> &rest_length) {
            validate_rest_length(rest_length);
            rest_length_ = rest_length;
        }

        void SetStiffness(const Quantity<NewtonPerMeter> &stiffness) {
            validate_stiffness(stiffness);
            stiffness_ = stiffness;
        }

        void SetDamping(const Quantity<KilogramPerSecond> &damping) {
            validate_damping(damping);
            damping_ = damping;
        }

    private:
        std::shared_ptr<Particle> p1_;
        std::shared_ptr<Particle> p2_;

        Quantity<NewtonPerMeter> stiffness_;
        Quantity<KilogramPerSecond> damping_;
        Quantity<Meter> rest_length_;

        static void validate_rest_length(const Quantity<Meter> &length) {
            if (const auto len = length.value; !std::isfinite(len) || len <= 0.0) {
                throw std::invalid_argument("Rest length must be a positive finite number");
            }
        }

        static void validate_stiffness(const Quantity<NewtonPerMeter> &stiffness) {
            if (stiffness.value < 0.0) {
                throw std::invalid_argument("Stiffness must be non-negative");
            }
        }

        static void validate_damping(const Quantity<KilogramPerSecond> &damping) {
            if (damping.value < 0.0) {
                throw std::invalid_argument("Damping must be non-negative");
            }
        }
    };
}

#endif //PHYTH_SPRING_H
