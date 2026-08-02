#ifndef PHYTH_CONSTRAINER_HPP
#define PHYTH_CONSTRAINER_HPP

#include <memory>
#include <utility>
#include <cmath>

#include "Particle.hpp"
#include "Phyth/Core/Interval.hpp"
#include "Phyth/Core/Unit.hpp"

namespace Phyth::Mechanics {

    /**
     * @brief Constrains a particle's distance from an anchor within an interval by corrct velocity
     *
     * If the particle moves outside [min, max], corrects its position
     * back to the nearest bound and adjusts velocity accordingly.
     */
    class DistanceConstrainer {
    public:
        DistanceConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle)
            : anchor_(std::move(anchor)), particle_(std::move(particle)),
              interval_((particle_->GetPosition() - anchor_->GetPosition()).Length()) {
            IntervalCheck();
        }

        DistanceConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                            const Interval<Meter>& interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), interval_(interval) {
            IntervalCheck();
        }

        /**
         * @brief Corrects position and velocity to satisfy the distance constraint
         * @param dt Time step used for velocity correction
         */
        void Correct(const Quantity<Second> dt) const {
            const Vector3<Quantity<Meter>> diff = particle_->GetPosition() - anchor_->GetPosition();
            const Quantity<Meter> length = diff.Length();
            if (interval_.Contains(length))
                return;

            const Vector3<Quantity<Meter>> target_offset = diff.Normalized() *
                                                            (length > interval_.GetMaximum()
                                                                 ? interval_.GetMaximum()
                                                                 : interval_.GetMinimum());

            particle_->SetVelocity(particle_->GetVelocity() - (diff - target_offset) / dt);
        }

    private:
        std::shared_ptr<Particle> anchor_;
        std::shared_ptr<Particle> particle_;

        Interval<Meter> interval_;

        void IntervalCheck() const {
            if (interval_.GetMinimum() < Quantity<Meter>(0)) {
                throw std::invalid_argument("The minimum of interval must be non-negative");
            }
        }
    };

    /**
     * @brief Constrains a particle's angular position around an anchor
     *
     * Restricts theta (polar angle from Z axis) and phi (azimuthal angle from X axis)
     * within specified intervals. Corrects position back to nearest valid angle.
     * Angles are wrapped to [-2pi, 2pi] before validation.
     */
    class AngleConstrainer {
    public:
        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle)
            : anchor_(std::move(anchor)), particle_(std::move(particle)),
              theta_interval_(std::get<1>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())),
              phi_interval_(std::get<2>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())) {
            IntervalCheck();
        }

        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                         const Interval<Radian>& theta_interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), theta_interval_(theta_interval),
              phi_interval_(std::get<2>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())) {
            IntervalCheck();
        }

        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                         const Interval<Radian>& theta_interval, const Interval<Radian>& phi_interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), theta_interval_(theta_interval),
              phi_interval_(phi_interval) {
            IntervalCheck();
        }

        /**
         * @brief Corrects position and velocity to satisfy angular constraints
         * @param dt Time step used for velocity correction
         */
        void Correct(const Quantity<Second> dt) const {
            const Vector3<Quantity<Meter>> diff = particle_->GetPosition() - anchor_->GetPosition();
            const Quantity<Meter> length = diff.Length();

            if (length == 0_m) {
                return;
            }

            const auto spherical = diff.ToSpherical();
            Quantity<Radian> theta = std::get<1>(spherical);
            Quantity<Radian> phi = std::get<2>(spherical);

            theta = Utils::fmod(theta + Quantity<Radian>(Consts::tau), Quantity<Radian>(Consts::varpi)) - Quantity<Radian>(Consts::tau);
            phi = Utils::fmod(phi + Quantity<Radian>(Consts::tau), Quantity<Radian>(Consts::varpi)) - Quantity<Radian>(Consts::tau);

            const Quantity<Radian> target_theta = theta_interval_.Contains(theta)
                                            ? theta
                                            : theta > theta_interval_.GetMaximum()
                                                  ? theta_interval_.GetMaximum()
                                                  : theta_interval_.GetMinimum();

            const Quantity<Radian> target_phi = phi_interval_.Contains(phi)
                                           ? phi
                                           : phi > phi_interval_.GetMaximum()
                                                 ? phi_interval_.GetMaximum()
                                                 : phi_interval_.GetMinimum();

            if (theta == target_theta && phi == target_phi) {
                return;
            }

            const Vector3<Quantity<Meter>> target_pos = anchor_->GetPosition() + Vector3{
                length * Utils::sin(target_theta) * Utils::cos(target_phi),
                length * Utils::sin(target_theta) * Utils::sin(target_phi),
                length * Utils::cos(target_theta)
            };

            const Vector3<Quantity<Meter>> position_error = diff - (target_pos - anchor_->GetPosition());
            particle_->SetVelocity(particle_->GetVelocity() - position_error / dt);
        }

    private:
        std::shared_ptr<Particle> anchor_;
        std::shared_ptr<Particle> particle_;

        Interval<Radian> theta_interval_;
        Interval<Radian> phi_interval_;

        void IntervalCheck() const {
            if (theta_interval_.GetMinimum() < Quantity<Radian>(-Consts::tau)) {
                throw std::invalid_argument("The minimum of theta interval must be greater than 2pi.");
            }
            if (phi_interval_.GetMinimum() < Quantity<Radian>(-Consts::tau)) {
                throw std::invalid_argument("The minimum of phi interval must be greater than 2pi.");
            }
            if (theta_interval_.GetMaximum() > Quantity<Radian>(Consts::tau)) {
                throw std::invalid_argument("The maximum of theta interval must be less than 2pi.");
            }
            if (phi_interval_.GetMaximum() > Quantity<Radian>(Consts::tau)) {
                throw std::invalid_argument("The maximum of phi interval must be less than 2pi.");
            }
        }
    };

}

#endif //PHYTH_CONSTRAINER_HPP