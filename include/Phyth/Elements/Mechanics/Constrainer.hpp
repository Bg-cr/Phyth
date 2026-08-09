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
     * @brief Constrains a particle's distance from an anchor within an interval
     *
     * DistanceConstrainer ensures that a particle stays within a specified
     * distance range [min, max] from an anchor point. If the particle moves
     * outside this range, it is corrected back to the nearest bound.
     *
     * The correction is applied to both position and velocity to maintain
     * consistency with the time integration scheme.
     *
     * Example:
     *   auto anchor = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
     *   auto particle = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
     *
     *   // Keep distance between 0.5m and 1.5m
     *   Interval<Meter> range(0.5_m, 1.5_m);
     *   DistanceConstrainer constraint(anchor, particle, range);
     *
     *   // In simulation loop:
     *   constraint.Correct();
     */
    class DistanceConstrainer {
    public:
        /**
         * @brief Construct with automatic interval set to initial distance
         *
         * The interval's min and max are both set to the initial distance,
         * effectively fixing the distance to a constant value (rigid constraint).
         *
         * @param anchor The fixed reference point
         * @param particle The particle to constrain
         */
        DistanceConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle)
            : anchor_(std::move(anchor)), particle_(std::move(particle)),
              interval_((particle_->GetPosition() - anchor_->GetPosition()).Length()) {
            IntervalCheck();
        }

        /**
         * @brief Construct with a custom distance interval
         *
         * @param anchor The fixed reference point
         * @param particle The particle to constrain
         * @param interval The allowed distance range [min, max] (meters)
         *
         * @throws std::invalid_argument If interval.min < 0
         */
        DistanceConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                            const Interval<Meter> &interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), interval_(interval) {
            IntervalCheck();
        }

        /**
         * @brief Correct velocity to satisfy the distance constraint (Velocity Constrain)
         *
         * If the current distance is outside [min, max], the particle is moved
         * to the nearest bound. Velocity is adjusted to compensate for the
         * position change over the current time step.
         *
         * The velocity correction is:
         *   v_new = v_old - (position_error) / dt
         *
         * This ensures that the constraint is satisfied after integration.
         */
        void Correct() const {
            const Vector3<Quantity<Meter> > diff = particle_->GetPosition() - anchor_->GetPosition();
            const Quantity<Meter> length = diff.Length();
            if (interval_.Contains(length))
                return;

            const Vector3<Quantity<Meter> > target_offset = diff.Normalized() *
                                                            (length > interval_.GetMaximum()
                                                                 ? interval_.GetMaximum()
                                                                 : interval_.GetMinimum());

            particle_->SetVelocity(particle_->GetVelocity() - (diff - target_offset) / Config::dt);
        }

    private:
        std::shared_ptr<Particle> anchor_;
        std::shared_ptr<Particle> particle_;

        Interval<Meter> interval_;

        /**
         * @brief Validate that the interval minimum is non-negative
         *
         * @throws std::invalid_argument If interval.min < 0 (negative distance is physically meaningless)
         */
        void IntervalCheck() const {
            if (interval_.GetMinimum() < Quantity<Meter>(0)) {
                throw std::invalid_argument("DistanceConstrainer: interval minimum must be non-negative");
            }
        }
    };

    /**
     * @brief Constrains a particle's angular position around an anchor
     *
     * AngleConstrainer restricts the polar angle (theta) and azimuthal angle (phi)
     * of a particle relative to an anchor point. Both angles are constrained
     * within user-specified intervals.
     *
     * Angles are wrapped to [-2pi, 2pi] before validation, so intervals should
     * be specified in this range.
     *
     * Example:
     *   auto anchor = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
     *   auto particle = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
     *
     *   // Restrict theta to [30deg, 60deg], phi to [-45deg, 45deg]
     *   Interval<Radian> theta_range(30.0_deg, 60.0_deg);
     *   Interval<Radian> phi_range(-45.0_deg, 45.0_deg);
     *   AngleConstrainer constraint(anchor, particle, theta_range, phi_range);
     *
     *   // In simulation loop:
     *   constraint.Correct(dt);
     *
     * Note: Theta is the polar angle from the Z axis (0 to pi in standard
     *       spherical coordinates). Phi is the azimuthal angle from the X axis.
     */
    class AngleConstrainer {
    public:
        /**
         * @brief Construct with automatic intervals set to initial angles
         *
         * Both theta and phi intervals are fixed to the particle's initial angles.
         *
         * @param anchor The fixed reference point
         * @param particle The particle to constrain
         */
        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle)
            : anchor_(std::move(anchor)), particle_(std::move(particle)),
              theta_interval_(std::get<1>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())),
              phi_interval_(std::get<2>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())) {
            IntervalCheck();
        }

        /**
         * @brief Construct with custom theta interval, phi fixed to initial value
         *
         * @param anchor The fixed reference point
         * @param particle The particle to constrain
         * @param theta_interval Allowed range for polar angle (radians)
         *
         * @throws std::invalid_argument If interval bounds exceed [-2pi, 2pi]
         */
        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                         const Interval<Radian> &theta_interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), theta_interval_(theta_interval),
              phi_interval_(std::get<2>((particle_->GetPosition() - anchor_->GetPosition()).ToSpherical())) {
            IntervalCheck();
        }

        /**
         * @brief Construct with custom theta and phi intervals
         *
         * @param anchor The fixed reference point
         * @param particle The particle to constrain
         * @param theta_interval Allowed range for polar angle (radians)
         * @param phi_interval Allowed range for azimuthal angle (radians)
         *
         * @throws std::invalid_argument If interval bounds exceed [-2pi, 2pi]
         */
        AngleConstrainer(std::shared_ptr<Particle> anchor, std::shared_ptr<Particle> particle,
                         const Interval<Radian> &theta_interval, const Interval<Radian> &phi_interval)
            : anchor_(std::move(anchor)), particle_(std::move(particle)), theta_interval_(theta_interval),
              phi_interval_(phi_interval) {
            IntervalCheck();
        }

        /**
         * @brief Correct velocity to satisfy angular constraints (Velocity Constrain)
         *
         * The particle's position is adjusted to the nearest valid angle if it
         * falls outside the specified intervals. Both theta and phi are wrapped
         * to [-2pi, 2pi] before validation.
         *
         * The correction uses the particle's current distance from the anchor
         * (preserved during angular correction) and adjusts velocity accordingly.
         *
         * @param dt Time step used for velocity correction
         *
         * If the particle is exactly at the anchor (length == 0), no correction is applied.
         */
        void Correct(const Quantity<Second> dt) const {
            const Vector3<Quantity<Meter> > diff = particle_->GetPosition() - anchor_->GetPosition();
            const Quantity<Meter> length = diff.Length();

            if (length == 0_m) {
                return;
            }

            const auto spherical = diff.ToSpherical();
            Quantity<Radian> theta = std::get<1>(spherical);
            Quantity<Radian> phi = std::get<2>(spherical);

            // Wrap angles to [-2pi, 2pi] for validation
            theta = Utils::fmod(theta + Quantity<Radian>(Consts::tau), Quantity<Radian>(Consts::varpi)) - Quantity<
                        Radian>(Consts::tau);
            phi = Utils::fmod(phi + Quantity<Radian>(Consts::tau), Quantity<Radian>(Consts::varpi)) - Quantity<Radian>(
                      Consts::tau);

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

            const Vector3<Quantity<Meter> > target_pos = anchor_->GetPosition() + Vector3{
                                                             length * Utils::sin(target_theta) * Utils::cos(target_phi),
                                                             length * Utils::sin(target_theta) * Utils::sin(target_phi),
                                                             length * Utils::cos(target_theta)
                                                         };

            const Vector3<Quantity<Meter> > position_error = diff - (target_pos - anchor_->GetPosition());
            particle_->SetVelocity(particle_->GetVelocity() - position_error / dt);
        }

    private:
        std::shared_ptr<Particle> anchor_;
        std::shared_ptr<Particle> particle_;

        Interval<Radian> theta_interval_;
        Interval<Radian> phi_interval_;

        /**
         * @brief Validate that angle intervals are within [-2pi, 2pi]
         *
         * @throws std::invalid_argument If any interval bound exceeds [-2pi, 2pi]
         *
         * Angles are wrapped to this range before validation, so intervals must
         * be specified within this range.
         */
        void IntervalCheck() const {
            if (theta_interval_.GetMinimum() < Quantity<Radian>(-Consts::tau)) {
                throw std::invalid_argument("AngleConstrainer: theta interval minimum must be >= -2pi");
            }
            if (phi_interval_.GetMinimum() < Quantity<Radian>(-Consts::tau)) {
                throw std::invalid_argument("AngleConstrainer: phi interval minimum must be >= -2pi");
            }
            if (theta_interval_.GetMaximum() > Quantity<Radian>(Consts::tau)) {
                throw std::invalid_argument("AngleConstrainer: theta interval maximum must be <= 2pi");
            }
            if (phi_interval_.GetMaximum() > Quantity<Radian>(Consts::tau)) {
                throw std::invalid_argument("AngleConstrainer: phi interval maximum must be <= 2pi");
            }
        }
    };
}

#endif // PHYTH_CONSTRAINER_HPP
