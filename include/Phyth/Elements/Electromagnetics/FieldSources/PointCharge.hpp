#ifndef PHYTH_POINT_CHARGE_HPP
#define PHYTH_POINT_CHARGE_HPP

#include "ChargeSource.hpp"
#include "MagneticSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"
#include "Phyth/Tools/TimeHistory.hpp"

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
        explicit PointCharge(const Quantity<Kilogram> mass,
                             const Vector3<Quantity<Meter> > &position,
                             const Quantity<Coulomb> charge_value)
            : Particle(mass, position),
              charge_value_(charge_value) {
            position_history_.Register(position_);
            velocity_history_.Register({});
            acceleration_history_.Register({});
        }

        void SetComputeForcesFunction(const std::function<void(PointCharge *)> &func) noexcept {
            Particle::SetComputeForcesFunction([func](Particle *p) {
                func(dynamic_cast<PointCharge *>(p));
            });
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return electric field obtained from Liénard-Wiechert formulas
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto r = point - position_;
            if (r.Length() < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at charge location");
            }

            if (fixed_) {
                return Consts::k_E * charge_value_ / r.LengthSquared() * r.Normalized();
            }
            Quantity<Second> retarded_dt;
            try {
                retarded_dt = SolveRetardedTimeOffset(point);
            } catch (std::runtime_error &) {
                Vector3<Quantity<Meter>> retarded_r = r;
                if (!position_history_.IsEmpty()) {
                    retarded_r = position_history_.GetValueByTime(0_s);
                    // The observation point has not yet seen the movement of charges.
                    // It is an electrostatic field of charges at p_0.
                }
                return Consts::k_E * charge_value_ / retarded_r.LengthSquared() * retarded_r.Normalized();
            }

            const auto retarded_pos = GetPositionAtOffset(retarded_dt);
            const auto retarded_v = GetVelocityAtOffset(retarded_dt);
            const auto retarded_a = GetAccelerationAtOffset(retarded_dt);

            const auto R = point - retarded_pos;
            const auto hat_R = R.Normalized();
            const auto beta = retarded_v / Consts::c;
            const auto beta_a = retarded_a / Consts::c;
            const auto gamma = 1 / Utils::sqrt(1_ - Utils::square(beta.Length()));

            const auto velocity_field = Consts::k_E * charge_value_ * (hat_R - beta)
                                        / Utils::square(gamma * R.Length())
                                        / Utils::cube(1_ - hat_R.Dot(beta));

            const auto acceleration_field = Consts::k_E * charge_value_
                                            * hat_R.Cross((hat_R - beta).Cross(beta_a))
                                            / Consts::c
                                            / Utils::cube(1_ - hat_R.Dot(beta))
                                            / R.Length();
            return velocity_field + acceleration_field;
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return magnetic field obtained from Liénard-Wiechert formulas
         */
        [[nodiscard]] Vector3<Quantity<Tesla> >
        GetMagneticFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            if (fixed_) {
                return {};
            }
            const auto r = point - position_;
            Quantity<Second> retarded_dt;
            try {
                retarded_dt = SolveRetardedTimeOffset(point);
            } catch (std::runtime_error &) {
                return Consts::mu_0 / Consts::varpi * charge_value_ * velocity_.Cross(r.Normalized()) / r.
                       LengthSquared();
            }
            const auto R = point - GetPositionAtOffset(retarded_dt);
            const auto hat_R = R.Normalized();

            return (hat_R / Consts::c).Cross(GetElectricFieldAt(point));
        }

        /**
         * @throw std::runtime_error when evaluated at distance < Config::epsilon
         *
         * @return electric potential obtained from Liénard-Wiechert formulas
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();

            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at charge location");
            }

            if (fixed_ || position_history_.IsEmpty()) {
                return Consts::k_E * charge_value_ / dist;
            }

            Quantity<Second> retarded_dt;
            try {
                retarded_dt = SolveRetardedTimeOffset(point);
            } catch (std::runtime_error &) {
                return Consts::k_E * charge_value_ / dist;
            }

            const auto retarded_pos = GetPositionAtOffset(retarded_dt);
            const auto retarded_v = GetVelocityAtOffset(retarded_dt);

            const auto R_vec = point - retarded_pos;
            const auto R = R_vec.Length();
            const auto hat_R = R_vec.Normalized();
            const auto beta = retarded_v / Consts::c;

            const auto denominator = R - hat_R.Dot(beta) * R;

            if (denominator < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Line-of-sight aligns with charge velocity at light speed");
            }

            return Consts::k_E * charge_value_ / denominator;
        }

        void ApplyElectricFieldForce(const Vector3<Quantity<NewtonPerCoulomb> > &field) {
            ApplyForce(charge_value_ * field);
        }

        void ApplyMagneticFieldForce(const Vector3<Quantity<Tesla> > &field) {
            ApplyForce(charge_value_ * velocity_.Cross(field));
        }

        void ApplyLorentzForce(const Vector3<Quantity<NewtonPerCoulomb> > &electric_field,
                               const Vector3<Quantity<Tesla> > &magnetic_field) {
            ApplyElectricFieldForce(electric_field);
            ApplyMagneticFieldForce(magnetic_field);
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return charge_value_;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            return (point - position_).Length();
        }

        void Integrate() noexcept override {
            Particle::Integrate();

            position_history_.Register(position_);
            velocity_history_.Register(velocity_);
            acceleration_history_.Register(external_force_ / mass_);
        }

    private:
        Quantity<Coulomb> charge_value_;

        TimeHistory<Vector3<Quantity<Meter> > > position_history_{};
        TimeHistory<Vector3<Quantity<MeterPerSecond> > > velocity_history_{};
        TimeHistory<Vector3<Quantity<MeterPerSecondSquared> > > acceleration_history_{};

        [[nodiscard]] Vector3<Quantity<Meter> >
        GetPositionAtOffset(const Quantity<Second> &offset) const {
            if (offset == 0_s) {
                return position_;
            }
            return position_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<MeterPerSecond> >
        GetVelocityAtOffset(const Quantity<Second> &offset) const {
            if (offset == 0_s) {
                return velocity_;
            }
            return velocity_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<MeterPerSecondSquared> >
        GetAccelerationAtOffset(const Quantity<Second> &offset) const {
            if (offset == 0_s) {
                return external_force_ / mass_;
            }
            return acceleration_history_.GetValueByOffset(offset);
        }


        [[nodiscard]] Quantity<Second>
        SolveRetardedTimeOffset(const Vector3<Quantity<Meter> > &point) const {
            if (position_history_.IsEmpty()) {
                return 0_s;
            }

            const Quantity<Meter> R0 = (point - position_).Length();
            Quantity<Second> dt_delay = R0 / Consts::c;

            const Quantity<Second> max_dt = static_cast<double>(position_history_.GetSize())
                                            * Config::dt;
            dt_delay = std::clamp(dt_delay, Quantity<Second>(0.0), max_dt);

            for (int iter = 0; iter < Config::max_iterations; ++iter) {
                const Vector3<Quantity<Meter> > pos_r = GetPositionAtOffset(dt_delay);
                const Vector3<Quantity<MeterPerSecond> > v_r = GetVelocityAtOffset(dt_delay);

                const Vector3<Quantity<Meter> > R_vec = point - pos_r;
                const Quantity<Meter> R = R_vec.Length();
                const Vector3<Scalar> hat_R = R_vec.Normalized();

                const Quantity<Meter> f = R - Consts::c * dt_delay;

                if (Utils::abs(f / Consts::c) < Quantity<Second>(Config::tolerance)) {
                    return dt_delay;
                }

                const Quantity<MeterPerSecond> df = -hat_R.Dot(v_r) - Consts::c;

                if (Utils::abs(df) < Quantity<MeterPerSecond>(Config::epsilon)) {
                    break;
                }

                const Quantity<Second> delta = f / df;
                dt_delay = dt_delay - delta;
                dt_delay = std::clamp(dt_delay, Quantity<Second>(0.0), max_dt);
            }

            return SolveRetardedTimeOffsetBisection(point, Quantity<Second>(0.0), max_dt);
        }

        [[nodiscard]] Quantity<Second>
        SolveRetardedTimeOffsetBisection(const Vector3<Quantity<Meter> > &point,
                                         Quantity<Second> dt_low,
                                         Quantity<Second> dt_high) const {
            auto ComputeResidual = [&](const Quantity<Second> dt) -> Quantity<Meter> {
                const Vector3<Quantity<Meter> > pos_r = GetPositionAtOffset(dt);
                const Quantity<Meter> R = (point - pos_r).Length();
                return R - Consts::c * dt;
            };

            Quantity<Meter> f_low = ComputeResidual(dt_low);
            Quantity<Meter> f_high = ComputeResidual(dt_high);

            if (f_low * f_high > 0_m2) {
                throw std::runtime_error("Signal has not reached observation point or no solution in interval");
            }

            for (int iter = 0; iter < Config::max_iterations; ++iter) {
                const Quantity<Second> dt_mid = (dt_low + dt_high) / 2;
                const Quantity<Meter> f_mid = ComputeResidual(dt_mid);

                if (Utils::abs(f_mid / Consts::c) < Quantity<Second>(Config::tolerance)) {
                    return dt_mid;
                }

                if (f_low * f_mid < 0_m2) {
                    dt_high = dt_mid;
                    f_high = f_mid;
                } else {
                    dt_low = dt_mid;
                    f_low = f_mid;
                }
            }

            throw std::runtime_error("Bisection failed to converge");
        }
    };
}

#endif //PHYTH_POINT_CHARGE_HPP