#ifndef PHYTH_MAGNETIC_DIPOLE_HPP
#define PHYTH_MAGNETIC_DIPOLE_HPP

#include "MagneticSource.hpp"
#include "ChargeSource.hpp"
#include "Phyth/Elements/Mechanics/Particle.hpp"
#include "Phyth/Tools/Quaternion.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"
#include "Phyth/Tools/TimeHistory.hpp"

namespace Phyth::Electromagnetics {

    class MagneticDipole : public Mechanics::Particle,
                           public MagneticSource,
                           public ChargeSource {
    public:
        explicit MagneticDipole(
            const Quantity<Kilogram> mass,
            const Vector3<Quantity<Meter>>& position,
            const Vector3<Quantity<AmpereMeterSquared>>& moment,
            const Quantity<KilogramMeterSquared>& moment_of_inertia = 0_kgm2
        ) : Particle(mass, position),
            moment_(moment),
            initial_moment_(moment),
            moment_of_inertia_(moment_of_inertia) {
            position_history_.Register(position_);
            velocity_history_.Register(velocity_);
            acceleration_history_.Register(external_force_ / mass_);
            moment_history_.Register(moment_);

            const auto p0 = ComputeRelativisticDipoleMoment(velocity_, moment_);
            electric_dipole_history_.Register(p0);
        }

        void SetComputeForcesFunction(const std::function<void(MagneticDipole *)> &func) noexcept {
            Particle::SetComputeForcesFunction([func](Particle *p) {
                func(dynamic_cast<MagneticDipole *>(p));
            });
        }

        [[nodiscard]] Vector3<Quantity<Tesla>>
        GetMagneticFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();

            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at dipole location");
            }

            if (fixed_ || position_history_.GetSize() < 4) {
                const auto hat_r = r / dist;
                constexpr auto coeff = Consts::mu_0 / (4 * Consts::pi);
                return coeff * (3 * moment_.Dot(hat_r) * hat_r - moment_) / Utils::cube(dist);
            }

            try {
                const Quantity<Second> retarded_dt = SolveRetardedTimeOffset(point);
                const auto retarded_pos = GetPositionAtOffset(retarded_dt);
                const auto retarded_m = GetMomentAtOffset(retarded_dt);
                const auto retarded_m_dot = ComputeMomentFirstDerivative(retarded_dt);
                const auto retarded_m_ddot = ComputeMomentSecondDerivative(retarded_dt);

                const auto R = point - retarded_pos;
                const auto R_mag = R.Length();
                if (R_mag < Quantity<Meter>(Config::epsilon)) {
                    throw std::runtime_error("Field diverges at dipole location");
                }
                const auto hat_R = R / R_mag;
                constexpr auto coeff = Consts::mu_0 / (4 * Consts::pi);

                const auto term1 = coeff * (3 * retarded_m.Dot(hat_R) * hat_R - retarded_m) / Utils::cube(R_mag);
                const auto term2 = coeff * (3 * retarded_m_dot.Dot(hat_R) * hat_R - retarded_m_dot) /
                                   (Consts::c * R_mag * R_mag);
                const auto term3 = coeff * retarded_m_ddot.Cross(hat_R).Cross(hat_R) / (Consts::c * Consts::c * R_mag);

                return term1 + term2 + term3;
            } catch (const std::runtime_error&) {
                const auto hat_r = r / dist;
                constexpr auto coeff = Consts::mu_0 / (4 * Consts::pi);
                return coeff * (3 * moment_.Dot(hat_r) * hat_r - moment_) / Utils::cube(dist);
            }
        }

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            const auto r = point - position_;
            const auto dist = r.Length();

            if (dist < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Field diverges at dipole location");
            }

            if (fixed_ || position_history_.GetSize() < 4 || velocity_.Length() < Quantity<MeterPerSecond>(Config::epsilon)) {
                const auto p = GetRelativisticDipoleMomentCurrent();
                const auto hat_r = r.Normalized();
                return Consts::k_E * (3 * p.Dot(hat_r) * hat_r - p) / Utils::cube(dist);
            }

            try {
                const Quantity<Second> retarded_dt = SolveRetardedTimeOffset(point);
                const auto retarded_pos = GetPositionAtOffset(retarded_dt);

                const auto retarded_p = GetElectricDipoleAtOffset(retarded_dt);
                const auto retarded_p_dot = ComputeElectricDipoleFirstDerivative(retarded_dt);
                const auto retarded_p_ddot = ComputeElectricDipoleSecondDerivative(retarded_dt);

                const auto R = point - retarded_pos;
                const auto R_mag = R.Length();
                if (R_mag < Quantity<Meter>(Config::epsilon)) {
                    throw std::runtime_error("Field diverges at dipole location");
                }
                const auto hat_R = R / R_mag;

                const auto term1 = Consts::k_E * (3 * retarded_p.Dot(hat_R) * hat_R - retarded_p) / Utils::cube(R_mag);
                const auto term2 = Consts::k_E * (3 * retarded_p_dot.Dot(hat_R) * hat_R - retarded_p_dot) /
                                   (Consts::c * R_mag * R_mag);
                const auto term3 = Consts::k_E * retarded_p_ddot.Cross(hat_R).Cross(hat_R) /
                                   (Consts::c * Consts::c * R_mag);

                return term1 + term2 + term3;
            } catch (const std::runtime_error&) {
                const auto p = GetRelativisticDipoleMomentCurrent();
                const auto hat_r = r.Normalized();
                return Consts::k_E * (3 * p.Dot(hat_r) * hat_r - p) / Utils::cube(dist);
            }
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            if ((point - position_).Length() < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges at dipole location");
            }

            if (fixed_ || position_history_.GetSize() < 4 || velocity_.Length() < Quantity<MeterPerSecond>(Config::epsilon)) {
                return 0_V;
            }

            try {
                const Quantity<Second> retarded_dt = SolveRetardedTimeOffset(point);
                const auto retarded_pos = GetPositionAtOffset(retarded_dt);
                const auto retarded_v = GetVelocityAtOffset(retarded_dt);
                const auto retarded_p = GetElectricDipoleAtOffset(retarded_dt);

                const auto R_vec = point - retarded_pos;
                const auto R = R_vec.Length();
                const auto hat_R = R_vec.Normalized();

                const auto beta = retarded_v / Consts::c;
                const auto denominator_factor = Scalar(1) - hat_R.Dot(beta);

                if (denominator_factor < Scalar(Config::epsilon)) {
                    throw std::runtime_error("Line-of-sight aligns with velocity at light speed");
                }

                return Consts::k_E * retarded_p.Dot(hat_R) / (R * R * denominator_factor);
            } catch (const std::runtime_error&) {
                return 0_V;
            }
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return 0_C;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const override {
            return (point - position_).Length();
        }

        [[nodiscard]] Vector3<Quantity<CoulombMeter>>
        GetElectricDipoleMoment() const {
            return GetRelativisticDipoleMomentCurrent();
        }

        void ApplyMagneticGradientForce(const MagneticSource& field_source) {
            const auto force = ComputeMagneticGradientForce(field_source);
            ApplyForce(force);
        }

        void ApplyMagneticTorque(const MagneticSource& field_source) {
            const auto B = field_source.GetMagneticFieldAt(position_);
            const auto torque = moment_.Cross(B).as<NewtonMeter>();
            ApplyTorque(torque);
        }

        void ApplyElectricGradientForce(const ChargeSource& field_source) {
            const auto force = ComputeElectricGradientForce(field_source);
            ApplyForce(force);
        }

        void ApplyElectricTorque(const ChargeSource& field_source) {
            const auto E = field_source.GetElectricFieldAt(position_);
            const auto p = GetRelativisticDipoleMomentCurrent();
            const auto torque = p.Cross(E).as<NewtonMeter>();
            ApplyTorque(torque);
        }

        void ApplyMagneticForce(const MagneticSource& field_source) {
            ApplyMagneticGradientForce(field_source);
            ApplyMagneticTorque(field_source);
        }

        void ApplyElectricForce(const ChargeSource& field_source) {
            ApplyElectricGradientForce(field_source);
            ApplyElectricTorque(field_source);
        }

        void ApplyElectromagneticForce(const MagneticSource& magnetic_source,
                                       const ChargeSource& electric_source) {
            ApplyMagneticForce(magnetic_source);
            ApplyElectricForce(electric_source);
        }

        void ApplyTorque(const Vector3<Quantity<NewtonMeter>>& torque) {
            external_torque_ += torque;
        }

        void Integrate(const Quantity<Second> dt) noexcept override {
            position_history_.SetDeltaTime(dt);
            velocity_history_.SetDeltaTime(dt);
            acceleration_history_.SetDeltaTime(dt);
            moment_history_.SetDeltaTime(dt);
            electric_dipole_history_.SetDeltaTime(dt);

            Particle::Integrate(dt);

            if (moment_of_inertia_ > 0_kgm2) {
                const auto alpha = external_torque_ / moment_of_inertia_;
                angular_velocity_ += alpha * dt;

                if (const auto delta_angle = angular_velocity_ * dt; delta_angle.Length() > 0_rad) {
                    const auto axis = delta_angle.Normalized();
                    const auto angle = delta_angle.Length();

                    const auto delta_q = Quaternion<Scalar>::FromAxisAngle(axis, angle);
                    orientation_ = delta_q * orientation_;
                    orientation_ = orientation_.Normalized();

                    moment_ = orientation_.Rotated(initial_moment_.Normalized()) * initial_moment_.Length();
                }

                external_torque_ = {};
            }

            position_history_.Register(position_);
            velocity_history_.Register(velocity_);
            acceleration_history_.Register(external_force_ / mass_);
            moment_history_.Register(moment_);

            const auto p = ComputeRelativisticDipoleMoment(velocity_, moment_);
            electric_dipole_history_.Register(p);
        }

        void SetMoment(const Vector3<Quantity<AmpereMeterSquared>>& m) {
            initial_moment_ = m;
            moment_ = orientation_.Rotated(initial_moment_.Normalized()) * initial_moment_.Length();
        }

        [[nodiscard]] const Vector3<Quantity<AmpereMeterSquared>>& GetMoment() const {
            return moment_;
        }

        [[nodiscard]] const Quaternion<Scalar>& GetOrientation() const {
            return orientation_;
        }

        void SetOrientation(const Quaternion<Scalar>& q) {
            orientation_ = q.Normalized();
            moment_ = orientation_.Rotated(initial_moment_.Normalized()) * initial_moment_.Length();
        }

        [[nodiscard]] const Vector3<Quantity<RadianPerSecond>>& GetAngularVelocity() const {
            return angular_velocity_;
        }

        void SetAngularVelocity(const Vector3<Quantity<RadianPerSecond>>& omega) {
            angular_velocity_ = omega;
        }

        [[nodiscard]] Quantity<KilogramMeterSquared> GetMomentOfInertia() const {
            return moment_of_inertia_;
        }

        void SetMomentOfInertia(const Quantity<KilogramMeterSquared>& I) {
            moment_of_inertia_ = I;
        }

    private:
        Vector3<Quantity<AmpereMeterSquared>> moment_;
        Vector3<Quantity<AmpereMeterSquared>> initial_moment_;

        Quantity<KilogramMeterSquared> moment_of_inertia_;
        Quaternion<Scalar> orientation_{};
        Vector3<Quantity<RadianPerSecond>> angular_velocity_{};
        Vector3<Quantity<NewtonMeter>> external_torque_{};

        TimeHistory<Vector3<Quantity<Meter>>> position_history_;
        TimeHistory<Vector3<Quantity<MeterPerSecond>>> velocity_history_;
        TimeHistory<Vector3<Quantity<MeterPerSecondSquared>>> acceleration_history_;
        TimeHistory<Vector3<Quantity<AmpereMeterSquared>>> moment_history_;

        TimeHistory<Vector3<Quantity<CoulombMeter>>> electric_dipole_history_;

        [[nodiscard]] static Vector3<Quantity<CoulombMeter>>
        ComputeRelativisticDipoleMoment(
            const Vector3<Quantity<MeterPerSecond>>& v,
            const Vector3<Quantity<AmpereMeterSquared>>& m
        ) {
            const auto beta = v / Consts::c;
            const auto beta2 = beta.LengthSquared();
            if (beta2 > Scalar(1 - Config::epsilon)) {
                return -v.Cross(m) / (Consts::c * Consts::c);
            }
            const auto gamma = Scalar(1) / Utils::sqrt(Scalar(1) - beta2);
            return -gamma * v.Cross(m) / (Consts::c * Consts::c);
        }

        [[nodiscard]] Vector3<Quantity<CoulombMeter>>
        GetRelativisticDipoleMomentCurrent() const {
            return ComputeRelativisticDipoleMoment(velocity_, moment_);
        }

        [[nodiscard]] Vector3<Quantity<Meter>>
        GetPositionAtOffset(const Quantity<Second>& offset) const {
            if (offset == 0_s) return position_;
            return position_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<MeterPerSecond>>
        GetVelocityAtOffset(const Quantity<Second>& offset) const {
            if (offset == 0_s) return velocity_;
            return velocity_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<MeterPerSecondSquared>>
        GetAccelerationAtOffset(const Quantity<Second>& offset) const {
            if (offset == 0_s) return external_force_ / mass_;
            return acceleration_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<AmpereMeterSquared>>
        GetMomentAtOffset(const Quantity<Second>& offset) const {
            if (offset == 0_s) return moment_;
            return moment_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Vector3<Quantity<CoulombMeter>>
        GetElectricDipoleAtOffset(const Quantity<Second>& offset) const {
            if (offset == 0_s) return GetRelativisticDipoleMomentCurrent();
            return electric_dipole_history_.GetValueByOffset(offset);
        }

        [[nodiscard]] Quantity<Second>
        SolveRetardedTimeOffset(const Vector3<Quantity<Meter>>& point) const {
            if (position_history_.IsEmpty()) {
                return 0_s;
            }

            const Quantity<Meter> R0 = (point - position_).Length();
            Quantity<Second> dt_delay = R0 / Consts::c;

            const Quantity<Second> max_dt = static_cast<double>(position_history_.GetSize())
                                            * position_history_.GetDeltaTime();
            dt_delay = std::clamp(dt_delay, Quantity<Second>(0.0), max_dt);

            for (int iter = 0; iter < Config::max_iterations; ++iter) {
                const auto pos_r = GetPositionAtOffset(dt_delay);
                const auto v_r = GetVelocityAtOffset(dt_delay);

                const Vector3<Quantity<Meter>> R_vec = point - pos_r;
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
        SolveRetardedTimeOffsetBisection(const Vector3<Quantity<Meter>>& point,
                                         Quantity<Second> dt_low,
                                         Quantity<Second> dt_high) const {
            auto ComputeResidual = [&](const Quantity<Second> dt) -> Quantity<Meter> {
                const auto pos_r = GetPositionAtOffset(dt);
                const Quantity<Meter> R = (point - pos_r).Length();
                return R - Consts::c * dt;
            };

            Quantity<Meter> f_low = ComputeResidual(dt_low);
            Quantity<Meter> f_high = ComputeResidual(dt_high);

            if (f_low * f_high > 0_m2) {
                return Utils::abs(f_low) < Utils::abs(f_high) ? dt_low : dt_high;
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

            return (dt_low + dt_high) / 2;
        }

        [[nodiscard]] Vector3<decltype(0_Am2 / 1_s)>
        ComputeMomentFirstDerivative(const Quantity<Second>& retarded_dt) const {
            const auto dt = moment_history_.GetDeltaTime();
            const auto size = moment_history_.GetSize();
            const auto idx = (retarded_dt / dt).to<size_t>();

            if (size < 2) return {};

            if (idx < 1) {
                const auto m0 = GetMomentAtOffset(retarded_dt);
                const auto m1 = GetMomentAtOffset(retarded_dt + dt);
                return (m1 - m0) / dt;
            }
            if (idx >= size - 1) {
                const auto m0 = GetMomentAtOffset(retarded_dt);
                const auto m1 = GetMomentAtOffset(retarded_dt - dt);
                return (m0 - m1) / dt;
            }
            const auto m_next = GetMomentAtOffset(retarded_dt + dt);
            const auto m_prev = GetMomentAtOffset(retarded_dt - dt);
            return (m_next - m_prev) / (2 * dt);
        }

        [[nodiscard]] Vector3<decltype(0_Am2 / 1_s / 1_s)>
        ComputeMomentSecondDerivative(const Quantity<Second>& retarded_dt) const {
            const auto dt = moment_history_.GetDeltaTime();
            const auto size = moment_history_.GetSize();
            const auto idx = (retarded_dt / dt).to<size_t>();

            if (size < 3) return {};

            if (idx < 2) {
                const auto m0 = GetMomentAtOffset(retarded_dt);
                const auto m1 = GetMomentAtOffset(retarded_dt + dt);
                const auto m2 = GetMomentAtOffset(retarded_dt + 2 * dt);
                return (m2 - 2 * m1 + m0) / (dt * dt);
            }
            if (idx >= size - 1) {
                const auto m0 = GetMomentAtOffset(retarded_dt);
                const auto m1 = GetMomentAtOffset(retarded_dt - dt);
                const auto m2 = GetMomentAtOffset(retarded_dt - 2 * dt);
                return (m0 - 2 * m1 + m2) / (dt * dt);
            }
            const auto m_prev = GetMomentAtOffset(retarded_dt + dt);
            const auto m_curr = GetMomentAtOffset(retarded_dt);
            const auto m_next = GetMomentAtOffset(retarded_dt - dt);
            return (m_next - 2 * m_curr + m_prev) / (dt * dt);
        }

        [[nodiscard]] Vector3<decltype(1_Cm / 1_s)>
        ComputeElectricDipoleFirstDerivative(const Quantity<Second>& retarded_dt) const {
            const auto dt = electric_dipole_history_.GetDeltaTime();
            const auto size = electric_dipole_history_.GetSize();
            const auto idx = (retarded_dt / dt).to<size_t>();

            if (size < 2) return {};

            if (idx < 1) {
                const auto p0 = GetElectricDipoleAtOffset(retarded_dt);
                const auto p1 = GetElectricDipoleAtOffset(retarded_dt + dt);
                return (p1 - p0) / dt;
            }
            if (idx >= size - 1) {
                const auto p0 = GetElectricDipoleAtOffset(retarded_dt);
                const auto p1 = GetElectricDipoleAtOffset(retarded_dt - dt);
                return (p0 - p1) / dt;
            }
            const auto p_next = GetElectricDipoleAtOffset(retarded_dt + dt);
            const auto p_prev = GetElectricDipoleAtOffset(retarded_dt - dt);
            return (p_next - p_prev) / (2 * dt);
        }

        [[nodiscard]] Vector3<decltype(1_Cm / 1_s / 1_s)>
        ComputeElectricDipoleSecondDerivative(const Quantity<Second>& retarded_dt) const {
            const auto dt = electric_dipole_history_.GetDeltaTime();
            const auto size = electric_dipole_history_.GetSize();
            const auto idx = (retarded_dt / dt).to<size_t>();

            if (size < 3) return {};

            if (idx < 2) {
                const auto p0 = GetElectricDipoleAtOffset(retarded_dt);
                const auto p1 = GetElectricDipoleAtOffset(retarded_dt + dt);
                const auto p2 = GetElectricDipoleAtOffset(retarded_dt + 2 * dt);
                return (p2 - 2 * p1 + p0) / (dt * dt);
            }
            if (idx >= size - 1) {
                const auto p0 = GetElectricDipoleAtOffset(retarded_dt);
                const auto p1 = GetElectricDipoleAtOffset(retarded_dt - dt);
                const auto p2 = GetElectricDipoleAtOffset(retarded_dt - 2 * dt);
                return (p0 - 2 * p1 + p2) / (dt * dt);
            }
            const auto p_prev = GetElectricDipoleAtOffset(retarded_dt + dt);
            const auto p_curr = GetElectricDipoleAtOffset(retarded_dt);
            const auto p_next = GetElectricDipoleAtOffset(retarded_dt - dt);
            return (p_next - 2 * p_curr + p_prev) / (dt * dt);
        }

        [[nodiscard]] Vector3<decltype(1_mps2 / 1_s)>
        ComputeAccelerationFirstDerivative(const Quantity<Second>& retarded_dt) const {
            const auto dt = acceleration_history_.GetDeltaTime();
            const auto size = acceleration_history_.GetSize();
            const auto idx = (retarded_dt / dt).to<size_t>();

            if (size < 2) return {};

            if (idx < 1) {
                const auto a0 = GetAccelerationAtOffset(retarded_dt);
                const auto a1 = GetAccelerationAtOffset(retarded_dt + dt);
                return (a1 - a0) / dt;
            }
            if (idx >= size - 1) {
                const auto a0 = GetAccelerationAtOffset(retarded_dt);
                const auto a1 = GetAccelerationAtOffset(retarded_dt - dt);
                return (a0 - a1) / dt;
            }
            const auto a_next = GetAccelerationAtOffset(retarded_dt + dt);
            const auto a_prev = GetAccelerationAtOffset(retarded_dt - dt);
            return (a_next - a_prev) / (2 * dt);
        }

        [[nodiscard]] Vector3<Quantity<Newton>>
        ComputeMagneticGradientForce(const MagneticSource& field_source) const {
            const auto eps_m = Quantity<Meter>(Config::epsilon);
            const auto m = moment_;

            Vector3<Quantity<Newton>> force;

            auto pos_x = position_;
            pos_x.x += eps_m;
            auto pos_x_neg = position_;
            pos_x_neg.x -= eps_m;
            const auto dB_dx = (field_source.GetMagneticFieldAt(pos_x) -
                                field_source.GetMagneticFieldAt(pos_x_neg)) / (2 * eps_m);
            force.x = m.Dot(dB_dx);

            auto pos_y = position_;
            pos_y.y += eps_m;
            auto pos_y_neg = position_;
            pos_y_neg.y -= eps_m;
            const auto dB_dy = (field_source.GetMagneticFieldAt(pos_y) -
                                field_source.GetMagneticFieldAt(pos_y_neg)) / (2 * eps_m);
            force.y = m.Dot(dB_dy);

            auto pos_z = position_;
            pos_z.z += eps_m;
            auto pos_z_neg = position_;
            pos_z_neg.z -= eps_m;
            const auto dB_dz = (field_source.GetMagneticFieldAt(pos_z) -
                                field_source.GetMagneticFieldAt(pos_z_neg)) / (2 * eps_m);
            force.z = m.Dot(dB_dz);

            return force;
        }

        [[nodiscard]] Vector3<Quantity<Newton>>
        ComputeElectricGradientForce(const ChargeSource& field_source) const {
            const auto eps_m = Quantity<Meter>(Config::epsilon);
            const auto p = GetRelativisticDipoleMomentCurrent();

            Vector3<Quantity<Newton>> force;

            auto pos_x = position_;
            pos_x.x += eps_m;
            auto pos_x_neg = position_;
            pos_x_neg.x -= eps_m;
            const auto dE_dx = (field_source.GetElectricFieldAt(pos_x) -
                                field_source.GetElectricFieldAt(pos_x_neg)) / (2 * eps_m);
            force.x = p.Dot(dE_dx);

            auto pos_y = position_;
            pos_y.y += eps_m;
            auto pos_y_neg = position_;
            pos_y_neg.y -= eps_m;
            const auto dE_dy = (field_source.GetElectricFieldAt(pos_y) -
                                field_source.GetElectricFieldAt(pos_y_neg)) / (2 * eps_m);
            force.y = p.Dot(dE_dy);

            auto pos_z = position_;
            pos_z.z += eps_m;
            auto pos_z_neg = position_;
            pos_z_neg.z -= eps_m;
            const auto dE_dz = (field_source.GetElectricFieldAt(pos_z) -
                                field_source.GetElectricFieldAt(pos_z_neg)) / (2 * eps_m);
            force.z = p.Dot(dE_dz);

            return force;
        }
    };

}

#endif // PHYTH_MAGNETIC_DIPOLE_HPP