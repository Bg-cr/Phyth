#ifndef PHYTH_LINE_CHARGE_HPP
#define PHYTH_LINE_CHARGE_HPP

#include "ChargeSource.hpp"

namespace Phyth::Electromagnetics {

    /**
     * @brief Uniform line charge along a straight segment
     *
     * Field is computed by numerical integration (default 64 segments).
     */
    class LineCharge : public ChargeSource {
    public:
        LineCharge(const Vector3<Quantity<Meter>>& start,
                   const Vector3<Quantity<Meter>>& end,
                   const Quantity<CoulombPerMeter> lambda)
            : start_(start), end_(end), lambda_(lambda) {
            direction_ = (end_ - start_).Normalized();
            length_ = (end_ - start_).Length();
        }

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const override {
            return NumericalIntegrateField(point, 64);
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const override {
            return AnalyticPotential(point);
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return lambda_ * length_;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const override {
            const auto v = point - start_;
            const auto t = Utils::clamp(v.Dot(direction_), 0_m, length_);
            const auto closest_point = start_ + t * direction_;
            return (point - closest_point).Length();
        }

    private:
        Vector3<Quantity<Meter>> start_, end_;
        Vector3<Scalar> direction_;
        Quantity<Meter> length_;
        Quantity<CoulombPerMeter> lambda_;

        /**
         * @brief Analytic potential for finite line charge
         *
         * @return k_E * lambda * ln((r1 + r2 + L) / (r1 + r2 - L))
         */
        [[nodiscard]] Quantity<Volt>
        AnalyticPotential(const Vector3<Quantity<Meter>>& point) const {
            const auto r1 = (point - start_).Length();
            const auto r2 = (point - end_).Length();

            if (r1 < Quantity<Meter>(Config::epsilon) || r2 < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("Potential diverges on line charge");
            }

            const auto numerator = r1 + r2 + length_;
            const auto denominator = r1 + r2 - length_;

            if (denominator < Quantity<Meter>(Config::epsilon)) {
                throw std::invalid_argument("The point is lied on the line charge");
            }

            return Consts::k_E * lambda_ * Utils::log(numerator / denominator);
        }

        /**
         * @brief Numerical integration of electric field along the line
         * @param point Evaluation point
         * @param segments Number of line segments
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb>>
        NumericalIntegrateField(const Vector3<Quantity<Meter>>& point, const int segments) const {
            Vector3 E{0_NpC, 0_NpC, 0_NpC};
            const auto dl = length_ / segments;

            for (int i = 0; i < segments; ++i) {
                const double t = (i + 0.5) / segments;
                const auto pos = start_ + t * (end_ - start_);
                const auto r = point - pos;
                const auto dist = r.Length();

                if (dist < Quantity<Meter>(Config::epsilon)) continue;

                const auto dq = lambda_ * dl;
                const auto dE = Consts::k_E * dq / (dist * dist) * r.Normalized();
                E += dE;
            }
            return E;
        }
    };

}

#endif //PHYTH_LINE_CHARGE_HPP