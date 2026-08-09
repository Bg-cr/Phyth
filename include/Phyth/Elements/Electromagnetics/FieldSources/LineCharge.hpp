#ifndef PHYTH_LINE_CHARGE_HPP
#define PHYTH_LINE_CHARGE_HPP

#include "ChargeSource.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Uniform line charge along a straight segment
     *
     * LineCharge represents a straight segment with uniform linear charge
     * density lambda (C/m). The electric field is computed by numerical
     * integration (64 segments by default), while the potential uses an
     * analytic closed-form solution.
     *
     * Example:
     *   auto start = Vector3<Quantity<Meter>>(-1.0_m, 0_m, 0_m);
     *   auto end = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     *   auto lambda = 1.0e-9_C/m;
     *
     *   LineCharge line(start, end, lambda);
     *   auto point = Vector3<Quantity<Meter>>(0_m, 0.1_m, 0_m);
     *   auto E = line.GetElectricFieldAt(point);
     *   auto V = line.GetElectricPotentialAt(point);
     *
     * Note: The field and potential diverge on the line charge itself.
     *       Querying at points on the line throws an exception.
     */
    class LineCharge : public ChargeSource {
    public:
        /**
         * @brief Construct a uniform line charge
         *
         * @param start Starting point of the line (meters)
         * @param end Ending point of the line (meters)
         * @param lambda Linear charge density (C/m)
         */
        LineCharge(const Vector3<Quantity<Meter> > &start,
                   const Vector3<Quantity<Meter> > &end,
                   const Quantity<CoulombPerMeter> lambda)
            : start_(start), end_(end), lambda_(lambda) {
            direction_ = (end_ - start_).Normalized();
            length_ = (end_ - start_).Length();
        }

        /**
         * @brief Compute the electric field at a point using numerical integration
         *
         * @param point Position to evaluate the field (meters)
         * @return Electric field intensity E (N/C)
         *
         * The field is integrated using the midpoint rule with 64 segments.
         * Points on the line are skipped (their contribution is ignored).
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegrateField(point, 64);
        }

        /**
         * @brief Compute the electric potential at a point using analytic formula
         *
         * The potential for a finite line charge is:
         *   V = k_E * lambda * ln((r1 + r2 + L) / (r1 + r2 - L))
         * where r1, r2 are distances to the endpoints and L is the segment length.
         *
         * @param point Position to evaluate the potential (meters)
         * @return Electric potential V (Volts)
         *
         * @throws std::runtime_error If point is on the line charge
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            return AnalyticPotential(point);
        }

        /**
         * @brief Get the total charge of the line
         *
         * @return lambda * length (Coulombs)
         */
        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return lambda_ * length_;
        }

        /**
         * @brief Get the minimum distance from a point to the line segment
         *
         * @param point Position to compute distance from
         * @return Distance to the closest point on the segment (meters)
         */
        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            const auto v = point - start_;
            const auto t = Utils::clamp(v.Dot(direction_), 0_m, length_);
            const auto closest_point = start_ + t * direction_;
            return (point - closest_point).Length();
        }

    private:
        Vector3<Quantity<Meter> > start_, end_;
        Vector3<Scalar> direction_;
        Quantity<Meter> length_;
        Quantity<CoulombPerMeter> lambda_;

        /**
         * @brief Analytic potential for finite line charge
         *
         * The formula is:
         *   V = k_E * lambda * ln((r1 + r2 + L) / (r1 + r2 - L))
         *
         * where r1 = |point - start|, r2 = |point - end|, L = segment length.
         *
         * @param point Position to evaluate the potential
         * @return Electric potential V (Volts)
         *
         * @throws std::runtime_error If point is on the line charge
         *
         * Note: The denominator (r1 + r2 - L) approaches zero when the point
         *       lies on the line segment, causing a logarithmic divergence.
         */
        [[nodiscard]] Quantity<Volt>
        AnalyticPotential(const Vector3<Quantity<Meter> > &point) const {
            const auto r1 = (point - start_).Length();
            const auto r2 = (point - end_).Length();

            if (r1 < Quantity<Meter>(Config::epsilon) || r2 < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("LineCharge: potential diverges on line charge");
            }

            const auto numerator = r1 + r2 + length_;
            const auto denominator = r1 + r2 - length_;

            if (denominator < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("LineCharge: point lies on the line charge");
            }

            return Consts::k_E * lambda_ * Utils::log(numerator / denominator);
        }

        /**
         * @brief Numerical integration of electric field along the line
         *
         * Uses the midpoint rule with the specified number of segments.
         * Each segment contributes dE = k_E * dq / r^2 * hat_r.
         *
         * @param point Evaluation point
         * @param segments Number of line segments for integration
         * @return Electric field intensity E (N/C)
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        NumericalIntegrateField(const Vector3<Quantity<Meter> > &point, const int segments) const {
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

#endif // PHYTH_LINE_CHARGE_HPP
