#ifndef PHYTH_CURRENT_SEGMENT_HPP
#define PHYTH_CURRENT_SEGMENT_HPP

#include "MagneticSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"
#include <cmath>
#include <stdexcept>

namespace Phyth::Electromagnetics {
    /**
     * @brief Magnetic field from a straight current-carrying wire segment
     */
    class CurrentSegment : public MagneticSource {
    public:
        /**
         * @throw std::invalid_argument if segment length is zero
         */
        CurrentSegment(const Vector3<Quantity<Meter> > &start,
                       const Vector3<Quantity<Meter> > &end,
                       const Quantity<Ampere> &current)
            : start_(start), end_(end), current_(current) {
            const auto delta = end_ - start_;
            length_ = delta.Length();

            if (length_ < Quantity<Meter>(Config::epsilon)) {
                throw std::invalid_argument("Zero-length segment");
            }

            direction_ = delta / length_;
        }

        [[nodiscard]] Vector3<Quantity<Tesla> >
        GetMagneticFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            if (IsPointOnWire(point)) {
                throw std::runtime_error("Point lies on wire");
            }

            return AnalyticSolution(point);
        }

        [[nodiscard]] Quantity<Ampere> GetCurrent() const {
            return current_;
        }

        [[nodiscard]] Quantity<Meter> GetLength() const {
            return length_;
        }

        [[nodiscard]] const Vector3<Quantity<Meter> > &GetStart() const {
            return start_;
        }

        [[nodiscard]] const Vector3<Quantity<Meter> > &GetEnd() const {
            return end_;
        }

    private:
        Vector3<Quantity<Meter> > start_;
        Vector3<Quantity<Meter> > end_;
        Vector3<Scalar> direction_;
        Quantity<Meter> length_;
        Quantity<Ampere> current_;

        /**
         * @brief Analytic closed-form solution for a finite straight wire
         *
         * @throw std::runtime_error when rho is zero (on wire axis)
         *
         * @return  (mu_0 * I / (4 * pi * rho)) *
         *              (z2/sqrt(rho^2 + z2^2) - z1/sqrt(rho^2 + z1^2)) * phi_hat
         */
        [[nodiscard]] Vector3<Quantity<Tesla> >
        AnalyticSolution(const Vector3<Quantity<Meter> > &point) const {
            const auto r_vec = point - start_;
            const auto proj = r_vec.Dot(direction_);

            const auto perp_vec = r_vec - proj * direction_;
            const auto rho = perp_vec.Length();

            if (rho < Config::epsilon * length_) {
                return Vector3{0_T, 0_T, 0_T};
            }

            const auto z1 = -proj;
            const auto z2 = length_ - proj;

            const auto rho2 = rho * rho;
            const auto term1 = z2 / Utils::sqrt(rho2 + z2 * z2);
            const auto term2 = z1 / Utils::sqrt(rho2 + z1 * z1);
            const auto angular_factor = term1 - term2;

            const auto B_mag = Consts::mu_0 * current_ / (4 * Consts::pi * rho) * angular_factor;

            const auto perp_hat = perp_vec.Normalized();
            const auto phi_hat = direction_.Cross(perp_hat);

            return B_mag * phi_hat;
        }

        /**
         * @brief Check if a point lies on the wire segment
         *
         * @return true if the point lies on the segment
         */
        [[nodiscard]] bool
        IsPointOnWire(const Vector3<Quantity<Meter> > &point) const {
            const auto r_vec = point - start_;
            const auto proj = r_vec.Dot(direction_);

            const auto tol = Quantity(
                Utils::max(Config::epsilon * length_, Quantity<Meter>(Config::epsilon))
            );

            if (proj < -tol || proj > length_ + tol) {
                return false;
            }

            const auto perp_vec = r_vec - proj * direction_;
            const auto rho = perp_vec.Length();

            return rho < tol;
        }

        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const {
            const auto v = point - start_;
            const auto proj = v.Dot(direction_);

            if (proj <= Quantity<Meter>(0)) {
                return v.Length();
            }

            if (proj >= length_) {
                return (point - end_).Length();
            }

            const auto closest = start_ + proj * direction_;
            return (point - closest).Length();
        }
    };
} // namespace Phyth::Electromagnetics

#endif // PHYTH_CURRENT_SEGMENT_HPP
