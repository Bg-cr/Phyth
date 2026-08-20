#ifndef PHYTH_SURFACE_FORCE_HPP
#define PHYTH_SURFACE_FORCE_HPP

#include "ForceSource.hpp"
#include "Phyth/Elements/Mechanics/IsotropicMaterial.hpp"
#include "Phyth/Physical/Constants.hpp"

#include <cassert>

namespace Phyth::Mechanics {
    /**
     * @brief Uniform surface traction on a rectangular patch.
     *
     * SurfaceForce represents a traction distribution on a flat rectangular
     * surface. The stress field is computed by numerical integration using
     * the midpoint rule with configurable resolution.
     *
     * Example:
     * @code
     *   auto center = Vector3<Quantity<Meter>>(0_m, 0_m, 0_m);
     *   auto u_axis = Vector3<Scalar>(1, 0, 0);
     *   auto v_axis = Vector3<Scalar>(0, 1, 0);
     *   auto traction = Vector3<Quantity<Pascal>>(0_Pa, 0_Pa, -1000_Pa);
     *
     *   SurfaceForce surface(center, u_axis, v_axis, 1.0_m, 1.0_m, traction, material);
     *   auto stress = surface.GetStressAt(point);
     * @endcode
     */
    class SurfaceForce : public ForceSource {
    public:
        /**
         * @brief Construct a rectangular surface traction.
         *
         * @param center Center of the rectangle
         * @param u_axis Direction of the width axis (unit vector)
         * @param v_axis Direction of the height axis (unit vector)
         * @param width Length along u_axis
         * @param height Length along v_axis
         * @param traction Uniform traction vector
         * @param material Material properties
         *
         * u_axis and v_axis must be orthogonal. This is checked at runtime.
         */
        explicit SurfaceForce(
            const Vector3<Quantity<Meter> > &center,
            const Vector3<Scalar> &u_axis,
            const Vector3<Scalar> &v_axis,
            const Quantity<Meter> width,
            const Quantity<Meter> height,
            const Vector3<Quantity<Pascal> > &traction,
            const IsotropicMaterial &material
        ) noexcept
            : center_(center),
              u_axis_(u_axis.Normalized()),
              v_axis_(v_axis.Normalized()),
              width_(width),
              height_(height),
              traction_(traction),
              material_(material) {
            assert(Utils::abs(u_axis_.Dot(v_axis_)) < Config::epsilon);
        }

        /**
         * @brief Compute stress tensor at a point using numerical integration.
         *
         * @param point Position to evaluate
         * @return Stress tensor in pascals
         */
        [[nodiscard]] Matrix3x3<Quantity<Pascal> > GetStressAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            return NumericalIntegrateStress(point, 32, 32);
        }

        /**
         * @brief Get the force vector at a point.
         *
         * Surface force is distributed, so returns zero at any single point.
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetForceAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            return {};
        }

        /**
         * @brief Get the total force of this source.
         *
         * @return traction * area
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetTotalForce() const override {
            return traction_ * width_ * height_;
        }

        /**
         * @brief Get the minimum distance from a point to the rectangular surface.
         *
         * @param point Position to compute distance from
         * @return Distance to the closest point on the surface
         */
        [[nodiscard]] Quantity<Meter> GetMinimumDistanceTo(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            const auto local = point - center_;
            auto u_proj = local.Dot(u_axis_);
            auto v_proj = local.Dot(v_axis_);

            u_proj = std::clamp(u_proj, -width_ / 2, width_ / 2);
            v_proj = std::clamp(v_proj, -height_ / 2, height_ / 2);

            const auto closest = center_ + u_proj * u_axis_ + v_proj * v_axis_;
            return (point - closest).Length();
        }

        /**
         * @brief Get the traction vector.
         */
        [[nodiscard]] Vector3<Quantity<Pascal> > GetTraction() const {
            return traction_;
        }

        /**
         * @brief Set the traction vector.
         */
        void SetTraction(const Vector3<Quantity<Pascal> > &traction) {
            traction_ = traction;
        }

    private:
        Vector3<Quantity<Meter> > center_;
        Vector3<Scalar> u_axis_, v_axis_;
        Quantity<Meter> width_, height_;
        Vector3<Quantity<Pascal> > traction_;
        IsotropicMaterial material_;

        /**
         * @brief Numerical integration of stress over the surface.
         *
         * Uses the midpoint rule with nu * nv segments.
         * Each segment contributes dsigma = Kelvin(dF) where dF = traction * dA.
         *
         * @param point Evaluation point
         * @param nu Number of segments along u axis
         * @param nv Number of segments along v axis
         * @return Stress tensor in pascals
         */
        [[nodiscard]] Matrix3x3<Quantity<Pascal> > NumericalIntegrateStress(
            const Vector3<Quantity<Meter> > &point,
            const int nu,
            const int nv
        ) const {
            Matrix3x3<Quantity<Pascal> > stress;

            const auto poisson_ratio = material_.GetPoissonRatio();
            const auto prefactor = -1.0 / (8.0 * Consts::pi * (1.0 - poisson_ratio));

            const auto du = width_ / nu;
            const auto dv = height_ / nv;
            const auto dA = du * dv;
            const auto dF = traction_ * dA;

            for (int i = 0; i < nu; ++i) {
                for (int j = 0; j < nv; ++j) {
                    const auto u = (i + 0.5) / nu * width_ - width_ / 2.0;
                    const auto v = (j + 0.5) / nv * height_ - height_ / 2.0;

                    const auto pos = center_ + u * u_axis_ + v * v_axis_;
                    const auto R = point - pos;
                    const auto R_mag = R.Length();

                    if (R_mag < Quantity<Meter>(Config::epsilon)) {
                        continue;
                    }

                    const auto R3 = R_mag * R_mag * R_mag;
                    const auto R5 = R3 * R_mag * R_mag;

                    const auto one_minus_2nu = 1.0 - 2.0 * poisson_ratio;
                    const auto dF_dot_R = dF.Dot(R);

                    for (int m = 0; m < 3; ++m) {
                        for (int n = 0; n < 3; ++n) {
                            const auto term1 = one_minus_2nu * (
                                                   dF[m] * R[n] + dF[n] * R[m]
                                                   - (m == n ? 1.0 : 0.0) * dF_dot_R
                                               ) / R3;

                            const auto term2 = 3.0 * dF_dot_R * R[m] * R[n] / R5;

                            stress(m, n) += prefactor * (term1 - term2);
                        }
                    }
                }
            }

            return stress;
        }
    };
}

#endif // PHYTH_SURFACE_FORCE_HPP
