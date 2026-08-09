#ifndef PHYTH_SURFACE_CHARGE_HPP
#define PHYTH_SURFACE_CHARGE_HPP

#include <cassert>

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Uniform surface charge distribution on a rectangular patch
     *
     * SurfaceCharge represents a flat rectangular surface with uniform
     * charge density sigma (C/m^2). The surface is defined by its center,
     * two orthogonal axes (u and v), and dimensions (width and height).
     *
     * Field and potential are computed by numerical integration over the
     * surface using the midpoint rule. Default resolution is 32x32 segments.
     *
     * Example:
     *   auto center = Vector3<Quantity<Meter>>(0_m, 0_m, 0_m);
     *   auto u_axis = Vector3<Scalar>(1, 0, 0);
     *   auto v_axis = Vector3<Scalar>(0, 1, 0);
     *   auto sigma = 1.0e-9_C/m2;
     *
     *   SurfaceCharge surface(center, u_axis, v_axis, 1.0_m, 1.0_m, sigma);
     *   auto point = Vector3<Quantity<Meter>>(0_m, 0_m, 0.1_m);
     *   auto E = surface.GetElectricFieldAt(point);
     *   auto V = surface.GetElectricPotentialAt(point);
     *
     * Note: The field is evaluated numerically, so accuracy depends on
     *       the integration resolution.
     */
    class SurfaceCharge : public ChargeSource {
    public:
        /**
         * @brief Construct a rectangular surface charge
         *
         * @param center Center of the rectangle (meters)
         * @param u_axis Direction of the width axis (must be unit vector)
         * @param v_axis Direction of the height axis (must be unit vector)
         * @param width Length along u_axis (meters)
         * @param height Length along v_axis (meters)
         * @param sigma Surface charge density (C/m^2)
         *
         * u_axis and v_axis must be orthogonal. This is checked at runtime
         * with an assert.
         */
        SurfaceCharge(const Vector3<Quantity<Meter> > &center,
                      const Vector3<Scalar> &u_axis,
                      const Vector3<Scalar> &v_axis,
                      const Quantity<Meter> width,
                      const Quantity<Meter> height,
                      const Quantity<CoulombPerSquareMeter> sigma)
            : center_(center), u_axis_(u_axis.Normalized()),
              v_axis_(v_axis.Normalized()),
              width_(width), height_(height), sigma_(sigma) {
            assert(Utils::abs(u_axis_.Dot(v_axis_)) < Config::epsilon);
        }

        /**
         * @brief Compute the electric field at a point using numerical integration
         *
         * @param point Position to evaluate the field (meters)
         * @return Electric field intensity E (N/C)
         *
         * The field is integrated using the midpoint rule with 32x32 segments.
         * Points on the surface are skipped (their contribution is ignored).
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegrateField(point, 32, 32);
        }

        /**
         * @brief Compute the electric potential at a point using numerical integration
         *
         * @param point Position to evaluate the potential (meters)
         * @return Electric potential V (Volts)
         *
         * The potential is integrated using the midpoint rule with 32x32 segments.
         * Points on the surface are skipped (their contribution is ignored).
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegratePotential(point, 32, 32);
        }

        /**
         * @brief Get the total charge of the surface
         *
         * @return sigma * width * height (Coulombs)
         */
        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return sigma_ * width_ * height_;
        }

        /**
         * @brief Get the minimum distance from a point to the rectangular surface
         *
         * @param point Position to compute distance from
         * @return Distance to the closest point on the surface (meters)
         *
         * The closest point is found by projecting the point onto the surface
         * plane and clamping to the rectangle bounds.
         */
        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            const auto local = point - center_;
            auto u_proj = local.Dot(u_axis_);
            auto v_proj = local.Dot(v_axis_);

            u_proj = std::clamp(u_proj, -width_ / 2, width_ / 2);
            v_proj = std::clamp(v_proj, -height_ / 2, height_ / 2);

            const auto closest = center_ + u_proj * u_axis_ + v_proj * v_axis_;
            return (point - closest).Length();
        }

        /** @return The surface charge density sigma (C/m^2) */
        [[nodiscard]] Quantity<CoulombPerSquareMeter>
        GetSurfaceDensity() const { return sigma_; }

    private:
        Vector3<Quantity<Meter> > center_;
        Vector3<Scalar> u_axis_, v_axis_;
        Quantity<Meter> width_, height_;
        Quantity<CoulombPerSquareMeter> sigma_;

        /**
         * @brief Numerical integration of electric field over the surface
         *
         * Uses the midpoint rule with nu * nv segments.
         * Each segment contributes dE = k_E * dq / r^2 * hat_r.
         *
         * @param point Evaluation point
         * @param nu Number of segments along u axis
         * @param nv Number of segments along v axis
         * @return Electric field intensity E (N/C)
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        NumericalIntegrateField(const Vector3<Quantity<Meter> > &point,
                                const int nu, const int nv) const {
            Vector3 E{0_NpC, 0_NpC, 0_NpC};
            const auto du = width_ / nu;
            const auto dv = height_ / nv;

            for (int i = 0; i < nu; ++i) {
                for (int j = 0; j < nv; ++j) {
                    const auto u = (i + 0.5) / nu * width_ - width_ / 2;
                    const auto v = (j + 0.5) / nv * height_ - height_ / 2;
                    const auto pos = center_ + u * u_axis_ + v * v_axis_;

                    const auto r = point - pos;
                    const auto dist = r.Length();
                    if (dist < Quantity<Meter>(Config::epsilon)) continue;

                    const auto dq = sigma_ * du * dv;
                    const auto dE = Consts::k_E * dq / (dist * dist) * r.Normalized();
                    E += dE;
                }
            }
            return E;
        }

        /**
         * @brief Numerical integration of electric potential over the surface
         *
         * Uses the midpoint rule with nu * nv segments.
         * Each segment contributes dV = k_E * dq / r.
         *
         * @param point Evaluation point
         * @param nu Number of segments along u axis
         * @param nv Number of segments along v axis
         * @return Electric potential V (Volts)
         */
        [[nodiscard]] Quantity<Volt>
        NumericalIntegratePotential(const Vector3<Quantity<Meter> > &point,
                                    const int nu, const int nv) const {
            auto V = 0_V;
            const auto du = width_ / nu;
            const auto dv = height_ / nv;

            for (int i = 0; i < nu; ++i) {
                for (int j = 0; j < nv; ++j) {
                    const auto u = (i + 0.5) / nu * width_ - width_ / 2;
                    const auto v = (j + 0.5) / nv * height_ - height_ / 2;
                    const auto pos = center_ + u * u_axis_ + v * v_axis_;

                    const auto r = point - pos;
                    const auto dist = r.Length();
                    if (dist < Quantity<Meter>(Config::epsilon)) continue;

                    const auto dq = sigma_ * du * dv;
                    V += Consts::k_E * dq / dist;
                }
            }
            return V;
        }
    };
}

#endif // PHYTH_SURFACE_CHARGE_HPP
