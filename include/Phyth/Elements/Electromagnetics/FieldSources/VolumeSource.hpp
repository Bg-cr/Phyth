#ifndef PHYTH_VOLUME_SOURCE_HPP
#define PHYTH_VOLUME_SOURCE_HPP

#include <cassert>

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Uniform volume charge distribution in a rectangular box
     *
     * VolumeCharge represents a rectangular box with uniform volume charge
     * density rho (C/m^3). The box is defined by its center, three orthogonal
     * axes, and dimensions along each axis.
     *
     * Field and potential are computed by numerical integration over the
     * volume using the midpoint rule. Default resolution is 16x16x16 segments.
     *
     * Example:
     *   auto center = Vector3<Quantity<Meter>>(0_m, 0_m, 0_m);
     *   auto ex = Vector3<Scalar>(1, 0, 0);
     *   auto ey = Vector3<Scalar>(0, 1, 0);
     *   auto ez = Vector3<Scalar>(0, 0, 1);
     *   auto rho = 1.0e-6_C/m3;
     *
     *   VolumeCharge volume(center, ex, ey, ez, 1.0_m, 1.0_m, 1.0_m, rho);
     *   auto point = Vector3<Quantity<Meter>>(0_m, 0_m, 0.5_m);
     *   auto E = volume.GetElectricFieldAt(point);
     *   auto V = volume.GetElectricPotentialAt(point);
     *
     * Note: The field is evaluated numerically, so accuracy depends on
     *       the integration resolution. Higher resolution improves accuracy
     *       but increases computation time.
     */
    class VolumeCharge : public ChargeSource {
    public:
        /**
         * @brief Construct a rectangular volume charge
         *
         * @param center Center of the box (meters)
         * @param axis_x Direction of the x-axis (must be unit vector)
         * @param axis_y Direction of the y-axis (must be unit vector)
         * @param axis_z Direction of the z-axis (must be unit vector)
         * @param size_x Length along axis_x (meters)
         * @param size_y Length along axis_y (meters)
         * @param size_z Length along axis_z (meters)
         * @param rho Volume charge density (C/m^3)
         *
         * The three axes must be mutually orthogonal. This is checked at
         * runtime with asserts.
         */
        VolumeCharge(const Vector3<Quantity<Meter> > &center,
                     const Vector3<Scalar> &axis_x,
                     const Vector3<Scalar> &axis_y,
                     const Vector3<Scalar> &axis_z,
                     const Quantity<Meter> size_x,
                     const Quantity<Meter> size_y,
                     const Quantity<Meter> size_z,
                     const Quantity<CoulombPerCubicMeter> rho)
            : center_(center), axis_x_(axis_x.Normalized()),
              axis_y_(axis_y.Normalized()), axis_z_(axis_z.Normalized()),
              size_x_(size_x), size_y_(size_y), size_z_(size_z), rho_(rho) {
            assert(Utils::abs(axis_x_.Dot(axis_y_)) < Config::epsilon);
            assert(Utils::abs(axis_y_.Dot(axis_z_)) < Config::epsilon);
            assert(Utils::abs(axis_z_.Dot(axis_x_)) < Config::epsilon);
        }

        /**
         * @brief Compute the electric field at a point using numerical integration
         *
         * @param point Position to evaluate the field (meters)
         * @return Electric field intensity E (N/C)
         *
         * The field is integrated using the midpoint rule with 16x16x16 segments.
         * Points inside the volume are skipped (their contribution is ignored).
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegrateField(point, 16, 16, 16);
        }

        /**
         * @brief Compute the electric potential at a point using numerical integration
         *
         * @param point Position to evaluate the potential (meters)
         * @return Electric potential V (Volts)
         *
         * The potential is integrated using the midpoint rule with 16x16x16 segments.
         * Points inside the volume are skipped (their contribution is ignored).
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegratePotential(point, 16, 16, 16);
        }

        /**
         * @brief Get the total charge of the volume
         *
         * @return rho * size_x * size_y * size_z (Coulombs)
         */
        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return rho_ * size_x_ * size_y_ * size_z_;
        }

        /**
         * @brief Get the minimum distance from a point to the rectangular box
         *
         * @param point Position to compute distance from
         * @return Distance to the closest point on the box surface or interior (meters)
         *
         * The closest point is found by projecting the point onto the box axes
         * and clamping to the box bounds. If the point is inside the box,
         * the minimum distance is 0.
         */
        [[nodiscard]] Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter> > &point) const override {
            const auto local = point - center_;
            auto px = local.Dot(axis_x_);
            auto py = local.Dot(axis_y_);
            auto pz = local.Dot(axis_z_);

            px = Utils::clamp(px, -size_x_ / 2, size_x_ / 2);
            py = Utils::clamp(py, -size_y_ / 2, size_y_ / 2);
            pz = Utils::clamp(pz, -size_z_ / 2, size_z_ / 2);

            const auto closest = center_ + px * axis_x_ + py * axis_y_ + pz * axis_z_;
            return (point - closest).Length();
        }

        /**
         * @brief Get the charge density at a point
         *
         * @param point Position to evaluate (meters)
         * @return rho if point is inside the box, 0 otherwise
         *
         * This overrides the default implementation to return the actual
         * charge density for points inside the volume.
         */
        [[nodiscard]] Quantity<CoulombPerCubicMeter>
        GetChargeDensityAt(const Vector3<Quantity<Meter> > &point) const override {
            const auto local = point - center_;
            const auto px = Utils::abs(local.Dot(axis_x_));
            const auto py = Utils::abs(local.Dot(axis_y_));
            const auto pz = Utils::abs(local.Dot(axis_z_));

            if (px <= size_x_ / 2 && py <= size_y_ / 2 && pz <= size_z_ / 2) {
                return rho_;
            }
            return 0_Cpm3;
        }

    private:
        Vector3<Quantity<Meter> > center_;
        Vector3<Scalar> axis_x_, axis_y_, axis_z_;
        Quantity<Meter> size_x_, size_y_, size_z_;
        Quantity<CoulombPerCubicMeter> rho_;

        /**
         * @brief Numerical integration of electric field over the volume
         *
         * Uses the midpoint rule with nx * ny * nz segments.
         * Each segment contributes dE = k_E * dq / r^2 * hat_r.
         *
         * @param point Evaluation point
         * @param nx Number of segments along x axis
         * @param ny Number of segments along y axis
         * @param nz Number of segments along z axis
         * @return Electric field intensity E (N/C)
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        NumericalIntegrateField(const Vector3<Quantity<Meter> > &point,
                                const int nx, const int ny, const int nz) const {
            Vector3 E{0_NpC, 0_NpC, 0_NpC};
            const auto dx = size_x_ / nx;
            const auto dy = size_y_ / ny;
            const auto dz = size_z_ / nz;
            const auto dV = dx * dy * dz;

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        const auto x = (i + 0.5) / nx * size_x_ - size_x_ / 2;
                        const auto y = (j + 0.5) / ny * size_y_ - size_y_ / 2;
                        const auto z = (k + 0.5) / nz * size_z_ - size_z_ / 2;
                        const auto pos = center_ + x * axis_x_ + y * axis_y_ + z * axis_z_;

                        const auto r = point - pos;
                        const auto dist = r.Length();
                        if (dist < Quantity<Meter>(Config::epsilon)) continue;

                        const auto dq = rho_ * dV;
                        const auto dE = Consts::k_E * dq / (dist * dist) * r.Normalized();
                        E += dE;
                    }
                }
            }
            return E;
        }

        /**
         * @brief Numerical integration of electric potential over the volume
         *
         * Uses the midpoint rule with nx * ny * nz segments.
         * Each segment contributes dV = k_E * dq / r.
         *
         * @param point Evaluation point
         * @param nx Number of segments along x axis
         * @param ny Number of segments along y axis
         * @param nz Number of segments along z axis
         * @return Electric potential V (Volts)
         */
        [[nodiscard]] Quantity<Volt>
        NumericalIntegratePotential(const Vector3<Quantity<Meter> > &point,
                                    const int nx, const int ny, const int nz) const {
            auto V = 0_V;
            const auto dx = size_x_ / nx;
            const auto dy = size_y_ / ny;
            const auto dz = size_z_ / nz;
            const auto dV = dx * dy * dz;

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        const auto x = (i + 0.5) / nx * size_x_ - size_x_ / 2;
                        const auto y = (j + 0.5) / ny * size_y_ - size_y_ / 2;
                        const auto z = (k + 0.5) / nz * size_z_ - size_z_ / 2;
                        const auto pos = center_ + x * axis_x_ + y * axis_y_ + z * axis_z_;

                        const auto r = point - pos;
                        const auto dist = r.Length();
                        if (dist < Quantity<Meter>(Config::epsilon)) continue;

                        V += Consts::k_E * rho_ * dV / dist;
                    }
                }
            }
            return V;
        }
    };
}

#endif // PHYTH_VOLUME_SOURCE_HPP
