#ifndef PHYTH_VOLUME_FORCE_HPP
#define PHYTH_VOLUME_FORCE_HPP

#include "ForceSource.hpp"
#include "Phyth/Elements/Mechanics/IsotropicMaterial.hpp"
#include "Phyth/Physical/Constants.hpp"

namespace Phyth::Mechanics {
    /**
     * @brief Uniform body force distributed over a rectangular box.
     *
     * VolumeForce represents a body force density acting throughout a
     * rectangular volume. The stress field is computed by numerical
     * integration using the midpoint rule with configurable resolution.
     *
     * Example:
     * @code
     *   auto center = Vector3<Quantity<Meter>>(0_m, 0_m, 0_m);
     *   auto force_density = Vector3<Quantity<NewtonPerMeterCubed>>(0, 0, -2700*9.81);
     *   VolumeForce volume(center, 1.0_m, 1.0_m, 1.0_m, force_density, material);
     *   auto stress = volume.GetStressAt(point);
     * @endcode
     */
    class VolumeForce : public ForceSource {
    public:
        /**
         * @brief Construct a rectangular volume force.
         *
         * @param center Center of the box
         * @param width Extent along x
         * @param height Extent along y
         * @param depth Extent along z
         * @param force_density Uniform body force density
         * @param material Material properties
         */
        explicit VolumeForce(
            const Vector3<Quantity<Meter> > &center,
            const Quantity<Meter> width,
            const Quantity<Meter> height,
            const Quantity<Meter> depth,
            const Vector3<Quantity<NewtonPerMeterCubed> > &force_density,
            const IsotropicMaterial &material
        ) noexcept
            : center_(center),
              width_(width),
              height_(height),
              depth_(depth),
              force_density_(force_density),
              material_(material) {
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
            return NumericalIntegrateStress(point, 16, 16, 16);
        }

        /**
         * @brief Get the force vector at a point.
         *
         * Volume force is distributed, so returns zero at any single point.
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetForceAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            return {};
        }

        /**
         * @brief Get the total force of this source.
         *
         * @return force_density * volume (Newtons)
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetTotalForce() const override {
            return force_density_ * width_ * height_ * depth_;
        }

        /**
         * @brief Get the body force density at a point.
         *
         * @param point Position to evaluate
         * @return Body force density at that point
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerMeterCubed> > GetBodyForceDensityAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            const auto half_w = width_ / 2.0;
            const auto half_h = height_ / 2.0;
            const auto half_d = depth_ / 2.0;

            if (const auto local = point - center_; Utils::abs(local.x) <= half_w &&
                                                    Utils::abs(local.y) <= half_h &&
                                                    Utils::abs(local.z) <= half_d) {
                return force_density_;
            }
            return {};
        }

        /**
         * @brief Get the minimum distance from a point to the rectangular volume.
         *
         * @param point Position to compute distance from
         * @return Distance to the closest point in the volume
         */
        [[nodiscard]] Quantity<Meter> GetMinimumDistanceTo(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            const auto half_w = width_ / 2.0;
            const auto half_h = height_ / 2.0;
            const auto half_d = depth_ / 2.0;

            const auto local = point - center_;
            const auto dx = std::clamp(local.x, -half_w, half_w);
            const auto dy = std::clamp(local.y, -half_h, half_h);
            const auto dz = std::clamp(local.z, -half_d, half_d);

            const auto closest = center_ + Vector3(dx, dy, dz);
            return (point - closest).Length();
        }

        /**
         * @brief Get the force density vector.
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerMeterCubed> > GetForceDensity() const {
            return force_density_;
        }

        /**
         * @brief Set the force density vector.
         */
        void SetForceDensity(const Vector3<Quantity<NewtonPerMeterCubed> > &force_density) {
            force_density_ = force_density;
        }

    private:
        Vector3<Quantity<Meter> > center_;
        Quantity<Meter> width_, height_, depth_;
        Vector3<Quantity<NewtonPerMeterCubed> > force_density_;
        IsotropicMaterial material_;

        /**
         * @brief Numerical integration of stress over the volume.
         *
         * Uses the midpoint rule with nx * ny * nz segments.
         * Each voxel contributes dsigma = Kelvin(dF) where dF = force_density * dV.
         */
        [[nodiscard]] Matrix3x3<Quantity<Pascal> > NumericalIntegrateStress(
            const Vector3<Quantity<Meter> > &point,
            const int nx,
            const int ny,
            const int nz
        ) const {
            Matrix3x3<Quantity<Pascal> > stress;

            const auto poisson_ratio = material_.GetPoissonRatio();
            const auto prefactor = -1.0 / (8.0 * Consts::pi * (1.0 - poisson_ratio));

            const auto dx = width_ / nx;
            const auto dy = height_ / ny;
            const auto dz = depth_ / nz;
            const auto dV = dx * dy * dz;
            const auto dF = force_density_ * dV;

            const auto half_w = width_ / 2.0;
            const auto half_h = height_ / 2.0;
            const auto half_d = depth_ / 2.0;

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        const auto x = (i + 0.5) / nx * width_ - half_w;
                        const auto y = (j + 0.5) / ny * height_ - half_h;
                        const auto z = (k + 0.5) / nz * depth_ - half_d;

                        const auto pos = center_ + Vector3(x, y, z);
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
            }

            return stress;
        }
    };
}

#endif // PHYTH_VOLUME_FORCE_HPP
