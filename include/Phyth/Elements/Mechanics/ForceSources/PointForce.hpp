#ifndef PHYTH_POINT_FORCE_HPP
#define PHYTH_POINT_FORCE_HPP

#include "ForceSource.hpp"
#include "Phyth/Elements/Mechanics/IsotropicMaterial.hpp"
#include "Phyth/Physical/Constants.hpp"

namespace Phyth::Mechanics {
    /**
     * @brief Concentrated force acting at a single point in an elastic medium.
     *
     * Produces a stress field via Kelvin's solution (Green's function for an
     * infinite elastic medium).
     *
     * Formula:
     *   sigma_ij(vec_r) = - (1 / (8 * pi * (1 - nu))) * [
     *       (1 - 2*nu) * (vec_F_i * vec_r_j + vec_F_j * vec_r_i - delta_ij * (vec_F dot vec_r)) / r^3
     *       - 3 * (vec_F dot vec_r) * vec_r_i * vec_r_j / r^5
     *   ]
     *
     * Example:
     * @code
     *   auto force = std::make_shared<PointForce>(
     *       origin,
     *       force,
     *       material
     *   );
     *   elastic_field.AddSource(force);
     * @endcode
     */
    class PointForce : public ForceSource {
    public:
        /**
         * @brief Construct a point force.
         *
         * @param position Location of the point force
         * @param force Force vector
         * @param material Material properties (uses Poisson's ratio only)
         */
        explicit PointForce(const Vector3<Quantity<Meter> > &position,
                            const Vector3<Quantity<Newton> > &force,
                            const IsotropicMaterial &material) noexcept
            : position_(position), force_(force), material_(material) {
        }

        /**
         * @brief Compute stress tensor at a point.
         *
         * @param point Position to evaluate
         * @return Stress tensor in pascals
         * @throws std::runtime_error If point is at the force location
         *
         * @note Sign convention: tension positive (positive stress = tensile).
         * @note If you prefer "compression positive", multiply the result by -1.
         */
        [[nodiscard]] Matrix3x3<Quantity<Pascal> > GetStressAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            const auto r = point - position_;
            const auto r_mag = r.Length();

            if (r_mag < Quantity<Meter>(Config::epsilon)) {
                throw std::runtime_error("PointForce: stress diverges at force location");
            }

            const auto poisson_ratio = material_.GetPoissonRatio();
            const Quantity<Joule> F_dot_r = force_.Dot(r);
            const Scalar kelvin_prefactor = -1.0 / (8.0 * Consts::pi * (1.0 - poisson_ratio));

            const auto r3 = r_mag * r_mag * r_mag;
            const auto r5 = r3 * r_mag * r_mag;

            const Scalar one_minus_2nu = 1.0 - 2.0 * poisson_ratio;
            const Scalar diagonal_coeff = 2.0 * one_minus_2nu;
            const Quantity<Joule> dipole_coeff = 3.0 * F_dot_r;

            Matrix3x3<Quantity<Pascal> > stress;

            // Diagonal components
            stress(0, 0) = kelvin_prefactor * (diagonal_coeff * force_.x * r.x / r3 - one_minus_2nu * F_dot_r / r3
                                               -
                                               dipole_coeff * r.x * r.x / r5);
            stress(1, 1) = kelvin_prefactor * (diagonal_coeff * force_.y * r.y / r3 - one_minus_2nu * F_dot_r / r3
                                               -
                                               dipole_coeff * r.y * r.y / r5);
            stress(2, 2) = kelvin_prefactor * (diagonal_coeff * force_.z * r.z / r3 - one_minus_2nu * F_dot_r / r3
                                               -
                                               dipole_coeff * r.z * r.z / r5);

            // Off-diagonal components (symmetric)
            stress(0, 1) = stress(1, 0) = kelvin_prefactor * (
                                              one_minus_2nu * (force_.x * r.y + force_.y * r.x) / r3 - dipole_coeff * r.
                                              x * r.y / r5);
            stress(0, 2) = stress(2, 0) = kelvin_prefactor * (
                                              one_minus_2nu * (force_.x * r.z + force_.z * r.x) / r3 - dipole_coeff * r.
                                              x * r.z / r5);
            stress(1, 2) = stress(2, 1) = kelvin_prefactor * (
                                              one_minus_2nu * (force_.y * r.z + force_.z * r.y) / r3 - dipole_coeff * r.
                                              y * r.z / r5);

            return stress;
        }

        /**
         * @brief Get the force vector at a point.
         *
         * Returns the force only at the exact force location, otherwise zero.
         *
         * @param point Position to evaluate.
         * @return Force vector at that point.
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetForceAt(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            if ((point - position_).Length() < Quantity<Meter>(Config::epsilon)) {
                return force_;
            }
            return {};
        }

        /**
         * @brief Get the total force of this source.
         *
         * @return Total force vector.
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetTotalForce() const override {
            return force_;
        }

        /**
         * @brief Get the minimum distance from a point to the force location.
         *
         * @param point Position to compute distance from.
         * @return Distance to the force location.
         */
        [[nodiscard]] Quantity<Meter> GetMinimumDistanceTo(
            const Vector3<Quantity<Meter> > &point
        ) const override {
            return (point - position_).Length();
        }

        /**
         * @brief Get the force location.
         *
         * @return Position of the point force.
         */
        [[nodiscard]] Vector3<Quantity<Meter> > GetPosition() const {
            return position_;
        }

        /**
         * @brief Get the force vector.
         *
         * @return Force vector.
         */
        [[nodiscard]] Vector3<Quantity<Newton> > GetForce() const {
            return force_;
        }

    private:
        Vector3<Quantity<Meter> > position_;
        Vector3<Quantity<Newton> > force_;
        IsotropicMaterial material_;
    };
}

#endif // PHYTH_POINT_FORCE_HPP
