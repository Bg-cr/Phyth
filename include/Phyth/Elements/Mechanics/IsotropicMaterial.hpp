#ifndef PHYTH_ISOTROPIC_MATERIAL_HPP
#define PHYTH_ISOTROPIC_MATERIAL_HPP

#include "Phyth/Core/Quantities.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Matrix3x3.hpp"

namespace Phyth {

    /**
     * @brief Linear elastic isotropic material model
     *
     * IsotropicMaterial implements the constitutive relationship for linear
     * elastic materials with direction-independent mechanical properties.
     * The material is fully defined by any two independent elastic constants,
     * with automatic derivation of all other parameters.
     *
     * The stress-strain relationship follows Hooke's law:
     *   sigma_ij = lambda * delta_ij * epsilon_kk + 2 * mu * epsilon_ij
     *
     * where lambda is Lame's first parameter and mu is the shear modulus.
     *
     * Supported material parameter pairs for construction:
     *   - (Young's modulus E, Poisson's ratio nu)
     *   - (Shear modulus G, Bulk modulus K)
     *
     * All quantities maintain physical dimensions through the Quantity system,
     * ensuring dimensional consistency at compile time.
     *
     * Example:
     *   // Steel-like material
     *   IsotropicMaterial steel(210.0e9_Pa, 0.3);
     *
     *   // Compute stress from strain
     *   Matrix3x3<Scalar> strain = ...;
     *   auto stress = steel.ComputeStress(strain);
     */
    class IsotropicMaterial {
    public:
        /**
         * @brief Default constructor
         *
         * Creates an uninitialized material. All elastic moduli are zero.
         * Use this only when the material will be assigned later.
         */
        IsotropicMaterial() = default;

        /**
         * @brief Construct from Young's modulus and Poisson's ratio
         *
         * This is the most common construction method for engineering materials.
         * All other elastic constants (shear modulus G, bulk modulus K,
         * and Lame's first parameter lambda) are derived automatically.
         *
         * Valid range: -1 < nu < 0.5 for thermodynamic stability
         *   - nu > 0.5: material becomes incompressible or auxetic
         *   - nu < -1: violates stability criteria
         *
         * @param youngs_modulus Young's modulus E (Pa)
         * @param poisson_ratio Poisson's ratio nu (dimensionless)
         */
        IsotropicMaterial(const Quantity<Pascal> youngs_modulus, const Scalar poisson_ratio)
            : youngs_modulus_(youngs_modulus), poisson_ratio_(poisson_ratio) {
            UpdateDerivedParameters();
        }

        /**
         * @brief Construct from shear modulus and bulk modulus
         *
         * This construction method is useful when the material properties
         * are derived from wave propagation measurements or high-pressure
         * experiments.
         *
         * All other constants (Young's modulus E, Poisson's ratio nu,
         * and Lame's first parameter lambda) are derived automatically.
         *
         * @param shear_modulus Shear modulus G (Pa), must be positive
         * @param bulk_modulus Bulk modulus K (Pa), must be positive
         */
        IsotropicMaterial(const Quantity<Pascal> shear_modulus, const Quantity<Pascal> bulk_modulus)
            : shear_modulus_(shear_modulus), bulk_modulus_(bulk_modulus) {
            UpdateDerivedParametersFromGK();
        }

        /**
         * @brief Get Young's modulus
         * @return Young's modulus E in pascals
         */
        [[nodiscard]] Quantity<Pascal> GetYoungsModulus() const { return youngs_modulus_; }

        /**
         * @brief Get Poisson's ratio
         * @return Poisson's ratio nu (dimensionless)
         */
        [[nodiscard]] Scalar GetPoissonRatio() const { return poisson_ratio_; }

        /**
         * @brief Get shear modulus
         * @return Shear modulus G in pascals
         */
        [[nodiscard]] Quantity<Pascal> GetShearModulus() const { return shear_modulus_; }

        /**
         * @brief Get bulk modulus
         * @return Bulk modulus K in pascals
         */
        [[nodiscard]] Quantity<Pascal> GetBulkModulus() const { return bulk_modulus_; }

        /**
         * @brief Get Lame's first parameter
         * @return Lame's first parameter lambda in pascals
         */
        [[nodiscard]] Quantity<Pascal> GetLameLambda() const { return lame_lambda_; }

        /**
         * @brief Compute stress tensor from strain tensor
         *
         * Implements the generalized Hooke's law for isotropic materials:
         *   sigma_ij = lambda * trace(epsilon) * delta_ij + 2 * mu * epsilon_ij
         *
         * The strain tensor is expected to represent infinitesimal strain
         * (small deformation theory). For finite strain, use appropriate
         * hyperelastic models instead.
         *
         * @param strain Strain tensor (dimensionless, symmetric)
         * @return Stress tensor in pascals (symmetric)
         *
         * Example:
         *   Matrix3x3<Scalar> strain;
         *   strain(0, 0) = 0.001;  // 0.1% axial strain
         *   auto stress = material.ComputeStress(strain);
         */
        [[nodiscard]] Matrix3x3<Quantity<Pascal> > ComputeStress(
            const Matrix3x3<Scalar> &strain
        ) const {
            const auto trace_strain = strain.Trace();

            Matrix3x3<Quantity<Pascal> > stress;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    stress(i, j) = lame_lambda_ * trace_strain * (i == j ? 1.0 : 0.0)
                                   + 2.0 * shear_modulus_ * strain(i, j);
                }
            }
            return stress;
        }

        /**
         * @brief Compute strain tensor from stress tensor
         *
         * Implements the inverse form of Hooke's law:
         *   epsilon_ij = (1 + nu) / E * sigma_ij - nu / E * trace(sigma) * delta_ij
         *
         * This is the compliance form of the constitutive relationship,
         * useful for displacement-based formulations.
         *
         * @param stress Stress tensor in pascals (symmetric)
         * @return Strain tensor (dimensionless, symmetric)
         *
         * Example:
         *   Matrix3x3<Quantity<Pascal>> stress;
         *   stress(0, 0) = 100.0e6_Pa;  // uniaxial tension
         *   auto strain = material.ComputeStrain(stress);
         */
        [[nodiscard]] Matrix3x3<Scalar> ComputeStrain(
            const Matrix3x3<Quantity<Pascal> > &stress
        ) const {
            const auto trace_stress = stress.Trace();

            Matrix3x3<Scalar> strain;
            const auto factor1 = (1.0 + poisson_ratio_) / youngs_modulus_;
            const auto factor2 = poisson_ratio_ / youngs_modulus_;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    strain(i, j) = factor1 * stress(i, j)
                                   - factor2 * trace_stress * (i == j ? 1.0 : 0.0);
                }
            }
            return strain;
        }

    private:
        Quantity<Pascal> youngs_modulus_;
        Scalar poisson_ratio_;

        Quantity<Pascal> shear_modulus_;
        Quantity<Pascal> bulk_modulus_;
        Quantity<Pascal> lame_lambda_;

        /**
         * @brief Update derived parameters from (E, nu)
         *
         * Relations:
         *   G = E / (2 * (1 + nu))
         *   K = E / (3 * (1 - 2 * nu))
         *   lambda = E * nu / ((1 + nu) * (1 - 2 * nu))
         *
         * Note: These formulas are singular for nu = 0.5 (incompressible limit)
         * and nu = -1 (stability limit).
         */
        void UpdateDerivedParameters() {
            shear_modulus_ = youngs_modulus_ / (2 * (1 + poisson_ratio_));
            bulk_modulus_ = youngs_modulus_ / (3 * (1 - 2 * poisson_ratio_));
            lame_lambda_ = youngs_modulus_ * poisson_ratio_ / ((1 + poisson_ratio_) * (1 - 2 * poisson_ratio_));
        }

        /**
         * @brief Update derived parameters from (G, K)
         *
         * Relations:
         *   E = 9 * G * K / (3 * K + G)
         *   nu = (3 * K - 2 * G) / (2 * (3 * K + G))
         *   lambda = K - 2/3 * G
         *
         * This formulation is more stable near the incompressible limit.
         */
        void UpdateDerivedParametersFromGK() {
            youngs_modulus_ = 9 * shear_modulus_ * bulk_modulus_ / (3 * bulk_modulus_ + shear_modulus_);
            poisson_ratio_ = (3 * bulk_modulus_ - 2 * shear_modulus_) / (2 * (3 * bulk_modulus_ + shear_modulus_));
            lame_lambda_ = bulk_modulus_ - 2.0 / 3.0 * shear_modulus_;
        }
    };
}

#endif //PHYTH_ISOTROPIC_MATERIAL_HPP