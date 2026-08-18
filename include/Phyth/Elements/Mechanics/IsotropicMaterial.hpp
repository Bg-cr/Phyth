#ifndef PHYTH_ISOTROPIC_MATERIAL_HPP
#define PHYTH_ISOTROPIC_MATERIAL_HPP

#include "Phyth/Core/Quantities.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Matrix3x3.hpp"

namespace Phyth {
    class IsotropicMaterial {
    public:
        IsotropicMaterial() = default;

        IsotropicMaterial(const Quantity<Pascal> youngs_modulus, const Scalar poisson_ratio)
            : youngs_modulus_(youngs_modulus), poisson_ratio_(poisson_ratio) {
            UpdateDerivedParameters();
        }

        IsotropicMaterial(const Quantity<Pascal> shear_modulus, const Quantity<Pascal> bulk_modulus)
            : shear_modulus_(shear_modulus), bulk_modulus_(bulk_modulus) {
            UpdateDerivedParametersFromGK();
        }

        [[nodiscard]] Quantity<Pascal> GetYoungsModulus() const { return youngs_modulus_; }
        [[nodiscard]] Scalar GetPoissonRatio() const { return poisson_ratio_; }
        [[nodiscard]] Quantity<Pascal> GetShearModulus() const { return shear_modulus_; }
        [[nodiscard]] Quantity<Pascal> GetBulkModulus() const { return bulk_modulus_; }
        [[nodiscard]] Quantity<Pascal> GetLameLambda() const { return lame_lambda_; }

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

        void UpdateDerivedParameters() {
            shear_modulus_ = youngs_modulus_ / (2 * (1 + poisson_ratio_));
            bulk_modulus_ = youngs_modulus_ / (3 * (1 - 2 * poisson_ratio_));
            lame_lambda_ = youngs_modulus_ * poisson_ratio_ / ((1 + poisson_ratio_) * (1 - 2 * poisson_ratio_));
        }

        void UpdateDerivedParametersFromGK() {
            youngs_modulus_ = 9 * shear_modulus_ * bulk_modulus_ / (3 * bulk_modulus_ + shear_modulus_);
            poisson_ratio_ = (3 * bulk_modulus_ - 2 * shear_modulus_) / (2 * (3 * bulk_modulus_ + shear_modulus_));
            lame_lambda_ = bulk_modulus_ - 2.0 / 3.0 * shear_modulus_;
        }
    };
}

#endif //PHYTH_ISOTROPIC_MATERIAL_HPP
