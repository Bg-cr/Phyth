#ifndef PHYTH_TETRAHEDRON_HPP
#define PHYTH_TETRAHEDRON_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "Phyth/Tools/Matrix.hpp"
#include "Phyth/Tools/Vector3.hpp"
#include "Phyth/Elements/Mechanics/IsotropicMaterial.hpp"

#include <array>

namespace Phyth::Mechanics {
    class Tetrahedron {
    public:
        Tetrahedron(const std::array<Vector3<Quantity<Meter> >, 4> &coords,
                    const IsotropicMaterial &material)
            : coords_(coords), material_(material) {
            ComputeGeometry();
        }

        [[nodiscard]] Matrix<12, 12, Quantity<NewtonPerMeter > > GetStiffness() const {
            return B_.Transpose() * D_ * B_ * volume_;
        }

        [[nodiscard]] Matrix<6, 12, Quantity<UnitRecT<Meter>>> GetBMatrix() const { return B_; }
        [[nodiscard]] Matrix<6, 6, Quantity<Pascal> > GetDMatrix() const { return D_; }
        [[nodiscard]] Quantity<MeterCubed> GetVolume() const { return volume_; }

    private:
        void ComputeGeometry() {
            const auto &p0 = coords_[0];
            const auto &p1 = coords_[1];
            const auto &p2 = coords_[2];
            const auto &p3 = coords_[3];

            const auto e1 = p1 - p0;
            const auto e2 = p2 - p0;
            const auto e3 = p3 - p0;

            const auto det = e1.Dot(e2.Cross(e3));
            volume_ = det / 6.0;

            // constant for linear tetrahedron
            // N_i(x) = a_i + b_i*x + c_i*y + d_i*z
            // (b_i, c_i, d_i) / (6V)
            const auto inv6V = 1.0 / (6.0 * volume_);

            const auto g1 = (p2 - p0).Cross(p3 - p0) * inv6V;
            const auto g2 = (p3 - p0).Cross(p1 - p0) * inv6V;
            const auto g3 = (p1 - p0).Cross(p2 - p0) * inv6V;
            const auto g0 = -(g1 + g2 + g3);

            // B matrix: 6 x 12
            // rows: [exx, eyy, ezz, gxy, gyz, gzx]
            // cols: [u0x, u0y, u0z, u1x, ..., u3z]
            B_ = Matrix<6, 12, Quantity<UnitRecT<Meter>>>();

            for (int i = 0; i < 4; ++i) {
                const auto g = i == 0 ? g0 : i == 1 ? g1 : i == 2 ? g2 : g3;
                const auto c = i * 3;

                B_(0, c + 0) = g.x;
                B_(1, c + 1) = g.y;
                B_(2, c + 2) = g.z;

                B_(3, c + 0) = g.y;
                B_(3, c + 1) = g.x;

                B_(4, c + 1) = g.z;
                B_(4, c + 2) = g.y;

                B_(5, c + 0) = g.z;
                B_(5, c + 2) = g.x;
            }

            ComputeDMatrix();
        }

        void ComputeDMatrix() {
            const auto lambda = material_.GetLameLambda();
            const auto mu = material_.GetShearModulus();

            D_ = Matrix<6, 6, Quantity<Pascal> >();

            D_(0, 0) = lambda + 2.0 * mu;
            D_(0, 1) = lambda;
            D_(0, 2) = lambda;
            D_(1, 0) = lambda;
            D_(1, 1) = lambda + 2.0 * mu;
            D_(1, 2) = lambda;
            D_(2, 0) = lambda;
            D_(2, 1) = lambda;
            D_(2, 2) = lambda + 2.0 * mu;

            D_(3, 3) = mu;
            D_(4, 4) = mu;
            D_(5, 5) = mu;
        }

        std::array<Vector3<Quantity<Meter> >, 4> coords_;
        IsotropicMaterial material_;
        Matrix<6, 12, Quantity<UnitRecT<Meter>>> B_;
        Matrix<6, 6, Quantity<Pascal> > D_;
        Quantity<MeterCubed> volume_;
    };
}

#endif // PHYTH_TETRAHEDRON_HPP
