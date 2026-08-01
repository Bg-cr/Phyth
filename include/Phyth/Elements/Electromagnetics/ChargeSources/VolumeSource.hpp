#ifndef PHYTH_VOLUME_SOURCE_HPP
#define PHYTH_VOLUME_SOURCE_HPP

#include <cassert>

#include "ChargeSource.hpp"

namespace Phyth::Electromagnetics {
    class VolumeCharge : public ChargeSource {
    public:
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

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegrateField(point, 16, 16, 16);
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegratePotential(point, 16, 16, 16);
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return rho_ * size_x_ * size_y_ * size_z_;
        }

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

#endif //PHYTH_VOLUME_SOURCE_HPP
