#ifndef PHYTH_SURFACE_CHARGE_HPP
#define PHYTH_SURFACE_CHARGE_HPP

#include <cassert>

#include "ChargeSource.hpp"
#include "Phyth/Physical/PhysicalConsts.hpp"

namespace Phyth::Electromagnetics {
    class SurfaceCharge : public ChargeSource {
    public:
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

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetElectricFieldAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegrateField(point, 32, 32);
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter> > &point) const override {
            return NumericalIntegratePotential(point, 32, 32);
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const override {
            return sigma_ * width_ * height_;
        }

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

        [[nodiscard]] Quantity<CoulombPerSquareMeter>
        GetSurfaceDensity() const { return sigma_; }

    private:
        Vector3<Quantity<Meter> > center_;
        Vector3<Scalar> u_axis_, v_axis_;
        Quantity<Meter> width_, height_;
        Quantity<CoulombPerSquareMeter> sigma_;

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

#endif //PHYTH_SURFACE_CHARGE_HPP
