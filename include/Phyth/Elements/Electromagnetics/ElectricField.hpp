#ifndef PHYTH_ELECTRIC_FIELD_HPP
#define PHYTH_ELECTRIC_FIELD_HPP

#include <algorithm>
#include <memory>
#include <vector>

#include "ChargeSources/ChargeSources.hpp"
#include "Phyth/Elements/Field.hpp"

namespace Phyth::Electromagnetics {
    class ElectricField : public Field<Vector3<Quantity<NewtonPerCoulomb> > > {
    public:
        ElectricField() = default;

        void AddSource(std::shared_ptr<ChargeSource> source) {
            if (source) {
                sources_.insert(std::move(source));
            }
        }

        template<typename... Args>
        void AddSources(Args &&... args) {
            (AddSource(std::forward<Args>(args)), ...);
        }

        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetValueAt(const Vector3<Quantity<Meter> > &point) const override {
            Vector3 total{0_NpC, 0_NpC, 0_NpC};
            for (const auto &source: sources_) {
                total += source->GetElectricFieldAt(point);
            }
            return total;
        }

        [[nodiscard]] Quantity<Volt>
        GetElectricPotential(const Vector3<Quantity<Meter> > &point) const {
            auto total = 0_V;
            for (const auto &source: sources_) {
                total += source->GetElectricPotentialAt(point);
            }
            return total;
        }

        [[nodiscard]] Quantity<CoulombPerCubicMeter>
        GetChargeDensityAt(const Vector3<Quantity<Meter> > &point) const {
            auto total = 0_Cpm3;
            for (const auto &source: sources_) {
                total += source->GetChargeDensityAt(point);
            }
            return total;
        }

        template <typename T, std::enable_if<is_charge_source_v<T>>>
        [[nodiscard]] std::unordered_set<std::shared_ptr<T>>
        GetChargeSources() const {
            std::unordered_set<std::shared_ptr<T> > sources;
            for (const auto &source: sources_) {
                if (auto transformed_source = std::dynamic_pointer_cast<T>(source)) {
                    sources.insert(transformed_source);
                }
            }
            return sources;
        }

        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const {
            auto total = 0_C;
            for (const auto &source: sources_) {
                total += source->GetTotalChargeValue();
            }
            return total;
        }

        void Clear() {
            sources_.clear();
        }

        void RemoveSource(const std::shared_ptr<ChargeSource> &source) {
            if (sources_.find(source) != sources_.end()) {
                sources_.erase(source);
            }
        }

        template <typename... Args>
        void RemoveSources(Args &&... args) {
            (RemoveSource(std::forward<Args>(args)), ...);
        }

        [[nodiscard]] std::unordered_set<std::shared_ptr<ChargeSource>> GetSources() const {
            return sources_;
        }

        std::string toString() const {
            return "ElectricField{"
                   "total_charge=" + QuantityToString(GetTotalChargeValue()) +
                   "}";
        }

    private:
        std::unordered_set<std::shared_ptr<ChargeSource> > sources_;
    };

    inline std::ostream &operator<<(std::ostream &os, const ElectricField &field) {
        os << field.toString();
        return os;
    }
}

#endif
