#ifndef PHYTH_ELECTRIC_FIELD_HPP
#define PHYTH_ELECTRIC_FIELD_HPP

#include <algorithm>
#include <memory>

#include "FieldSources/ChargeSource.hpp"
#include "Phyth/Elements/Field.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Electric field composed of multiple charge sources
     *
     * ElectricField aggregates multiple ChargeSource objects and computes
     * the total electric field, potential, and charge density at any point
     * by summing their contributions.
     *
     * This class implements the Field interface, returning a Vector3 of
     * NewtonPerCoulomb quantities (electric field intensity E).
     *
     * Example:
     *   ElectricField field;
     *
     *   auto point_charge = std::make_shared<PointCharge>(1.0_C, Vector3{0_m, 0_m, 0_m});
     *   field.AddSource(point_charge);
     *
     *   auto point = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     *   auto E = field.GetValueAt(point);           // Total E-field
     *   auto V = field.GetElectricPotential(point); // Total potential
     */
    class ElectricField : public Field<Vector3<Quantity<NewtonPerCoulomb> > > {
    public:
        ElectricField() = default;

        /**
         * @brief Add a single charge source to the field
         *
         * @param source Shared pointer to the source to add
         *
         * Duplicate sources are automatically handled by the unordered_set.
         * Null pointers are silently ignored.
         */
        void AddSource(std::shared_ptr<ChargeSource> source) {
            if (source) {
                sources_.insert(std::move(source));
            }
        }

        /**
         * @brief Add multiple charge sources
         *
         * @tparam Args Parameter pack of shared_ptr<ChargeSource>
         * @param args The sources to add
         *
         * Example:
         *   field.AddSources(charge1, charge2, dipole);
         */
        template<typename... Args>
        void AddSources(Args &&... args) {
            (AddSource(std::forward<Args>(args)), ...);
        }

        /**
         * @brief Compute the total electric field at a point
         *
         * Sums the electric field contribution from all sources.
         *
         * @param point Position in 3D space (meters)
         * @return Total electric field intensity E (N/C) as a Vector3
         */
        [[nodiscard]] Vector3<Quantity<NewtonPerCoulomb> >
        GetValueAt(const Vector3<Quantity<Meter> > &point) const override {
            Vector3 total{0_NpC, 0_NpC, 0_NpC};
            for (const auto &source: sources_) {
                total += source->GetElectricFieldAt(point);
            }
            return total;
        }

        /**
         * @brief Compute the total electric potential at a point
         *
         * Sums the electric potential contribution from all sources.
         *
         * @param point Position in 3D space (meters)
         * @return Total electric potential V (Volts)
         */
        [[nodiscard]] Quantity<Volt>
        GetElectricPotential(const Vector3<Quantity<Meter> > &point) const {
            auto total = 0_V;
            for (const auto &source: sources_) {
                total += source->GetElectricPotentialAt(point);
            }
            return total;
        }

        /**
         * @brief Compute the total charge density at a point
         *
         * Sums the charge density contribution from all sources.
         *
         * @param point Position in 3D space (meters)
         * @return Total charge density (C/m^3)
         */
        [[nodiscard]] Quantity<CoulombPerCubicMeter>
        GetChargeDensityAt(const Vector3<Quantity<Meter> > &point) const {
            auto total = 0_Cpm3;
            for (const auto &source: sources_) {
                total += source->GetChargeDensityAt(point);
            }
            return total;
        }

        /**
         * @brief Get all sources of a specific charge source type
         *
         * @tparam T The charge source type to filter by (must satisfy is_charge_source_v)
         * @return A set of shared_ptr<T> containing all sources that are of type T
         *
         * Example:
         *   auto point_charges = field.GetChargeSources<PointCharge>();
         *   // Only returns sources that are PointCharge instances
         */
        template<typename T, typename = std::enable_if<is_charge_source_v<T> > >
        [[nodiscard]] std::unordered_set<std::shared_ptr<T> >
        GetChargeSources() const {
            std::unordered_set<std::shared_ptr<T> > sources;
            for (const auto &source: sources_) {
                if (auto transformed_source = std::dynamic_pointer_cast<T>(source)) {
                    sources.insert(transformed_source);
                }
            }
            return sources;
        }

        /**
         * @brief Get the total charge of all sources combined
         *
         * Sums the total charge value from each source.
         *
         * @return Total charge (Coulombs)
         */
        [[nodiscard]] Quantity<Coulomb> GetTotalChargeValue() const {
            auto total = 0_C;
            for (const auto &source: sources_) {
                total += source->GetTotalChargeValue();
            }
            return total;
        }

        /** @brief Remove all sources from the field */
        void Clear() {
            sources_.clear();
        }

        /**
         * @brief Remove a specific source from the field
         *
         * @param source Shared pointer to the source to remove
         *
         * If the source is not present, this function does nothing.
         */
        void RemoveSource(const std::shared_ptr<ChargeSource> &source) {
            if (sources_.find(source) != sources_.end()) {
                sources_.erase(source);
            }
        }

        /**
         * @brief Remove multiple sources
         *
         * @tparam Args Parameter pack of shared_ptr<ChargeSource>
         * @param args The sources to remove
         */
        template<typename... Args>
        void RemoveSources(Args &&... args) {
            (RemoveSource(std::forward<Args>(args)), ...);
        }

        /** @brief Get all sources currently in the field */
        [[nodiscard]] std::unordered_set<std::shared_ptr<ChargeSource> > GetSources() const {
            return sources_;
        }

        /**
         * @brief Convert the field to a string representation
         *
         * @return String containing the total charge value
         */
        std::string toString() const {
            return "ElectricField{"
                   "total_charge=" + QuantityToString(GetTotalChargeValue()) +
                   "}";
        }

    private:
        std::unordered_set<std::shared_ptr<ChargeSource> > sources_;
    };

    /**
     * @brief Stream insertion operator for ElectricField
     *
     * Output format: "ElectricField{total_charge=X C}"
     */
    inline std::ostream &operator<<(std::ostream &os, const ElectricField &field) {
        os << field.toString();
        return os;
    }
}

#endif // PHYTH_ELECTRIC_FIELD_HPP
