#ifndef PHYTH_MAGNETIC_FIELD_HPP
#define PHYTH_MAGNETIC_FIELD_HPP

#include <algorithm>
#include <memory>

#include "FieldSources/MagneticSource.hpp"
#include "Phyth/Elements/Field.hpp"

namespace Phyth::Electromagnetics {
    /**
     * @brief Magnetic field composed of multiple magnetic sources
     *
     * MagneticField aggregates multiple MagneticSource objects and computes
     * the total magnetic field at any point by summing their contributions.
     *
     * This class implements the Field interface, returning a Vector3 of
     * Tesla quantities (magnetic flux density B).
     *
     * Example:
     *   MagneticField field;
     *
     *   auto dipole = std::make_shared<MagneticDipole>(...);
     *   field.AddSource(dipole);
     *
     *   auto point = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     *   auto B = field.GetValueAt(point);  // Total B-field at that point
     */
    class MagneticField : public Field<Vector3<Quantity<Tesla> > > {
    public:
        MagneticField() = default;

        /**
         * @brief Add a single magnetic source to the field
         *
         * @param source Shared pointer to the source to add
         *
         * Duplicate sources are automatically handled by the unordered_set.
         * Null pointers are silently ignored.
         */
        void AddSource(std::shared_ptr<MagneticSource> source) {
            if (source) {
                sources_.insert(std::move(source));
            }
        }

        /**
         * @brief Add multiple magnetic sources
         *
         * @tparam Args Parameter pack of shared_ptr<MagneticSource>
         * @param args The sources to add
         *
         * Example:
         *   field.AddSources(dipole1, dipole2, coil);
         */
        template<typename... Args>
        void AddSources(Args &&... args) {
            (AddSource(std::forward<Args>(args)), ...);
        }

        /**
         * @brief Compute the total magnetic field at a point
         *
         * Sums the magnetic field contribution from all sources.
         *
         * @param point Position in 3D space (meters)
         * @return Total magnetic flux density B (Tesla) as a Vector3
         */
        [[nodiscard]] Vector3<Quantity<Tesla> >
        GetValueAt(const Vector3<Quantity<Meter> > &point) const override {
            Vector3 total{0_T, 0_T, 0_T};
            for (const auto &source: sources_) {
                total += source->GetMagneticFieldAt(point);
            }
            return total;
        }

        /**
         * @brief Get all sources of a specific type
         *
         * @tparam T The magnetic source type to filter by (must satisfy is_magnetic_source_v)
         *
         * @return A set of shared_ptr<T> containing all sources that are of type T
         *
         * Example:
         *   auto dipoles = field.GetMagneticSources<MagneticDipole>();
         *   // Only returns sources that are MagneticDipole instances
         */
        template<typename T, typename = std::enable_if<is_magnetic_source_v<T> > >
        [[nodiscard]] std::unordered_set<std::shared_ptr<T> >
        GetMagneticSources() const {
            std::unordered_set<std::shared_ptr<T> > sources;
            for (const auto &source: sources_) {
                if (auto transformed_source = std::dynamic_pointer_cast<T>(source)) {
                    sources.insert(transformed_source);
                }
            }
            return sources;
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
        void RemoveSource(const std::shared_ptr<MagneticSource> &source) {
            if (sources_.find(source) != sources_.end()) {
                sources_.erase(source);
            }
        }

        /**
         * @brief Remove multiple sources
         *
         * @tparam Args Parameter pack of shared_ptr<MagneticSource>
         * @param args The sources to remove
         */
        template<typename... Args>
        void RemoveSources(Args &&... args) {
            (RemoveSource(std::forward<Args>(args)), ...);
        }

        /** @brief Get all sources currently in the field */
        [[nodiscard]] std::unordered_set<std::shared_ptr<MagneticSource> > GetSources() const {
            return sources_;
        }

    private:
        std::unordered_set<std::shared_ptr<MagneticSource> > sources_;
    };
}

#endif // PHYTH_MAGNETIC_FIELD_HPP
