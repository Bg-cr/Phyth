#ifndef PHYTH_MAGNETIC_FIELD_HPP
#define PHYTH_MAGNETIC_FIELD_HPP

#include <algorithm>
#include <memory>

#include "FieldSources/MagneticSource.hpp"
#include "Phyth/Elements/Field.hpp"

namespace Phyth::Electromagnetics {
    class MagneticField : public Field<Vector3<Quantity<Tesla>>> {
    public:
        MagneticField() = default;

        void AddSource(std::shared_ptr<MagneticSource> source) {
            if (source) {
                sources_.insert(std::move(source));
            }
        }

        template<typename... Args>
        void AddSources(Args&&... args) {
            (AddSource(std::forward<Args>(args)), ...);
        }

        [[nodiscard]] Vector3<Quantity<Tesla>>
        GetValueAt(const Vector3<Quantity<Meter>>& point) const override {
            Vector3 total{0_T, 0_T, 0_T};
            for (const auto& source : sources_) {
                total += source->GetMagneticFieldAt(point);
            }
            return total;
        }

        template <typename T, typename = std::enable_if<is_magnetic_source_v<T>>>
        [[nodiscard]] std::unordered_set<std::shared_ptr<T>>
        GetMagneticSources() const {
            std::unordered_set<std::shared_ptr<T>> sources;
            for (const auto& source : sources_) {
                if (auto transformed_source = std::dynamic_pointer_cast<T>(source)) {
                    sources.insert(transformed_source);
                }
            }
            return sources;
        }

        void Clear() {
            sources_.clear();
        }

        void RemoveSource(const std::shared_ptr<MagneticSource>& source) {
            if (sources_.find(source) != sources_.end()) {
                sources_.erase(source);
            }
        }

        template <typename... Args>
        void RemoveSources(Args&&... args) {
            (RemoveSource(std::forward<Args>(args)), ...);
        }

        [[nodiscard]] std::unordered_set<std::shared_ptr<MagneticSource>> GetSources() const {
            return sources_;
        }

    private:
        std::unordered_set<std::shared_ptr<MagneticSource>> sources_;
    };

}

#endif //PHYTH_MAGNETIC_FIELD_HPP
