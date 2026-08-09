#ifndef PHYTH_FIELD_HPP
#define PHYTH_FIELD_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Tools/Vector3.hpp"

#include <memory>
#include <unordered_set>

namespace Phyth {
    /**
     * @brief Base class for physical fields (scalar or vector)
     *
     * Field is an abstract interface for any physical quantity that varies
     * over space, such as:
     *   - Scalar fields: temperature, pressure, electric potential
     *   - Vector fields: electric field, magnetic field, velocity field
     *
     * @tparam T The value type returned by the field (e.g., Quantity<Voltage>
     *           for scalar potential, or Vector3<Quantity<ElectricFieldIntensity>>
     *           for vector fields)
     *
     * Example:
     *   class ElectricField : public Field<Vector3<Quantity<ElectricFieldIntensity>>> {
     *       [[nodiscard]] value_type GetValueAt(const Vector3<Quantity<Meter>>& point) const override {
     *           // Compute E-field at the given position
     *           return {0_V/m, 0_V/m, -1.0_V/m}; // uniform field for demo
     *       }
     *   };
     *
     *   class TemperatureField : public Field<Quantity<Kelvin>> {
     *       [[nodiscard]] value_type GetValueAt(const Vector3<Quantity<Meter>>& point) const override {
     *           return 300.0_K; // uniform temperature
     *       }
     *   };
     */
    template<typename T>
    class Field {
    public:
        using value_type = T;

        virtual ~Field() = default;

        /**
         * @brief Get the field value at a specific point in space
         *
         * @param point Position in 3D space (Meters)
         * @return The field value at that point
         *
         * This is the core method of any field. Implementations must handle
         * coordinate bounds and interpolation as needed.
         */
        [[nodiscard]] virtual T GetValueAt(const Vector3<Quantity<Meter> > &point) const = 0;

        /**
         * @brief Update the field state (for time-varying fields)
         *
         * This method can be overridden to advance time-dependent fields
         * (e.g., electromagnetic waves, changing boundary conditions).
         * The default implementation does nothing.
         */
        virtual void Update() {
        }
    };

    /**
     * @brief Implementation helper for is_field trait
     *
     * Detects if a type inherits from Field<something>.
     */
    template<typename, typename = void>
    struct is_field_impl : std::false_type {
    };

    template<typename T>
    struct is_field_impl<T, std::void_t<
                decltype(static_cast<Field<typename T::value_type> *>(std::declval<T *>()))
            > > : std::is_base_of<Field<typename T::value_type>, T> {
    };

    /**
     * @brief Type trait to check if a type is a Field instantiation
     *
     * Example:
     *   class MyField : public Field<Quantity<Meter>> {};
     *   static_assert(is_field_v<MyField>);  // true
     *   static_assert(is_field_v<int>);      // false
     */
    template<typename T>
    struct is_field : is_field_impl<T> {
    };

    /**
     * @brief Convenience variable template for is_field
     */
    template<typename T>
    inline constexpr bool is_field_v = is_field<T>::value;

    /**
     * @brief Combine multiple fields of the same type
     *
     * CombineField aggregates multiple Field objects and presents them as
     * a single field. The combined behavior is determined by subclasses
     * that override GetValueAt.
     *
     * @tparam FieldT The field type to combine (must inherit from Field)
     *
     * Example:
     *   class CombinedElectricField : public CombineField<ElectricField> {
     *       [[nodiscard]] value_type GetValueAt(const Vector3<Quantity<Meter>>& point) const override {
     *           value_type total = 0_V/m;
     *           for (const auto& field : fields_) {
     *               total += field->GetValueAt(point);
     *           }
     *           return total;
     *       }
     *   };
     *   CombinedElectricField combined;
     *   combined.AddField(pointChargeField);
     *   combined.AddField(uniformField);
     */
    template<typename FieldT>
    class CombineField : public Field<typename FieldT::value_type> {
        static_assert(is_field_v<FieldT>, "FieldT must be a Field.");

    public:
        /**
         * @brief Add a field to the combination
         *
         * @param field Shared pointer to the field to add
         *
         * Fields are stored in an unordered_set, so duplicates are prevented
         * and ownership is shared.
         */
        void AddField(std::shared_ptr<FieldT> field) {
            fields_.insert(std::move(field));
        }

    protected:
        std::unordered_set<std::shared_ptr<FieldT> > fields_;
    };
}

#endif // PHYTH_FIELD_HPP
