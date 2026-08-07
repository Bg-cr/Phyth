#ifndef PHYTH_FIELD_HPP
#define PHYTH_FIELD_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Tools/Vector3.hpp"

#include <memory>
#include <unordered_set>

namespace Phyth {

    /**
     * @brief Base class for physical fields (scalar or vector)
     * @tparam T The value type returned by the field
     */
    template <typename T>
    class Field {
    public:
        using value_type = T;

        virtual ~Field() = default;

        [[nodiscard]] virtual T GetValueAt(const Vector3<Quantity<Meter>>& point) const = 0;
        virtual void Update() {}
    };

    template <typename, typename = void>
    struct is_field_impl : std::false_type {};

    template <typename T>
    struct is_field_impl<T, std::void_t<
        decltype(static_cast<Field<typename T::value_type>*>(std::declval<T*>()))
    >> : std::is_base_of<Field<typename T::value_type>, T> {};

    /**
     * @brief Type trait to check if T is a Field
     */
    template <typename T>
    struct is_field : is_field_impl<T> {};

    template <typename T>
    inline constexpr bool is_field_v = is_field<T>::value;

    /**
     * @brief Manage multiple fields of the same type
     * @tparam FieldT The field type to combine
     */
    template <typename FieldT>
    class CombineField : public Field<typename FieldT::value_type> {
        static_assert(is_field_v<FieldT>, "FieldT must be a Field.");
    public:
        void AddField(std::shared_ptr<FieldT> field) {
            fields_.insert(std::move(field));
        }
    protected:
        std::unordered_set<std::shared_ptr<FieldT>> fields_;
    };

} // namespace Phyth

#endif //PHYTH_FIELD_HPP