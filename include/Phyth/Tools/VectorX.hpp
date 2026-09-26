#ifndef PHYTH_VECTOR_X_HPP
#define PHYTH_VECTOR_X_HPP

#include "Phyth/Core/Quantity.hpp"
#include "QuantityFuncs.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

namespace Phyth {
    template<typename T>
    class VectorX {
        static_assert(is_quantity_v<T>, "VectorX only supports Quantity types");

    public:
        VectorX() = default;

        explicit VectorX(std::size_t n, T value = T(0))
            : data_(n, value) {
        }

        explicit VectorX(std::vector<T> data)
            : data_(std::move(data)) {
        }

        [[nodiscard]] std::size_t Size() const { return data_.size(); }
        [[nodiscard]] bool Empty() const { return data_.empty(); }
        void Resize(std::size_t n, T value = T(0)) { data_.resize(n, value); }

        T &operator[](std::size_t i) { return data_[i]; }
        const T &operator[](std::size_t i) const { return data_[i]; }

        T &At(std::size_t i) {
            if (i >= data_.size()) throw std::out_of_range("VectorX::At");
            return data_[i];
        }

        const T &At(std::size_t i) const {
            if (i >= data_.size()) throw std::out_of_range("VectorX::At");
            return data_[i];
        }

        template<typename U>
        auto operator+(const VectorX<U> &other) const {
            CheckSize(other);
            using R = decltype(data_[0] + other[0]);
            VectorX<R> result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] + other[i];
            return result;
        }

        template<typename U>
        auto operator-(const VectorX<U> &other) const {
            CheckSize(other);
            using R = decltype(data_[0] - other[0]);
            VectorX<R> result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] - other[i];
            return result;
        }

        auto operator-() const {
            using R = decltype(-data_[0]);
            VectorX<R> result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = -data_[i];
            return result;
        }

        template<typename U, typename = std::enable_if_t<!is_vector_x_v<U>> >
        auto operator*(U scalar) const {
            using R = decltype(data_[0] * scalar);
            VectorX<R> result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] * scalar;
            return result;
        }

        VectorX &operator+=(const VectorX &other) {
            CheckSize(other);
            for (std::size_t i = 0; i < data_.size(); ++i)
                data_[i] += other[i];
            return *this;
        }

        VectorX &operator-=(const VectorX &other) {
            CheckSize(other);
            for (std::size_t i = 0; i < data_.size(); ++i)
                data_[i] -= other[i];
            return *this;
        }

        VectorX &operator*=(T scalar) {
            for (auto &v: data_) v *= scalar;
            return *this;
        }

        template<typename U>
        auto Dot(const VectorX<U> &other) const {
            CheckSize(other);
            using R = decltype(data_[0] * other[0]);
            R result = R(0);
            for (std::size_t i = 0; i < data_.size(); ++i)
                result += data_[i] * other[i];
            return result;
        }

        [[nodiscard]] auto Length() const {
            return Utils::sqrt(Dot(*this));
        }

        [[nodiscard]] auto LengthSquared() const {
            return Dot(*this);
        }

        [[nodiscard]] auto Normalized() const {
            auto len = Length();
            using R = decltype(data_[0] / len);
            VectorX<R> result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i] / len;
            return result;
        }

        template<typename U>
        bool operator==(const VectorX<U> &other) const {
            if (data_.size() != other.Size()) return false;
            for (std::size_t i = 0; i < data_.size(); ++i)
                if (!(data_[i] == other[i])) return false;
            return true;
        }

        template<typename U>
        bool operator!=(const VectorX<U> &other) const {
            return !(*this == other);
        }

        template<typename UnitT>
        [[nodiscard]] VectorX<Quantity<UnitT> > As() const {
            VectorX<Quantity<UnitT> > result(data_.size());
            for (std::size_t i = 0; i < data_.size(); ++i)
                result[i] = data_[i].template As<UnitT>();
            return result;
        }

        [[nodiscard]] const std::vector<T> &Data() const { return data_; }
        [[nodiscard]] std::vector<T> &Data() { return data_; }

    private:
        template<typename U>
        void CheckSize(const VectorX<U> &other) const {
            if (data_.size() != other.Size())
                throw std::invalid_argument("VectorX size mismatch");
        }

        std::vector<T> data_;
    };

    template<typename T>
    struct is_vector_x : std::false_type {
    };

    template<typename T>
    struct is_vector_x<VectorX<T> > : std::true_type {
    };

    template<typename T>
    inline constexpr bool is_vector_x_v = is_vector_x<T>::value;

    template<typename T, typename U, typename = std::enable_if_t<!is_vector_x_v<T>> >
    auto operator*(T scalar, const VectorX<U> &v) {
        return v * scalar;
    }
}

#endif // PHYTH_VECTOR_X_HPP
