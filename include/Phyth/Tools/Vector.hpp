#ifndef PHYTH_VECTOR_HPP
#define PHYTH_VECTOR_HPP

#include "Phyth/Core/Quantity.hpp"
#include "QuantityFuncs.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <stdexcept>
#include <ostream>
#include <type_traits>

namespace Phyth {
    template<std::size_t N, typename T>
    struct Vector {
        static_assert(is_quantity_v<T>, "Vector only supports Quantity types");
        static_assert(N != 0, "Does not support zero dimensional space");

        std::array<T, N> datas;

        constexpr Vector() = default;

        template<typename... Args,
            std::enable_if_t<
                sizeof...(Args) == N &&
                std::conjunction_v<std::is_convertible<Args, T>...>,
                int> = 0>
        constexpr Vector(Args... comps) : datas{comps...} {
        }

        constexpr T &At(const size_t i) {
            return datas.at(i);
        }

        constexpr const T &At(const size_t i) const {
            return datas.at(i);
        }

        constexpr T &operator[](const size_t i) {
            return datas[i];
        }

        constexpr const T &operator[](const size_t i) const {
            return datas[i];
        }

        template<typename U>
        constexpr auto operator+(const Vector<N, U> &other) const {
            using ResultScalar = decltype(datas[0] + other.datas[0]);
            return make_binary<ResultScalar>(
                other, std::make_index_sequence<N>{},
                [](const auto &a, const auto &b) { return a + b; });
        }

        template<typename U>
        constexpr auto operator-(const Vector<N, U> &other) const {
            using ResultScalar = decltype(datas[0] - other.datas[0]);
            return make_binary<ResultScalar>(
                other, std::make_index_sequence<N>{},
                [](const auto &a, const auto &b) { return a - b; });
        }

        constexpr auto operator-() const {
            using ResultScalar = decltype(-datas[0]);
            return make_unary<ResultScalar>(
                std::make_index_sequence<N>{},
                [](const auto &a) { return -a; });
        }

        template<typename U,
            typename = std::enable_if_t<!std::is_same_v<std::decay_t<U>, Vector>> >
        constexpr auto operator*(U scalar) const {
            using ResultScalar = decltype(datas[0] * scalar);
            return make_scalar<ResultScalar>(
                scalar, std::make_index_sequence<N>{},
                [](const auto &a, const auto &s) { return a * s; });
        }

        template<typename U>
        constexpr auto operator/(U scalar) const {
            using ResultScalar = decltype(datas[0] / scalar);
            return make_scalar<ResultScalar>(
                scalar, std::make_index_sequence<N>{},
                [](const auto &a, const auto &s) { return a / s; });
        }

        template<typename U>
        constexpr Vector &operator+=(const Vector<N, U> &other) {
            for (std::size_t i = 0; i < N; ++i)
                datas[i] += other.datas[i];
            return *this;
        }

        template<typename U>
        constexpr Vector &operator-=(const Vector<N, U> &other) {
            for (std::size_t i = 0; i < N; ++i)
                datas[i] -= other.datas[i];
            return *this;
        }

        template<typename U>
        constexpr Vector &operator*=(U scalar) {
            for (std::size_t i = 0; i < N; ++i)
                datas[i] *= scalar;
            return *this;
        }

        template<typename U>
        constexpr Vector &operator/=(U scalar) {
            for (std::size_t i = 0; i < N; ++i)
                datas[i] /= scalar;
            return *this;
        }

        constexpr bool operator==(const Vector &other) const {
            for (std::size_t i = 0; i < N; ++i)
                if (datas[i] != other.datas[i])
                    return false;
            return true;
        }

        constexpr bool operator!=(const Vector &other) const {
            return !(*this == other);
        }

        [[nodiscard]] auto LengthSquared() const {
            using ResultScalar = decltype(datas[0] * datas[0]);
            ResultScalar result{0};
            for (std::size_t i = 0; i < N; ++i)
                result += datas[i] * datas[i];
            return result;
        }

        [[nodiscard]] auto Length() const {
            return Utils::sqrt(LengthSquared());
        }

        [[nodiscard]] auto Normalized() const {
            auto len = Length();
            if (len == decltype(len){0}) {
                throw std::runtime_error("Normalize the zero vector");
            }
            using ResultScalar = decltype(datas[0] / len);
            return make_scalar<ResultScalar>(
                len, std::make_index_sequence<N>{},
                [](const auto &a, const auto &l) { return a / l; });
        }

        template<typename U>
        [[nodiscard]] constexpr auto Dot(const Vector<N, U> &other) const {
            using ResultScalar = decltype(datas[0] * other.datas[0]);
            ResultScalar result{0};
            for (std::size_t i = 0; i < N; ++i)
                result += datas[i] * other.datas[i];
            return result;
        }

        template<typename U>
        [[nodiscard]] auto Project(const Vector<N, U> &other) const {
            return Dot(other) / other.LengthSquared() * other;
        }

        template<typename UnitT>
        [[nodiscard]] auto As() const {
            using ResultScalar = decltype(datas[0].template As<UnitT>());
            return make_as<ResultScalar, UnitT>(std::make_index_sequence<N>{});
        }

        friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
            os << "(";
            for (std::size_t i = 0; i < N; ++i) {
                if (i > 0) os << ", ";
                os << v.datas[i];
            }
            os << ")";
            return os;
        }

    private:
        template<typename ResultScalar, typename U, std::size_t... I, typename F>
        constexpr auto make_binary(const Vector<N, U> &other,
                                   std::index_sequence<I...>, F f) const {
            return Vector<N, ResultScalar>(f(datas[I], other.datas[I])...);
        }

        template<typename ResultScalar, std::size_t... I, typename F>
        constexpr auto make_unary(std::index_sequence<I...>, F f) const {
            return Vector<N, ResultScalar>(f(datas[I])...);
        }

        template<typename ResultScalar, typename U, std::size_t... I, typename F>
        constexpr auto make_scalar(U scalar,
                                   std::index_sequence<I...>, F f) const {
            return Vector<N, ResultScalar>(f(datas[I], scalar)...);
        }

        template<typename ResultScalar, typename UnitT, std::size_t... I>
        constexpr auto make_as(std::index_sequence<I...>) const {
            return Vector<N, ResultScalar>(datas[I].template As<UnitT>()...);
        }
    };

    template<typename>
    struct is_vector : std::false_type {
    };

    template<std::size_t N, typename QuantityT>
    struct is_vector<Vector<N, QuantityT> > : std::true_type {
    };

    template<typename VecT>
    inline constexpr auto is_vector_v = is_vector<VecT>::value;

    template<std::size_t N, typename T>
    struct is_vector_of_n : std::false_type {
    };

    template<std::size_t N, typename T>
    struct is_vector_of_n<N, Vector<N, T> > : std::true_type {
    };

    template<std::size_t N, typename T>
    inline constexpr bool is_vector_of_n_v = is_vector_of_n<N, T>::value;

    template<typename T, std::size_t N, typename U,
        typename = std::enable_if_t<!is_vector_of_n_v<N, T>> >
    constexpr auto operator*(T scalar, const Vector<N, U> &v) {
        return v * scalar;
    }
}

#endif //PHYTH_VECTOR_HPP
