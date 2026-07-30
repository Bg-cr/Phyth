#ifndef PHYTH_VECTOR3_HPP
#define PHYTH_VECTOR3_HPP

#include "../Core/Dimension.hpp"
#include "../Core/Quantity.hpp"
#include "QuantityFuncs.hpp"
#include "../Core/Units.hpp"

#include <ostream>

namespace Phyth {
    template<typename T>
    struct Vector3 {
        static_assert(is_quantity_v<T>, "Vector3 only supports Quantity or arithmetic types!");
        T x, y, z;

        constexpr Vector3() : x(0), y(0), z(0) {
        }

        constexpr Vector3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {
        }

        static constexpr Vector3 FromPolar(T r, const Scalar theta_rad) {
            return Vector3(
                r * Utils::cos(theta_rad),
                r * Utils::sin(theta_rad),
                T(0)
            );
        }

        static constexpr Vector3 FromSpherical(T r, const Scalar theta_rad, const Scalar phi_rad) {
            return Vector3(
                r * Utils::sin(theta_rad) * Utils::cos(phi_rad),
                r * Utils::sin(theta_rad) * Utils::sin(phi_rad),
                r * Utils::cos(theta_rad)
            );
        }

        [[nodiscard]] auto ToPolar() const {
            Vector3 r = Length();
            auto theta = Utils::atan2(y / r, x / r);
            return std::pair<decltype(r), Scalar>{r, theta};
        }

        [[nodiscard]] auto ToSpherical() const {
            auto r = Length();
            auto theta = Utils::acos(z / r);
            auto phi = Utils::atan2(y / r, x / r);
            return std::tuple<decltype(r), Scalar, Scalar>{r, theta, phi};
        }

        constexpr T &operator[](int i) {
            return (&x)[i];
        }

        constexpr const T &operator[](int i) const {
            return (&x)[i];
        }

        template<typename U>
        constexpr auto operator+(const Vector3<U> &other) const {
            return Vector3<decltype(x + other.x)>(
                x + other.x,
                y + other.y,
                z + other.z
            );
        }

        template<typename U>
        constexpr auto operator-(const Vector3<U> &other) const {
            return Vector3<decltype(x - other.x)>(
                x - other.x,
                y - other.y,
                z - other.z
            );
        }

        constexpr auto operator-() const {
            return Vector3(-x, -y, -z);
        }

        template<typename U>
        constexpr auto operator*(U scalar) const {
            return Vector3<decltype(x * scalar)>(
                x * scalar,
                y * scalar,
                z * scalar
            );
        }

        template<typename U>
        constexpr auto operator/(U scalar) const {
            return Vector3<decltype(x / scalar)>(
                x / scalar,
                y / scalar,
                z / scalar
            );
        }

        template<typename U>
        constexpr Vector3 &operator+=(const Vector3<U> &other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        template<typename U>
        constexpr Vector3 &operator-=(const Vector3<U> &other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        template<typename U>
        constexpr Vector3 &operator*=(U scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        template<typename U>
        constexpr Vector3 &operator/=(U scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        template<typename U>
        constexpr bool operator==(const Vector3<U> &other) const {
            return x == other.x && y == other.y && z == other.z;
        }

        template<typename U>
        constexpr bool operator!=(const Vector3<U> &other) const {
            return !(*this == other);
        }

        template<typename U>
        constexpr auto Dot(const Vector3<U> &other) const -> decltype(x * other.x) {
            return x * other.x + y * other.y + z * other.z;
        }

        template<typename U>
        constexpr auto Cross(const Vector3<U> &other) const {
            return Vector3<decltype(x * other.x)>(
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            );
        }

        [[nodiscard]] auto Length() const {
            return Utils::sqrt(x * x + y * y + z * z);
        }

        [[nodiscard]] auto LengthSquared() const {
            return x * x + y * y + z * z;
        }

        [[nodiscard]] Vector3<Scalar> Normalized() const {
            auto len = Length();
            return Vector3<decltype(x / len)>(
                x / len,
                y / len,
                z / len
            );
        }

        template<typename U>
        [[nodiscard]] auto Project(const Vector3<U> &other) const {
            return Dot(other) / other.LengthSquared() * other;
        }

        friend std::ostream &operator<<(std::ostream &os, const Vector3 &v) {
            os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
            return os;
        }
    };

    template <typename>
    struct is_vector3 : std::false_type {};

    template <typename QuantityT>
    struct is_vector3<Vector3<QuantityT>> : std::true_type {};

    template <typename VecT>
    inline constexpr auto is_vector3_v = is_vector3<VecT>::value;

    template<typename T, typename U,
        typename = std::enable_if_t<!is_vector3_v<T>> >
    constexpr auto operator*(T scalar, const Vector3<U> &v) {
        return Vector3<decltype(scalar * v.x)>(
            scalar * v.x,
            scalar * v.y,
            scalar * v.z
        );
    }

    template<typename Dimension>
    using Vec3 = Vector3<Quantity<Dimension> >;

    using Vec3Len = Vec3<Meter>;
    using Vec3Vel = Vec3<MeterPerSecond>;
    using Vec3Acc = Vec3<MeterPerSecondSquared>;
}

#endif //PHYTH_VECTOR3_HPP
