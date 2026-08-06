#ifndef PHYTH_QUATERNION_HPP
#define PHYTH_QUATERNION_HPP

#include "Phyth/Core/Dimension.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Quantities.hpp"
#include "QuantityFuncs.hpp"
#include "Vector3.hpp"

#include <ostream>

namespace Phyth {
    /**
     * @brief Quaternion with Scalar components
     *
     * @tparam T Component type
     */
    template<typename T>
    struct Quaternion {
        static_assert(is_quantity_v<T>,
                      "Quaternion only supports Quantity types!");

        T w, x, y, z;

        /**
         * @brief Default constructor - identity quaternion
         */
        constexpr Quaternion() : w(1), x(0), y(0), z(0) {
        }

        /**
         * @brief Construct from scalar and vector parts
         * @param w_ Scalar part
         * @param x_ x-component of vector part
         * @param y_ y-component of vector part
         * @param z_ z-component of vector part
         */
        constexpr Quaternion(T w_, T x_, T y_, T z_) : w(w_), x(x_), y(y_), z(z_) {
        }

        /**
         * @brief Construct from scalar and vector parts
         * @param w_ Scalar part
         * @param v_ Vector part
         */
        constexpr Quaternion(T w_, const Vector3<T>& v_) : w(w_), x(v_.x), y(v_.y), z(v_.z) {
        }

        /**
         * @brief Construct from vector part only (scalar part = 0)
         * @param v_ Vector part
         */
        constexpr explicit Quaternion(const Vector3<T>& v_) : w(0), x(v_.x), y(v_.y), z(v_.z) {
        }

        /**
         * @brief Construct from axis and angle
         *
         * Creates a unit quaternion representing rotation by angle theta
         * around the given axis.
         *
         * @param axis Unit vector axis of rotation
         * @param theta Rotation angle in radians
         * @return Unit quaternion representing the rotation
         */
        static Quaternion FromAxisAngle(const Vector3<T>& axis, const Scalar theta) {
            const auto half_theta = theta / 2;
            const auto s = Utils::sin(half_theta);
            const auto c = Utils::cos(half_theta);
            return Quaternion(c, axis * s);
        }

        /**
         * @brief Construct from Euler angles (ZYX convention)
         *
         * Creates a quaternion from yaw, pitch, and roll angles.
         *
         * @param yaw Rotation around Z axis (radians)
         * @param pitch Rotation around Y axis (radians)
         * @param roll Rotation around X axis (radians)
         * @return Unit quaternion representing the combined rotation
         */
        static Quaternion FromEuler(const Scalar yaw, const Scalar pitch, const Scalar roll) {
            const auto cy = Utils::cos(yaw / 2);
            const auto sy = Utils::sin(yaw / 2);
            const auto cp = Utils::cos(pitch / 2);
            const auto sp = Utils::sin(pitch / 2);
            const auto cr = Utils::cos(roll / 2);
            const auto sr = Utils::sin(roll / 2);

            return Quaternion(
                cr * cp * cy + sr * sp * sy,
                sr * cp * cy - cr * sp * sy,
                cr * sp * cy + sr * cp * sy,
                cr * cp * sy - sr * sp * cy
            );
        }

        [[nodiscard]] constexpr Vector3<T> ToVector() const {
            return Vector3<T>(x, y, z);
        }

        [[nodiscard]] constexpr T ToQuantity() const {
            return w;
        }

        [[nodiscard]] auto LengthSquared() const {
            return w * w + x * x + y * y + z * z;
        }

        [[nodiscard]] auto Length() const {
            return Utils::sqrt(LengthSquared());
        }

        [[nodiscard]] constexpr Quaternion Conjugated() const {
            return Quaternion(w, -x, -y, -z);
        }

        [[nodiscard]] auto Inversed() const {
            const auto length_sq = LengthSquared();
            const auto conj = Conjugated();
            return Quaternion<decltype(w / length_sq)>(
                conj.w / length_sq,
                conj.x / length_sq,
                conj.y / length_sq,
                conj.z / length_sq
            );
        }

        [[nodiscard]] auto Normalized() const {
            const auto length = Length();
            return Quaternion<decltype(w / length)>(
                w / length,
                x / length,
                y / length,
                z / length
            );
        }

        [[nodiscard]] bool IsUnit() const {
            return Utils::abs(Scalar(1) - LengthSquared()) < Config::tolerance;
        }

        [[nodiscard]] auto Rotated(const Vector3<T>& v) const {
            const auto qv = *this * Quaternion(v);
            const auto result = qv * Inversed();
            return result.ToVector();
        }

        template<typename U>
        constexpr auto operator*(const Quaternion<U>& other) const {
            using R = decltype(w * other.w);
            return Quaternion<R>(
                w * other.w - x * other.x - y * other.y - z * other.z,
                w * other.x + x * other.w + y * other.z - z * other.y,
                w * other.y - x * other.z + y * other.w + z * other.x,
                w * other.z + x * other.y - y * other.x + z * other.w
            );
        }

        template<typename U>
        constexpr auto operator+(const Quaternion<U>& other) const {
            return Quaternion<decltype(w + other.w)>(
                w + other.w,
                x + other.x,
                y + other.y,
                z + other.z
            );
        }

        template<typename U>
        constexpr auto operator-(const Quaternion<U>& other) const {
            return Quaternion<decltype(w - other.w)>(
                w - other.w,
                x - other.x,
                y - other.y,
                z - other.z
            );
        }

        constexpr auto operator-() const {
            return Quaternion(-w, -x, -y, -z);
        }

        template<typename U>
        constexpr auto operator*(Quantity<U> scalar) const {
            return Quaternion<decltype(w * scalar)>(
                w * scalar,
                x * scalar,
                y * scalar,
                z * scalar
            );
        }

        template<typename U>
        constexpr auto operator/(U scalar) const {
            return Quaternion<decltype(w / scalar)>(
                w / scalar,
                x / scalar,
                y / scalar,
                z / scalar
            );
        }

        template<typename U>
        constexpr Quaternion& operator+=(const Quaternion<U>& other) {
            w += other.w;
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        template<typename U>
        constexpr Quaternion& operator-=(const Quaternion<U>& other) {
            w -= other.w;
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        template<typename U>
        constexpr Quaternion& operator*=(U scalar) {
            w *= scalar;
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }


        template<typename U>
        constexpr Quaternion& operator/=(U scalar) {
            w /= scalar;
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        template<typename U>
        constexpr bool operator==(const Quaternion<U>& other) const {
            return w == other.w && x == other.x && y == other.y && z == other.z;
        }

        template<typename U>
        constexpr bool operator!=(const Quaternion<U>& other) const {
            return !(*this == other);
        }

        template <typename UnitT>
        [[nodiscard]] Quaternion<Quantity<UnitT>> as() const {
            return {
                w.template as<UnitT>(),
                x.template as<UnitT>(),
                y.template as<UnitT>(),
                z.template as<UnitT>()
            };
        }

        template<typename U>
        static auto Slerp(const Quaternion& q1, const Quaternion<U>& q2, const Scalar t) {
            auto cos_angle = q1.Dot(q2);
            auto q2_adjusted = q2;
            if (cos_angle < 0) {
                q2_adjusted = -q2_adjusted;
                cos_angle = -cos_angle;
            }

            const auto cos_angle_clamped = Utils::clamp(cos_angle, -1_, 1_);
            const auto angle = Utils::acos(cos_angle_clamped);

            if (angle < 1e-6_rad) {
                return (q1 * (1_ - t) + q2_adjusted * t).Normalized();
            }

            const auto sin_angle = Utils::sin(angle);
            const auto a = Utils::sin((1_ - t) * angle) / sin_angle;
            const auto b = Utils::sin(t * angle) / sin_angle;

            return (q1 * a + q2_adjusted * b).Normalized();
        }

        template<typename U>
        [[nodiscard]] constexpr auto Dot(const Quaternion<U>& other) const {
            return w * other.w + x * other.x + y * other.y + z * other.z;
        }

        constexpr T& operator[](int i) {
            return (&w)[i];
        }

        constexpr const T& operator[](int i) const {
            return (&w)[i];
        }

        friend std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
            os << "(" << q.w << " + " << q.x << "i + " << q.y << "j + " << q.z << "k)";
            return os;
        }
    };

    template<typename>
    struct is_quaternion : std::false_type {};

    template<typename T>
    struct is_quaternion<Quaternion<T>> : std::true_type {};

    template<typename T>
    inline constexpr auto is_quaternion_v = is_quaternion<T>::value;

    template<typename T, typename U>
    constexpr auto operator*(const Quantity<T> scalar, const Quaternion<U>& q) {
        return Quaternion<decltype(scalar * q.w)>(
            scalar * q.w,
            scalar * q.x,
            scalar * q.y,
            scalar * q.z
        );
    }
}

#endif // PHYTH_QUATERNION_HPP