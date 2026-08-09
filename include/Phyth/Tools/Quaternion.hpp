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
     * @brief Quaternion with Quantity components
     *
     * A quaternion represents rotations in 3D space. It consists of a scalar
     * part (w) and a vector part (x, y, z).
     *
     * This implementation uses Quantity types for all components, ensuring that
     * quaternion operations are dimensionally consistent. For rotation
     * quaternions, all components should be dimensionless (Scalar type).
     *
     * @tparam T Component type (must be a Quantity instantiation, typically Scalar)
     *
     * Example:
     *   using Rotation = Quaternion<Scalar>;
     *
     *   // Rotate a vector by 90 degrees around Z axis
     *   auto q = Rotation::FromAxisAngle(Vector3<Scalar>(0, 0, 1), 90.0_deg);
     *   auto v = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     *   auto rotated = q.Rotated(v);  // (0_m, 1_m, 0_m)
     */
    template<typename T>
    struct Quaternion {
        static_assert(is_quantity_v<T>,
                      "Quaternion only supports Quantity types");

        /** Scalar (real) part */
        T w;
        /** X component of the vector (imaginary) part */
        T x;
        /** Y component of the vector (imaginary) part */
        T y;
        /** Z component of the vector (imaginary) part */
        T z;

        /**
         * @brief Default constructor - identity quaternion (1 + 0i + 0j + 0k)
         */
        constexpr Quaternion() : w(1), x(0), y(0), z(0) {
        }

        /**
         * @brief Construct from scalar and vector components
         *
         * @param w_ Scalar part
         * @param x_ X component of vector part
         * @param y_ Y component of vector part
         * @param z_ Z component of vector part
         */
        constexpr Quaternion(T w_, T x_, T y_, T z_) : w(w_), x(x_), y(y_), z(z_) {
        }

        /**
         * @brief Construct from scalar and vector parts
         *
         * @param w_ Scalar part
         * @param v_ Vector part (x, y, z)
         */
        constexpr Quaternion(T w_, const Vector3<T> &v_) : w(w_), x(v_.x), y(v_.y), z(v_.z) {
        }

        /**
         * @brief Construct from vector part only (scalar part = 0)
         *
         * This creates a pure quaternion (0 + v), useful for representing
         * vectors in quaternion form.
         *
         * @param v_ Vector part
         */
        constexpr explicit Quaternion(const Vector3<T> &v_) : w(0), x(v_.x), y(v_.y), z(v_.z) {
        }

        /**
         * @brief Create a unit quaternion from axis and angle
         *
         * @param axis Unit vector axis of rotation (must be normalized)
         * @param theta Rotation angle in radians (dimensionless)
         * @return Unit quaternion representing the rotation
         *
         * Example:
         *   auto axis = Vector3<Scalar>(0, 0, 1);
         *   auto q = Quaternion<Scalar>::FromAxisAngle(axis, 90.0_deg);
         */
        static Quaternion FromAxisAngle(const Vector3<T> &axis, const Scalar theta) {
            const auto half_theta = theta / 2;
            const auto s = Utils::sin(half_theta);
            const auto c = Utils::cos(half_theta);
            return Quaternion(c, axis * s);
        }

        /**
         * @brief Create a unit quaternion from Euler angles (ZYX convention)
         *
         * Converts yaw-pitch-roll Euler angles to a quaternion.
         * The rotation order is: first roll (X), then pitch (Y), then yaw (Z).
         *
         * @param yaw Rotation around Z axis (radians, dimensionless)
         * @param pitch Rotation around Y axis (radians, dimensionless)
         * @param roll Rotation around X axis (radians, dimensionless)
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

        /** @brief Extract the vector part as a Vector3 */
        [[nodiscard]] constexpr Vector3<T> ToVector() const {
            return Vector3<T>(x, y, z);
        }

        /** @brief Extract the scalar part */
        [[nodiscard]] constexpr T ToQuantity() const {
            return w;
        }

        /**
         * @brief Squared magnitude
         *
         * Returns w^2 + x^2 + y^2 + z^2
         */
        [[nodiscard]] auto LengthSquared() const {
            return w * w + x * x + y * y + z * z;
        }

        /** @brief Magnitude (sqrt of LengthSquared) */
        [[nodiscard]] auto Length() const {
            return Utils::sqrt(LengthSquared());
        }

        /**
         * @brief Conjugate: (w + xi + yj + zk) -> (w - xi - yj - zk)
         *
         * For unit quaternions, the conjugate is also the inverse.
         */
        [[nodiscard]] constexpr Quaternion Conjugated() const {
            return Quaternion(w, -x, -y, -z);
        }

        /**
         * @brief Inverse: q^(-1) = conjugate / |q|^2
         *
         * For unit quaternions, inverse equals conjugate.
         */
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

        /** @brief Normalized quaternion (unit quaternion) */
        [[nodiscard]] auto Normalized() const {
            const auto length = Length();
            return Quaternion<decltype(w / length)>(
                w / length,
                x / length,
                y / length,
                z / length
            );
        }

        /**
         * @brief Check if the quaternion is a unit quaternion
         *
         * Returns true if |w^2 + x^2 + y^2 + z^2 - 1| < Config::tolerance
         */
        [[nodiscard]] bool IsUnit() const {
            return Utils::abs(Scalar(1) - LengthSquared()) < Config::tolerance;
        }

        /**
         * @brief Rotate a vector by this quaternion
         *
         * The rotation is performed as: q * (0 + v) * q^(-1)
         *
         * @param v The vector to rotate
         * @return The rotated vector
         *
         * Example:
         *   auto q = Quaternion<Scalar>::FromAxisAngle(axis, 90.0_deg);
         *   auto v = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
         *   auto rotated = q.Rotated(v);
         */
        [[nodiscard]] auto Rotated(const Vector3<T> &v) const {
            const auto qv = *this * Quaternion(v);
            const auto result = qv * Inversed();
            return result.ToVector();
        }

        /**
         * @brief Quaternion multiplication
         *
         * The product of two quaternions represents the composition of
         * their rotations.
         */
        template<typename U>
        constexpr auto operator*(const Quaternion<U> &other) const {
            using R = decltype(w * other.w);
            return Quaternion<R>(
                w * other.w - x * other.x - y * other.y - z * other.z,
                w * other.x + x * other.w + y * other.z - z * other.y,
                w * other.y - x * other.z + y * other.w + z * other.x,
                w * other.z + x * other.y - y * other.x + z * other.w
            );
        }

        /** @brief Quaternion addition */
        template<typename U>
        constexpr auto operator+(const Quaternion<U> &other) const {
            return Quaternion<decltype(w + other.w)>(
                w + other.w,
                x + other.x,
                y + other.y,
                z + other.z
            );
        }

        /** @brief Quaternion subtraction */
        template<typename U>
        constexpr auto operator-(const Quaternion<U> &other) const {
            return Quaternion<decltype(w - other.w)>(
                w - other.w,
                x - other.x,
                y - other.y,
                z - other.z
            );
        }

        /** @brief Unary negation */
        constexpr auto operator-() const {
            return Quaternion(-w, -x, -y, -z);
        }

        /**
         * @brief Scalar multiplication
         *
         * @param scalar The scalar multiplier (must be a Quantity, typically Scalar)
         * @return Quaternion with all components multiplied by scalar
         */
        template<typename U>
        constexpr auto operator*(Quantity<U> scalar) const {
            return Quaternion<decltype(w * scalar)>(
                w * scalar,
                x * scalar,
                y * scalar,
                z * scalar
            );
        }

        /**
         * @brief Scalar division
         *
         * @param scalar The scalar divisor (must be a Quantity, typically Scalar)
         * @return Quaternion with all components divided by scalar
         */
        template<typename U>
        constexpr auto operator/(U scalar) const {
            return Quaternion<decltype(w / scalar)>(
                w / scalar,
                x / scalar,
                y / scalar,
                z / scalar
            );
        }

        /** @brief Addition assignment */
        template<typename U>
        constexpr Quaternion &operator+=(const Quaternion<U> &other) {
            w += other.w;
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        /** @brief Subtraction assignment */
        template<typename U>
        constexpr Quaternion &operator-=(const Quaternion<U> &other) {
            w -= other.w;
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        /** @brief Scalar multiplication assignment */
        template<typename U>
        constexpr Quaternion &operator*=(U scalar) {
            w *= scalar;
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        /** @brief Scalar division assignment */
        template<typename U>
        constexpr Quaternion &operator/=(U scalar) {
            w /= scalar;
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        /** @brief Equality comparison */
        template<typename U>
        constexpr bool operator==(const Quaternion<U> &other) const {
            return w == other.w && x == other.x && y == other.y && z == other.z;
        }

        /** @brief Inequality comparison */
        template<typename U>
        constexpr bool operator!=(const Quaternion<U> &other) const {
            return !(*this == other);
        }

        /**
         * @brief Convert all components to a different unit of the same dimension
         *
         * @tparam UnitT The target unit type
         * @return Quaternion with each component converted to UnitT
         */
        template<typename UnitT>
        [[nodiscard]] Quaternion<Quantity<UnitT> > as() const {
            return {
                w.template as<UnitT>(),
                x.template as<UnitT>(),
                y.template as<UnitT>(),
                z.template as<UnitT>()
            };
        }

        /**
         * @brief Spherical linear interpolation between two quaternions
         *
         * Interpolates along the shortest path on the unit sphere.
         * Handles negative dot products to maintain the shortest path.
         *
         * @tparam U Component type of the second quaternion
         * @param q1 Start quaternion
         * @param q2 End quaternion
         * @param t Interpolation parameter [0, 1] (dimensionless scalar)
         * @return Interpolated quaternion (unit quaternion)
         *
         * Example:
         *   auto q1 = Quaternion<Scalar>();  // identity
         *   auto q2 = Quaternion<Scalar>::FromAxisAngle(axis, 90.0_deg);
         *   auto q_mid = Quaternion<Scalar>::Slerp(q1, q2, 0.5);  // 45-degree rotation
         */
        template<typename U>
        static auto Slerp(const Quaternion &q1, const Quaternion<U> &q2, const Scalar t) {
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

        /**
         * @brief Dot product of two quaternions
         *
         * @tparam U Component type of the other quaternion
         * @param other The other quaternion
         * @return The dot product as a Quantity
         */
        template<typename U>
        [[nodiscard]] constexpr auto Dot(const Quaternion<U> &other) const {
            return w * other.w + x * other.x + y * other.y + z * other.z;
        }

        /**
         * @brief Index-based access (non-const)
         *
         * @param i Index: 0 = w, 1 = x, 2 = y, 3 = z
         * @return Reference to the component
         */
        constexpr T &operator[](int i) {
            return (&w)[i];
        }

        /**
         * @brief Index-based access (const)
         *
         * @param i Index: 0 = w, 1 = x, 2 = y, 3 = z
         * @return Const reference to the component
         */
        constexpr const T &operator[](int i) const {
            return (&w)[i];
        }

        /**
         * @brief Stream insertion operator for Quaternion
         *
         * Output format: "(w + xi + yj + zk)" where each component uses
         * Quantity's stream formatting.
         */
        friend std::ostream &operator<<(std::ostream &os, const Quaternion &q) {
            os << "(" << q.w << " + " << q.x << "i + " << q.y << "j + " << q.z << "k)";
            return os;
        }
    };

    /**
     * @brief Trait: check if a type is a Quaternion instantiation
     */
    template<typename>
    struct is_quaternion : std::false_type {
    };

    template<typename T>
    struct is_quaternion<Quaternion<T> > : std::true_type {
    };

    /**
     * @brief Convenience variable template for is_quaternion
     */
    template<typename T>
    inline constexpr auto is_quaternion_v = is_quaternion<T>::value;

    /**
     * @brief Scalar * Quaternion multiplication (commutative)
     *
     * @tparam T Scalar type (must be a Quantity, typically Scalar)
     * @tparam U Quaternion component type
     * @param scalar The scalar multiplier
     * @param q The quaternion
     * @return Quaternion with all components multiplied by scalar
     */
    template<typename T, typename U>
    constexpr auto operator*(const Quantity<T> scalar, const Quaternion<U> &q) {
        return Quaternion<decltype(scalar * q.w)>(
            scalar * q.w,
            scalar * q.x,
            scalar * q.y,
            scalar * q.z
        );
    }
}

#endif // PHYTH_QUATERNION_HPP
