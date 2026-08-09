#ifndef PHYTH_VECTOR3_HPP
#define PHYTH_VECTOR3_HPP

#include "Phyth/Core/Dimension.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "QuantityFuncs.hpp"

#include <ostream>

namespace Phyth {
    /**
     * @brief 3D vector with Quantity components
     *
     * Vector3 is a homogeneous 3D vector where all components share the same
     * Quantity type. This ensures that vector operations preserve physical
     * dimensions correctly.
     *
     * @tparam T Component type (must be a Quantity instantiation)
     *
     * Example:
     *   Vector3<Quantity<Meter>> position(1.0_m, 2.0_m, 3.0_m);
     *   Vector3<Quantity<MeterPerSecond>> velocity(1.0_m/s, 0.0_m/s, 0.0_m/s);
     *   auto displacement = velocity * 2.0_s;  // Vector3<Quantity<Meter>>
     */
    template<typename T>
    struct Vector3 {
        static_assert(is_quantity_v<T>, "Vector3 only supports Quantity types");

        /** X component */
        T x;
        /** Y component */
        T y;
        /** Z component */
        T z;

        /** @brief Construct a zero vector (0, 0, 0) */
        constexpr Vector3() : x(0), y(0), z(0) {
        }

        /**
         * @brief Construct a vector with specified components
         *
         * @param x_ X component
         * @param y_ Y component
         * @param z_ Z component
         */
        constexpr Vector3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {
        }

        /**
         * @brief Construct a 2D polar vector (z = 0)
         *
         * @param r Radial distance (must be a Quantity)
         * @param theta_rad Angle from the positive X axis, in radians (dimensionless)
         * @return Vector3 with (r*cos(theta), r*sin(theta), 0)
         *
         * Example:
         *   auto pos = Vector3<Quantity<Meter>>::FromPolar(5.0_m, 30.0_deg);
         *   // pos ~ (4.33_m, 2.5_m, 0_m)
         */
        static constexpr Vector3 FromPolar(T r, const Scalar theta_rad) {
            return Vector3(
                r * Utils::cos(theta_rad),
                r * Utils::sin(theta_rad),
                T(0)
            );
        }

        /**
         * @brief Construct a 3D spherical vector
         *
         * @param r Radial distance (must be a Quantity)
         * @param theta_rad Polar angle from Z axis, in radians (dimensionless)
         * @param phi_rad Azimuthal angle in XY plane from X axis, in radians
         * @return Vector3 with spherical coordinates converted to Cartesian
         *
         * Example:
         *   auto pos = Vector3<Quantity<Meter>>::FromSpherical(
         *       5.0_m, 60.0_deg, 45.0_deg
         *   );
         */
        static constexpr Vector3 FromSpherical(T r, const Scalar theta_rad, const Scalar phi_rad) {
            return Vector3(
                r * Utils::sin(theta_rad) * Utils::cos(phi_rad),
                r * Utils::sin(theta_rad) * Utils::sin(phi_rad),
                r * Utils::cos(theta_rad)
            );
        }

        /**
         * @brief Convert to 2D polar coordinates (z component is ignored)
         *
         * @return std::pair containing (r, theta) where:
         *         - r is the length in the XY plane (Quantity)
         *         - theta is the angle from X axis (Scalar, dimensionless)
         *
         * Note: This function name might be misleading. It returns XY-plane polar
         *       coordinates, not a full 3D conversion. Consider using ToSpherical()
         *       for full 3D conversion.
         */
        [[nodiscard]] auto ToPolar() const {
            Vector3 r = Length();
            auto theta = Utils::atan2(y / r, x / r);
            return std::pair<decltype(r), Scalar>{r, theta};
        }

        /**
         * @brief Convert to spherical coordinates
         *
         * @return std::tuple containing (r, theta, phi) where:
         *         - r is the Euclidean length (Quantity)
         *         - theta is the polar angle from Z axis (Scalar, dimensionless)
         *         - phi is the azimuthal angle from X axis (Scalar, dimensionless)
         */
        [[nodiscard]] auto ToSpherical() const {
            auto r = Length();
            auto theta = Utils::acos(z / r);
            auto phi = Utils::atan2(y / r, x / r);
            return std::tuple<decltype(r), Scalar, Scalar>{r, theta, phi};
        }

        /**
         * @brief Index-based access (non-const)
         *
         * @param i Index: 0 = x, 1 = y, 2 = z
         * @return Reference to the component
         *
         * Example:
         *   vec[0] = 1.0_m;  // modifies x
         */
        constexpr T &operator[](int i) {
            return (&x)[i];
        }

        /**
         * @brief Index-based access (const)
         *
         * @param i Index: 0 = x, 1 = y, 2 = z
         * @return Const reference to the component
         */
        constexpr const T &operator[](int i) const {
            return (&x)[i];
        }

        /**
         * @brief Vector addition
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to add
         * @return Vector3 with component type T + U (promoted if needed)
         */
        template<typename U>
        constexpr auto operator+(const Vector3<U> &other) const {
            return Vector3<decltype(x + other.x)>(
                x + other.x,
                y + other.y,
                z + other.z
            );
        }

        /**
         * @brief Vector subtraction
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to subtract
         * @return Vector3 with component type T - U (promoted if needed)
         */
        template<typename U>
        constexpr auto operator-(const Vector3<U> &other) const {
            return Vector3<decltype(x - other.x)>(
                x - other.x,
                y - other.y,
                z - other.z
            );
        }

        /**
         * @brief Unary negation
         *
         * @return Vector3 with all components negated
         */
        constexpr auto operator-() const {
            return Vector3(-x, -y, -z);
        }

        /**
         * @brief Scalar multiplication (vector * scalar)
         *
         * @tparam U Scalar type (must not be a Vector3)
         * @param scalar The scalar multiplier
         * @return Vector3 with each component multiplied by scalar
         *
         * Note: This is intentionally constrained to non-Vector3 types to avoid
         *       ambiguity with element-wise multiplication.
         */
        template<typename U, typename = std::enable_if_t<!std::is_same_v<std::decay_t<U>, Vector3>> >
        constexpr auto operator*(U scalar) const {
            return Vector3<decltype(x * scalar)>(
                x * scalar,
                y * scalar,
                z * scalar
            );
        }

        /**
         * @brief Scalar division (vector / scalar)
         *
         * @tparam U Scalar type
         * @param scalar The scalar divisor
         * @return Vector3 with each component divided by scalar
         */
        template<typename U>
        constexpr auto operator/(U scalar) const {
            return Vector3<decltype(x / scalar)>(
                x / scalar,
                y / scalar,
                z / scalar
            );
        }

        /**
         * @brief Addition assignment
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to add
         * @return Reference to this vector
         */
        template<typename U>
        constexpr Vector3 &operator+=(const Vector3<U> &other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        /**
         * @brief Subtraction assignment
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to subtract
         * @return Reference to this vector
         */
        template<typename U>
        constexpr Vector3 &operator-=(const Vector3<U> &other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        /**
         * @brief Scalar multiplication assignment
         *
         * @tparam U Scalar type
         * @param scalar The scalar multiplier
         * @return Reference to this vector
         */
        template<typename U>
        constexpr Vector3 &operator*=(U scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        /**
         * @brief Scalar division assignment
         *
         * @tparam U Scalar type
         * @param scalar The scalar divisor
         * @return Reference to this vector
         */
        template<typename U>
        constexpr Vector3 &operator/=(U scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        /**
         * @brief Equality comparison
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to compare against
         * @return true if all components are equal
         */
        template<typename U>
        constexpr bool operator==(const Vector3<U> &other) const {
            return x == other.x && y == other.y && z == other.z;
        }

        /**
         * @brief Inequality comparison
         *
         * @tparam U Quantity type of the other vector
         * @param other Vector to compare against
         * @return true if any component differs
         */
        template<typename U>
        constexpr bool operator!=(const Vector3<U> &other) const {
            return !(*this == other);
        }

        /**
         * @brief Dot product
         *
         * @tparam U Quantity type of the other vector
         * @param other The other vector
         * @return The dot product as a Quantity
         *
         * Example:
         *   auto a = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
         *   auto b = Vector3<Quantity<Meter>>(0_m, 2.0_m, 0_m);
         *   auto dot = a.Dot(b);  // (0 m^2, 0 m^2, 0 m^2) (Meter*Meter = Area)
         */
        template<typename U>
        constexpr auto Dot(const Vector3<U> &other) const -> decltype(x * other.x) {
            return x * other.x + y * other.y + z * other.z;
        }

        /**
         * @brief Cross product
         *
         * @tparam U Quantity type of the other vector
         * @param other The other vector
         * @return Vector3 containing the cross product
         *
         * Example:
         *   auto a = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
         *   auto b = Vector3<Quantity<Meter>>(0_m, 1.0_m, 0_m);
         *   auto cross = a.Cross(b);  // (0 m^2, 0 m^2, 1 m^2)
         */
        template<typename U>
        constexpr auto Cross(const Vector3<U> &other) const {
            return Vector3<decltype(x * other.x)>(
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            );
        }

        /**
         * @brief Euclidean length (magnitude)
         *
         * @return sqrt(x^2 + y^2 + z^2) as a Quantity
         */
        [[nodiscard]] auto Length() const {
            return Utils::sqrt(x * x + y * y + z * z);
        }

        /**
         * @brief Squared Euclidean length
         *
         * @return x^2 + y^2 + z^2 as a Quantity
         *
         * Use this instead of Length() when you only need to compare magnitudes,
         * as it avoids the square root computation.
         */
        [[nodiscard]] auto LengthSquared() const {
            return x * x + y * y + z * z;
        }

        /**
         * @brief Normalized vector (unit vector in the same direction)
         *
         * @return Vector3 with magnitude 1
         *
         * @note The component type becomes Quantity<Unit<DimDivT<DimT, DimT>>>
         *       which is dimensionless. This is intentional.
         */
        [[nodiscard]] Vector3<Scalar> Normalized() const {
            auto len = Length();
            return Vector3<decltype(x / len)>(
                x / len,
                y / len,
                z / len
            );
        }

        /**
         * @brief Project this vector onto another vector
         *
         * @tparam U Quantity type of the target vector
         * @param other The vector to project onto
         * @return The projection of this onto other
         *
         * Example:
         *   auto a = Vector3<Quantity<Meter>>(1.0_m, 1.0_m, 0_m);
         *   auto b = Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
         *   auto proj = a.Project(b);  // (1.0 m, 0 m, 0 m)
         */
        template<typename U>
        [[nodiscard]] auto Project(const Vector3<U> &other) const {
            return Dot(other) / other.LengthSquared() * other;
        }

        /**
         * @brief Convert all components to a different unit of the same dimension
         *
         * @tparam UnitT The target unit type
         * @return Vector3 with each component converted to UnitT
         *
         * Example:
         *   auto pos = Vector3<Quantity<Kilometer>>(1.0_km, 0_km, 0_km);
         *   auto in_meters = pos.as<Meter>();  // (1000 m, 0 m, 0 m)
         */
        template<typename UnitT>
        [[nodiscard]] Vector3<Quantity<UnitT> > as() const {
            return {x.template as<UnitT>(), y.template as<UnitT>(), z.template as<UnitT>()};
        }

        /**
         * @brief Stream insertion operator for Vector3
         *
         * Output format: "(x, y, z)" where each component uses Quantity's
         * stream formatting.
         *
         * Example output: "(1 m, 2 m, 3 m)"
         */
        friend std::ostream &operator<<(std::ostream &os, const Vector3 &v) {
            os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
            return os;
        }
    };

    /**
     * @brief Trait: check if a type is a Vector3 instantiation
     */
    template<typename>
    struct is_vector3 : std::false_type {
    };

    template<typename QuantityT>
    struct is_vector3<Vector3<QuantityT> > : std::true_type {
    };

    /**
     * @brief Convenience variable template for is_vector3
     */
    template<typename VecT>
    inline constexpr auto is_vector3_v = is_vector3<VecT>::value;

    /**
     * @brief Scalar multiplication (scalar * vector)
     *
     * This is the symmetric counterpart to Vector3::operator*(U scalar).
     *
     * @tparam T Scalar type (must not be a Vector3)
     * @tparam U Quantity type of the vector components
     * @param scalar The scalar multiplier
     * @param v The vector
     * @return Vector3 with each component multiplied by scalar
     *
     * Example:
     *   auto v = 2.0 * Vector3<Quantity<Meter>>(1.0_m, 2.0_m, 3.0_m);
     *   // v = (2 m, 4 m, 6 m)
     */
    template<typename T, typename U,
        typename = std::enable_if_t<!is_vector3_v<T>> >
    constexpr auto operator*(T scalar, const Vector3<U> &v) {
        return Vector3<decltype(scalar * v.x)>(
            scalar * v.x,
            scalar * v.y,
            scalar * v.z
        );
    }
}

#endif // PHYTH_VECTOR3_HPP
