#ifndef PHYTH_MATRIX3X3_HPP
#define PHYTH_MATRIX3X3_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Vector3.hpp"
#include <ostream>

namespace Phyth {
    /**
     * @brief 3x3 matrix with Quantity components
     *
     * Matrix3x3 stores three column vectors, enabling matrix-vector operations
     * while preserving physical dimensions.
     *
     * @tparam T Component type (must be a Quantity instantiation)
     *
     * Example:
     *   Matrix3x3<Quantity<Meter>> m;  // 3x3 matrix with meter components
     *   auto result = m * Vector3<Quantity<Meter>>(1.0_m, 0_m, 0_m);
     */
    template<typename T>
    struct Matrix3x3 {
        static_assert(is_quantity_v<T>, "Matrix3x3 only supports Quantity types");

        Vector3<T> x_basis;
        Vector3<T> y_basis;
        Vector3<T> z_basis;

        /** @brief Construct a zero matrix */
        constexpr Matrix3x3() : x_basis(0, 0, 0), y_basis(0, 0, 0), z_basis(0, 0, 0) {
        }

        /**
         * @brief Construct from three column vectors
         *
         * @param col0 First column
         * @param col1 Second column
         * @param col2 Third column
         */
        constexpr Matrix3x3(const Vector3<T> &col0, const Vector3<T> &col1, const Vector3<T> &col2)
            : x_basis(col0), y_basis(col1), z_basis(col2) {
        }

        /**
         * @brief Construct from individual components (row-major ordering)
         *
         * @param m00 Row 0, Col 0
         * @param m01 Row 0, Col 1
         * @param m02 Row 0, Col 2
         * @param m10 Row 1, Col 0
         * @param m11 Row 1, Col 1
         * @param m12 Row 1, Col 2
         * @param m20 Row 2, Col 0
         * @param m21 Row 2, Col 1
         * @param m22 Row 2, Col 2
         */
        constexpr Matrix3x3(
            T m00, T m01, T m02,
            T m10, T m11, T m12,
            T m20, T m21, T m22
        ) : x_basis(m00, m10, m20),
            y_basis(m01, m11, m21),
            z_basis(m02, m12, m22) {
        }

        /** @brief Identity matrix */
        static constexpr Matrix3x3 Identity() {
            return Matrix3x3(
                T(1), T(0), T(0),
                T(0), T(1), T(0),
                T(0), T(0), T(1)
            );
        }

        /** @brief Transpose */
        [[nodiscard]] constexpr Matrix3x3 Transpose() const {
            return Matrix3x3(
                x_basis.x, y_basis.x, z_basis.x,
                x_basis.y, y_basis.y, z_basis.y,
                x_basis.z, y_basis.z, z_basis.z
            );
        }

        /** @brief Trace (sum of diagonal elements) */
        [[nodiscard]] constexpr T Trace() const {
            return x_basis.x + y_basis.y + z_basis.z;
        }

        /** @brief Symmetric part: (M + M^T) / 2 */
        [[nodiscard]] constexpr auto Symmetric() const {
            return (*this + this->Transpose()) * T(0.5);
        }

        /** @brief Antisymmetric part: (M - M^T) / 2 */
        [[nodiscard]] constexpr auto Antisymmetric() const {
            return (*this - this->Transpose()) * T(0.5);
        }

        /** @brief Matrix addition */
        template<typename U>
        constexpr auto operator+(const Matrix3x3<U> &other) const {
            return Matrix3x3<decltype(x_basis + other.x_basis)>(
                x_basis + other.x_basis,
                y_basis + other.y_basis,
                z_basis + other.z_basis
            );
        }

        /** @brief Matrix subtraction */
        template<typename U>
        constexpr auto operator-(const Matrix3x3<U> &other) const {
            return Matrix3x3<decltype(x_basis - other.x_basis)>(
                x_basis - other.x_basis,
                y_basis - other.y_basis,
                z_basis - other.z_basis
            );
        }

        /** @brief Unary negation */
        constexpr auto operator-() const {
            return Matrix3x3(-x_basis, -y_basis, -z_basis);
        }

        /** @brief Scalar multiplication (matrix * scalar) */
        template<typename U>
        constexpr auto operator*(U scalar) const {
            return Matrix3x3<decltype(x_basis * scalar)>(
                x_basis * scalar,
                y_basis * scalar,
                z_basis * scalar
            );
        }

        /** @brief Matrix-vector multiplication: M * v */
        template<typename U>
        constexpr auto operator*(const Vector3<U> &v) const {
            // M * v = x_basis * v.x + y_basis * v.y + z_basis * v.z
            return x_basis * v.x + y_basis * v.y + z_basis * v.z;
        }

        /** @brief Matrix-matrix multiplication: A * B */
        template<typename U>
        constexpr auto operator*(const Matrix3x3<U> &other) const {
            // A * B: column j of result = A * (column j of B)
            return Matrix3x3<decltype(*this * other.x_basis)>(
                *this * other.x_basis,
                *this * other.y_basis,
                *this * other.z_basis
            );
        }

        /** @brief Addition assignment */
        template<typename U>
        constexpr Matrix3x3 &operator+=(const Matrix3x3<U> &other) {
            x_basis += other.x_basis;
            y_basis += other.y_basis;
            z_basis += other.z_basis;
            return *this;
        }

        /** @brief Subtraction assignment */
        template<typename U>
        constexpr Matrix3x3 &operator-=(const Matrix3x3<U> &other) {
            x_basis -= other.x_basis;
            y_basis -= other.y_basis;
            z_basis -= other.z_basis;
            return *this;
        }

        /** @brief Scalar multiplication assignment */
        template<typename U>
        constexpr Matrix3x3 &operator*=(U scalar) {
            x_basis *= scalar;
            y_basis *= scalar;
            z_basis *= scalar;
            return *this;
        }

        /** @brief Equality comparison */
        template<typename U>
        constexpr bool operator==(const Matrix3x3<U> &other) const {
            return x_basis == other.x_basis && y_basis == other.y_basis && z_basis == other.z_basis;
        }

        /** @brief Inequality comparison */
        template<typename U>
        constexpr bool operator!=(const Matrix3x3<U> &other) const {
            return !(*this == other);
        }

        /** @brief Convert all components to a different unit of the same dimension */
        template<typename UnitT>
        [[nodiscard]] Matrix3x3<Quantity<UnitT> > as() const {
            return Matrix3x3<Quantity<UnitT> >(
                x_basis.template as<UnitT>(),
                y_basis.template as<UnitT>(),
                z_basis.template as<UnitT>()
            );
        }

        /** @brief Stream insertion operator */
        friend std::ostream &operator<<(std::ostream &os, const Matrix3x3 &m) {
            os << "[" << m.x_basis << ",\n"
                    << " " << m.y_basis << ",\n"
                    << " " << m.z_basis << "]";
            return os;
        }
    };

    /** @brief Scalar * Matrix */
    template<typename T, typename U>
    constexpr auto operator*(T scalar, const Matrix3x3<U> &m) {
        return m * scalar;
    }
}

#endif // PHYTH_MATRIX3X3_HPP
