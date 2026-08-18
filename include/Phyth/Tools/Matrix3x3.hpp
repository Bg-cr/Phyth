#ifndef PHYTH_MATRIX3X3_HPP
#define PHYTH_MATRIX3X3_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Vector3.hpp"
#include <ostream>
#include <stdexcept>

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

        /**
         * @brief Construct a zero matrix
         *
         * Initializes all elements to zero using Quantity's default constructor.
         */
        constexpr Matrix3x3() : x_basis(0, 0, 0), y_basis(0, 0, 0), z_basis(0, 0, 0) {
        }

        /**
         * @brief Construct from three column vectors
         *
         * @param col0 First column vector
         * @param col1 Second column vector
         * @param col2 Third column vector
         */
        constexpr Matrix3x3(const Vector3<T> &col0, const Vector3<T> &col1, const Vector3<T> &col2)
            : x_basis(col0), y_basis(col1), z_basis(col2) {
        }

        /**
         * @brief Construct from individual components (row-major ordering)
         *
         * This constructor takes components in row-major order but stores them
         * internally in column-major order.
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
         *
         * Example:
         *   Matrix3x3 m(1, 2, 3,
         *               4, 5, 6,
         *               7, 8, 9);
         *   // m(1, 0) == 4, m(0, 1) == 2
         */
        constexpr Matrix3x3(
            T m00, T m01, T m02,
            T m10, T m11, T m12,
            T m20, T m21, T m22
        ) : x_basis(m00, m10, m20),
            y_basis(m01, m11, m21),
            z_basis(m02, m12, m22) {
        }

        /**
         * @brief Create an identity matrix
         *
         * Returns a 3x3 identity matrix with ones on the diagonal and zeros elsewhere.
         *
         * @return Matrix3x3 Identity matrix
         *
         * Example:
         *   auto I = Matrix3x3<Quantity<Meter>>::Identity();
         */
        static constexpr Matrix3x3 Identity() {
            return Matrix3x3(
                T(1), T(0), T(0),
                T(0), T(1), T(0),
                T(0), T(0), T(1)
            );
        }

        /**
         * @brief Compute the transpose of the matrix
         *
         * @return Matrix3x3 Transposed matrix
         *
         * Example:
         *   auto m = Matrix3x3(...);
         *   auto mT = m.Transpose();
         */
        [[nodiscard]] constexpr Matrix3x3 Transpose() const {
            return Matrix3x3(
                x_basis.x, y_basis.x, z_basis.x,
                x_basis.y, y_basis.y, z_basis.y,
                x_basis.z, y_basis.z, z_basis.z
            );
        }

        /**
         * @brief Compute the trace of the matrix
         *
         * The trace is the sum of the diagonal elements: M[0][0] + M[1][1] + M[2][2]
         *
         * @return T Sum of diagonal elements
         *
         * Example:
         *   auto tr = m.Trace();  // m00 + m11 + m22
         */
        [[nodiscard]] constexpr T Trace() const {
            return x_basis.x + y_basis.y + z_basis.z;
        }

        /**
         * @brief Compute the symmetric part of the matrix
         *
         * Returns (M + M^T) / 2
         *
         * @return Matrix3x3 Symmetric part
         *
         * Example:
         *   auto sym = m.Symmetric();  // (m + m^T) / 2
         */
        [[nodiscard]] constexpr auto Symmetric() const {
            return (*this + this->Transpose()) * T(0.5);
        }

        /**
         * @brief Compute the antisymmetric part of the matrix
         *
         * Returns (M - M^T) / 2
         *
         * @return Matrix3x3 Antisymmetric part
         *
         * Example:
         *   auto anti = m.Antisymmetric();  // (m - m^T) / 2
         */
        [[nodiscard]] constexpr auto Antisymmetric() const {
            return (*this - this->Transpose()) * T(0.5);
        }

        /**
         * @brief Access element at row i, column j with bounds checking (0-based indexing)
         *
         * @param i Row index (0, 1, 2)
         * @param j Column index (0, 1, 2)
         * @return const T& Element value
         * @throws std::out_of_range if index is out of bounds
         *
         * Example:
         *   auto val = m.at(1, 2);  // Get element at row 1, column 2
         */
        [[nodiscard]] const T& at(const int i, const int j) const {
            if (i < 0 || i > 2 || j < 0 || j > 2) {
                throw std::out_of_range("Matrix3x3::at: index out of range");
            }
            switch (j) {
                case 0: return x_basis[i];
                case 1: return y_basis[i];
                case 2: return z_basis[i];
                default: throw std::out_of_range("Matrix3x3::at: invalid column");
            }
        }

        /**
         * @brief Access element at row i, column j with bounds checking (0-based indexing) - mutable version
         *
         * @param i Row index (0, 1, 2)
         * @param j Column index (0, 1, 2)
         * @return T& Element value
         * @throws std::out_of_range if index is out of bounds
         *
         * Example:
         *   m.at(0, 1) = 5;  // Set element at row 0, column 1 to 5
         */
        T& at(const int i, const int j) {
            if (i < 0 || i > 2 || j < 0 || j > 2) {
                throw std::out_of_range("Matrix3x3::at: index out of range");
            }
            switch (j) {
                case 0: return x_basis[i];
                case 1: return y_basis[i];
                case 2: return z_basis[i];
                default: throw std::out_of_range("Matrix3x3::at: invalid column");
            }
        }

        /**
         * @brief Access element at row i, column j (0-based indexing) - no bounds checking
         *
         * @param i Row index (0, 1, 2)
         * @param j Column index (0, 1, 2)
         * @return const T& Element value
         *
         * @note Matrix is stored in column-major order internally
         * @warning No bounds checking is performed. Use at() for safe access.
         *
         * Example:
         *   auto val = m(1, 2);  // Get element at row 1, column 2
         */
        [[nodiscard]] constexpr const T& operator()(const int i, const int j) const {
            if (j == 0) return x_basis[i];
            if (j == 1) return y_basis[i];
            return z_basis[i];
        }

        /**
         * @brief Access element at row i, column j (0-based indexing) - mutable version
         *
         * @param i Row index (0, 1, 2)
         * @param j Column index (0, 1, 2)
         * @return T& Element value
         *
         * @note Matrix is stored in column-major order internally
         * @warning No bounds checking is performed. Use at() for safe access.
         *
         * Example:
         *   m(0, 1) = 5;  // Set element at row 0, column 1 to 5
         */
        constexpr T& operator()(const int i, const int j) {
            if (j == 0) return x_basis[i];
            if (j == 1) return y_basis[i];
            return z_basis[i];
        }

        /**
         * @brief Matrix addition
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to add
         * @return Matrix3x3<decltype(x_basis + other.x_basis)> Sum of the two matrices
         *
         * Example:
         *   auto sum = m1 + m2;
         */
        template<typename U>
        constexpr auto operator+(const Matrix3x3<U> &other) const {
            return Matrix3x3<decltype(x_basis + other.x_basis)>(
                x_basis + other.x_basis,
                y_basis + other.y_basis,
                z_basis + other.z_basis
            );
        }

        /**
         * @brief Matrix subtraction
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to subtract
         * @return Matrix3x3<decltype(x_basis - other.x_basis)> Difference of the two matrices
         *
         * Example:
         *   auto diff = m1 - m2;
         */
        template<typename U>
        constexpr auto operator-(const Matrix3x3<U> &other) const {
            return Matrix3x3<decltype(x_basis - other.x_basis)>(
                x_basis - other.x_basis,
                y_basis - other.y_basis,
                z_basis - other.z_basis
            );
        }

        /**
         * @brief Unary negation
         *
         * @return Matrix3x3 Negated matrix
         *
         * Example:
         *   auto neg = -m;
         */
        constexpr auto operator-() const {
            return Matrix3x3(-x_basis, -y_basis, -z_basis);
        }

        /**
         * @brief Scalar multiplication (matrix * scalar)
         *
         * @tparam U Scalar type
         * @param scalar Scalar value to multiply with
         * @return Matrix3x3<decltype(x_basis * scalar)> Matrix scaled by scalar
         *
         * Example:
         *   auto scaled = m * 2.0;
         */
        template<typename U>
        constexpr auto operator*(U scalar) const {
            return Matrix3x3<decltype(x_basis * scalar)>(
                x_basis * scalar,
                y_basis * scalar,
                z_basis * scalar
            );
        }

        /**
         * @brief Matrix-vector multiplication: M * v
         *
         * @tparam U Quantity type of the vector
         * @param v Vector to multiply with
         * @return auto Resulting vector
         *
         * Example:
         *   auto v = m * Vector3<Quantity<Meter>>(1, 2, 3);
         */
        template<typename U>
        constexpr auto operator*(const Vector3<U> &v) const {
            // M * v = x_basis * v.x + y_basis * v.y + z_basis * v.z
            return x_basis * v.x + y_basis * v.y + z_basis * v.z;
        }

        /**
         * @brief Matrix-matrix multiplication: A * B
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to multiply with
         * @return Matrix3x3<decltype(*this * other.x_basis)> Product of the two matrices
         *
         * Example:
         *   auto product = m1 * m2;
         */
        template<typename U>
        constexpr auto operator*(const Matrix3x3<U> &other) const {
            // A * B: column j of result = A * (column j of B)
            return Matrix3x3<decltype(*this * other.x_basis)>(
                *this * other.x_basis,
                *this * other.y_basis,
                *this * other.z_basis
            );
        }

        /**
         * @brief Addition assignment
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to add
         * @return Matrix3x3& Reference to this matrix
         *
         * Example:
         *   m += other;
         */
        template<typename U>
        constexpr Matrix3x3 &operator+=(const Matrix3x3<U> &other) {
            x_basis += other.x_basis;
            y_basis += other.y_basis;
            z_basis += other.z_basis;
            return *this;
        }

        /**
         * @brief Subtraction assignment
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to subtract
         * @return Matrix3x3& Reference to this matrix
         *
         * Example:
         *   m -= other;
         */
        template<typename U>
        constexpr Matrix3x3 &operator-=(const Matrix3x3<U> &other) {
            x_basis -= other.x_basis;
            y_basis -= other.y_basis;
            z_basis -= other.z_basis;
            return *this;
        }

        /**
         * @brief Scalar multiplication assignment
         *
         * @tparam U Scalar type
         * @param scalar Scalar value to multiply with
         * @return Matrix3x3& Reference to this matrix
         *
         * Example:
         *   m *= 2.0;
         */
        template<typename U>
        constexpr Matrix3x3 &operator*=(U scalar) {
            x_basis *= scalar;
            y_basis *= scalar;
            z_basis *= scalar;
            return *this;
        }

        /**
         * @brief Equality comparison
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to compare against
         * @return true if all corresponding components are equal
         *
         * Example:
         *   if (m1 == m2) { ... }
         */
        template<typename U>
        constexpr bool operator==(const Matrix3x3<U> &other) const {
            return x_basis == other.x_basis && y_basis == other.y_basis && z_basis == other.z_basis;
        }

        /**
         * @brief Inequality comparison
         *
         * @tparam U Quantity type of the other matrix
         * @param other Matrix to compare against
         * @return true if any corresponding component differs
         *
         * Example:
         *   if (m1 != m2) { ... }
         */
        template<typename U>
        constexpr bool operator!=(const Matrix3x3<U> &other) const {
            return !(*this == other);
        }

        /**
         * @brief Convert all components to a different unit of the same dimension
         *
         * @tparam UnitT Target unit type
         * @return Matrix3x3<Quantity<UnitT>> Matrix with components converted to the new unit
         *
         * Example:
         *   auto m_mm = m.as<Millimeter>();
         */
        template<typename UnitT>
        [[nodiscard]] Matrix3x3<Quantity<UnitT> > as() const {
            return Matrix3x3<Quantity<UnitT> >(
                x_basis.template as<UnitT>(),
                y_basis.template as<UnitT>(),
                z_basis.template as<UnitT>()
            );
        }

        /**
         * @brief Stream insertion operator
         *
         * @param os Output stream
         * @param m Matrix to output
         * @return std::ostream& Reference to the output stream
         *
         * Example:
         *   std::cout << m << std::endl;
         */
        friend std::ostream &operator<<(std::ostream &os, const Matrix3x3 &m) {
            os << "[" << m.x_basis << ",\n"
                    << " " << m.y_basis << ",\n"
                    << " " << m.z_basis << "]";
            return os;
        }
    };

    /**
     * @brief Scalar * Matrix multiplication
     *
     * @tparam T Scalar type
     * @tparam U Quantity type of the matrix
     * @param scalar Scalar value
     * @param m Matrix to multiply
     * @return auto Matrix scaled by scalar
     *
     * Example:
     *   auto scaled = 2.0 * m;
     */
    template<typename T, typename U>
    constexpr auto operator*(T scalar, const Matrix3x3<U> &m) {
        return m * scalar;
    }
}

#endif // PHYTH_MATRIX3X3_HPP