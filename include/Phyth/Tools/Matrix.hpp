#ifndef PHYTH_MATRIX_HPP
#define PHYTH_MATRIX_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Vector3.hpp"
#include <array>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Phyth {
    template<std::size_t Row, std::size_t Col, typename T>
    struct Matrix {
        static_assert(is_quantity_v<T>, "Matrix only supports Quantity types");
        static_assert(Row > 0 && Col > 0, "Matrix dimensions must be positive");

        std::array<Vector<Row, T>, Col> datas;

        constexpr Matrix()
            : Matrix(MakeZero(std::make_index_sequence<Col>{})) {}

        template<typename... Args,
            std::enable_if_t<
                sizeof...(Args) == Col &&
                std::conjunction_v<std::is_convertible<Args, Vector<Row, T> >...>,
                int> = 0>
        constexpr Matrix(Args... basics) : datas{basics...} {}

        template<typename... Args,
            std::enable_if_t<
                sizeof...(Args) == Row * Col &&
                std::conjunction_v<std::is_convertible<Args, T>...>,
                int> = 0>
        constexpr Matrix(Args... elems)
            : Matrix(std::array<T, Row * Col>{static_cast<T>(elems)...},
                     std::make_index_sequence<Col>{}) {}

        static constexpr Matrix Identity() {
            static_assert(Row == Col, "Identity matrix requires a square matrix");
            return IdentityImpl(std::make_index_sequence<Col>{});
        }

        [[nodiscard]] constexpr auto Transpose() const {
            return TransposeImpl(std::make_index_sequence<Row>{},
                                 std::make_index_sequence<Col>{});
        }

        [[nodiscard]] constexpr T Trace() const {
            static_assert(Row == Col, "Trace requires a square matrix");
            return TraceImpl(std::make_index_sequence<Row>{});
        }

        [[nodiscard]] constexpr auto Symmetric() const {
            return (*this + this->Transpose()) * T(0.5);
        }

        [[nodiscard]] constexpr auto Antisymmetric() const {
            return (*this - this->Transpose()) * T(0.5);
        }

        [[nodiscard]] const T &At(const std::size_t i, const std::size_t j) const {
            if (i >= Row || j >= Col) {
                throw std::out_of_range("Matrix::At: index out of range");
            }
            return datas[j][i];
        }

        T &At(const std::size_t i, const std::size_t j) {
            if (i >= Row || j >= Col) {
                throw std::out_of_range("Matrix::At: index out of range");
            }
            return datas[j][i];
        }

        [[nodiscard]] constexpr const T &operator()(const std::size_t i, const std::size_t j) const {
            return datas[j][i];
        }

        constexpr T &operator()(const std::size_t i, const std::size_t j) {
            return datas[j][i];
        }

        template<typename U>
        constexpr auto operator+(const Matrix<Row, Col, U> &other) const {
            return AddImpl(other, std::make_index_sequence<Col>{});
        }

        template<typename U>
        constexpr auto operator-(const Matrix<Row, Col, U> &other) const {
            return SubImpl(other, std::make_index_sequence<Col>{});
        }

        constexpr auto operator-() const {
            return NegImpl(std::make_index_sequence<Col>{});
        }

        template<typename U>
        constexpr auto operator*(U scalar) const {
            return ScaleImpl(scalar, std::make_index_sequence<Col>{});
        }

        template<typename U>
        constexpr auto operator*(const Vector<Row, U> &v) const {
            using ResultVec = decltype(datas[0] * v[0]);
            ResultVec result = datas[0] * v[0];
            for (std::size_t c = 1; c < Col; ++c) {
                result = result + datas[c] * v[c];
            }
            return result;
        }

        template<typename U, std::size_t OtherCol>
        constexpr auto operator*(const Matrix<Col, OtherCol, U> &other) const {
            using ResultMat = Matrix<Row, OtherCol, decltype(*this * other.datas[0])>;
            ResultMat result;
            for (std::size_t c = 0; c < OtherCol; ++c) {
                result.datas[c] = *this * other.datas[c];
            }
            return result;
        }

        template<typename U>
        constexpr Matrix &operator+=(const Matrix<Row, Col, U> &other) {
            for (std::size_t c = 0; c < Col; ++c) {
                datas[c] += other.datas[c];
            }
            return *this;
        }

        template<typename U>
        constexpr Matrix &operator-=(const Matrix<Row, Col, U> &other) {
            for (std::size_t c = 0; c < Col; ++c) {
                datas[c] -= other.datas[c];
            }
            return *this;
        }

        template<typename U>
        constexpr Matrix &operator*=(U scalar) {
            for (std::size_t c = 0; c < Col; ++c) {
                datas[c] *= scalar;
            }
            return *this;
        }

        template<typename U>
        constexpr bool operator==(const Matrix<Row, Col, U> &other) const {
            for (std::size_t c = 0; c < Col; ++c) {
                if (!(datas[c] == other.datas[c])) return false;
            }
            return true;
        }

        template<typename U>
        constexpr bool operator!=(const Matrix<Row, Col, U> &other) const {
            return !(*this == other);
        }

        template<typename UnitT>
        [[nodiscard]] Matrix<Row, Col, Quantity<UnitT> > As() const {
            return AsImpl<UnitT>(std::make_index_sequence<Col>{});
        }

        friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
            os << "[";
            for (std::size_t r = 0; r < Row; ++r) {
                os << (r == 0 ? "" : "\n ");
                for (std::size_t c = 0; c < Col; ++c) {
                    if (c != 0) os << " ";
                    os << m.datas[c][r];
                }
            }
            os << "]";
            return os;
        }

    private:
        template<std::size_t... Cs>
        static constexpr Matrix MakeZero(std::index_sequence<Cs...>) {
            return Matrix((static_cast<void>(Cs), MakeZeroCol(std::make_index_sequence<Row>{}))...);
        }

        template<std::size_t... Rs>
        static constexpr Vector<Row, T> MakeZeroCol(std::index_sequence<Rs...>) {
            return Vector<Row, T>((static_cast<void>(Rs), T(0))...);
        }

        template<std::size_t... Cs>
        static constexpr Matrix IdentityImpl(std::index_sequence<Cs...>) {
            return Matrix(MakeIdentityCol<Cs>(std::make_index_sequence<Row>{})...);
        }

        template<std::size_t C, std::size_t... Rs>
        static constexpr Vector<Row, T> MakeIdentityCol(std::index_sequence<Rs...>) {
            return Vector<Row, T>((Rs == C ? T(1) : T(0))...);
        }

        template<std::size_t... Cs>
        constexpr Matrix(const std::array<T, Row * Col> &flat,
                         std::index_sequence<Cs...>)
            : datas{MakeCol<Cs>(flat, std::make_index_sequence<Row>{})...} {}

        template<std::size_t C, std::size_t... Rs>
        static constexpr Vector<Row, T>
        MakeCol(const std::array<T, Row * Col> &flat, std::index_sequence<Rs...>) {
            return Vector<Row, T>(flat[Rs * Col + C]...);
        }

        template<std::size_t... Rs, std::size_t... Cs>
        constexpr auto TransposeImpl(std::index_sequence<Rs...>,
                                     std::index_sequence<Cs...>) const {
            return Matrix<Col, Row, T>(
                MakeTransposeCol<Rs>(std::make_index_sequence<Col>{})...);
        }

        template<std::size_t R, std::size_t... Cs>
        constexpr Vector<Col, T> MakeTransposeCol(std::index_sequence<Cs...>) const {
            return Vector<Col, T>(datas[Cs][R]...);
        }

        template<std::size_t... Rs>
        constexpr T TraceImpl(std::index_sequence<Rs...>) const {
            return (datas[Rs][Rs] + ...);
        }

        template<typename U, std::size_t... Cs>
        constexpr auto AddImpl(const Matrix<Row, Col, U> &other,
                               std::index_sequence<Cs...>) const {
            return Matrix<Row, Col, decltype(datas[0] + other.datas[0])>(
                (datas[Cs] + other.datas[Cs])...);
        }

        template<typename U, std::size_t... Cs>
        constexpr auto SubImpl(const Matrix<Row, Col, U> &other,
                               std::index_sequence<Cs...>) const {
            return Matrix<Row, Col, decltype(datas[0] - other.datas[0])>(
                (datas[Cs] - other.datas[Cs])...);
        }

        template<std::size_t... Cs>
        constexpr auto NegImpl(std::index_sequence<Cs...>) const {
            return Matrix<Row, Col, decltype(-datas[0])>((-datas[Cs])...);
        }

        template<typename U, std::size_t... Cs>
        constexpr auto ScaleImpl(U scalar, std::index_sequence<Cs...>) const {
            return Matrix<Row, Col, decltype(datas[0] * scalar)>(
                (datas[Cs] * scalar)...);
        }

        template<typename UnitT, std::size_t... Cs>
        constexpr Matrix<Row, Col, Quantity<UnitT> >
        AsImpl(std::index_sequence<Cs...>) const {
            return Matrix<Row, Col, Quantity<UnitT> >(
                (datas[Cs].template As<UnitT>())...);
        }
    };

    template<typename T, typename U, std::size_t Row, std::size_t Col>
    constexpr auto operator*(T scalar, const Matrix<Row, Col, U> &m) {
        return m * scalar;
    }
}

#endif // PHYTH_MATRIX_HPP