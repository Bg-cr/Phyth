#ifndef PHYTH_MATRIX_X_HPP
#define PHYTH_MATRIX_X_HPP

#include "VectorX.hpp"
#include "Phyth/Core/Quantity.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

namespace Phyth {
    template<typename T>
    class MatrixX {
        static_assert(is_quantity_v<T>, "MatrixX only supports Quantity types");

    public:
        MatrixX() : rows_(0), cols_(0) {
        }

        MatrixX(std::size_t rows, std::size_t cols, T value = T(0))
            : rows_(rows), cols_(cols), data_(rows * cols, value) {
        }

        [[nodiscard]] std::size_t Rows() const { return rows_; }
        [[nodiscard]] std::size_t Cols() const { return cols_; }

        void Resize(std::size_t rows, std::size_t cols, T value = T(0)) {
            rows_ = rows;
            cols_ = cols;
            data_.resize(rows * cols, value);
        }

        T &operator()(std::size_t i, std::size_t j) { return data_[i * cols_ + j]; }
        const T &operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }

        T &At(std::size_t i, std::size_t j) {
            if (i >= rows_ || j >= cols_) throw std::out_of_range("MatrixX::At");
            return data_[i * cols_ + j];
        }

        const T &At(std::size_t i, std::size_t j) const {
            if (i >= rows_ || j >= cols_) throw std::out_of_range("MatrixX::At");
            return data_[i * cols_ + j];
        }

        template<typename U>
        auto operator+(const MatrixX<U> &other) const {
            CheckSize(other);
            using R = decltype(data_[0] + other(0, 0));
            MatrixX<R> result(rows_, cols_);
            for (std::size_t k = 0; k < data_.size(); ++k)
                result.Data()[k] = data_[k] + other.Data()[k];
            return result;
        }

        template<typename U>
        auto operator-(const MatrixX<U> &other) const {
            CheckSize(other);
            using R = decltype(data_[0] - other(0, 0));
            MatrixX<R> result(rows_, cols_);
            for (std::size_t k = 0; k < data_.size(); ++k)
                result.Data()[k] = data_[k] - other.Data()[k];
            return result;
        }

        auto operator-() const {
            using R = decltype(-data_[0]);
            MatrixX<R> result(rows_, cols_);
            for (std::size_t k = 0; k < data_.size(); ++k)
                result.Data()[k] = -data_[k];
            return result;
        }

        template<typename U, typename = std::enable_if_t<!is_matrix_x_v<U>> >
        auto operator*(U scalar) const {
            using R = decltype(data_[0] * scalar);
            MatrixX<R> result(rows_, cols_);
            for (std::size_t k = 0; k < data_.size(); ++k)
                result.Data()[k] = data_[k] * scalar;
            return result;
        }

        template<typename U>
        auto operator*(const VectorX<U> &x) const {
            if (x.Size() != cols_)
                throw std::invalid_argument("MatrixX * VectorX size mismatch");
            using R = decltype(data_[0] * x[0]);
            VectorX<R> result(rows_, R(0));
            for (std::size_t i = 0; i < rows_; ++i) {
                R sum = R(0);
                for (std::size_t j = 0; j < cols_; ++j)
                    sum += (*this)(i, j) * x[j];
                result[i] = sum;
            }
            return result;
        }

        template<typename U>
        auto operator*(const MatrixX<U> &other) const {
            if (cols_ != other.Rows())
                throw std::invalid_argument("MatrixX * MatrixX size mismatch");
            using R = decltype(data_[0] * other(0, 0));
            MatrixX<R> result(rows_, other.Cols());
            for (std::size_t i = 0; i < rows_; ++i)
                for (std::size_t j = 0; j < other.Cols(); ++j) {
                    R sum = R(0);
                    for (std::size_t k = 0; k < cols_; ++k)
                        sum += (*this)(i, k) * other(k, j);
                    result(i, j) = sum;
                }
            return result;
        }

        MatrixX &operator+=(const MatrixX &other) {
            CheckSize(other);
            for (std::size_t k = 0; k < data_.size(); ++k)
                data_[k] += other.data_[k];
            return *this;
        }

        MatrixX &operator-=(const MatrixX &other) {
            CheckSize(other);
            for (std::size_t k = 0; k < data_.size(); ++k)
                data_[k] -= other.data_[k];
            return *this;
        }

        MatrixX &operator*=(T scalar) {
            for (auto &v: data_) v *= scalar;
            return *this;
        }

        template<typename U>
        bool operator==(const MatrixX<U> &other) const {
            if (rows_ != other.Rows() || cols_ != other.Cols()) return false;
            for (std::size_t k = 0; k < data_.size(); ++k)
                if (!(data_[k] == other.Data()[k])) return false;
            return true;
        }

        template<typename U>
        bool operator!=(const MatrixX<U> &other) const {
            return !(*this == other);
        }

        [[nodiscard]] MatrixX Transpose() const {
            MatrixX result(cols_, rows_);
            for (std::size_t i = 0; i < rows_; ++i)
                for (std::size_t j = 0; j < cols_; ++j)
                    result(j, i) = (*this)(i, j);
            return result;
        }

        [[nodiscard]] T Trace() const {
            if (rows_ != cols_)
                throw std::invalid_argument("MatrixX::Trace requires square matrix");
            T result = T(0);
            for (std::size_t i = 0; i < rows_; ++i)
                result += (*this)(i, i);
            return result;
        }

        [[nodiscard]] auto Symmetric() const {
            return (*this + this->Transpose()) * T(0.5);
        }

        [[nodiscard]] auto Antisymmetric() const {
            return (*this - this->Transpose()) * T(0.5);
        }

        template<typename UnitT>
        [[nodiscard]] MatrixX<Quantity<UnitT> > As() const {
            MatrixX<Quantity<UnitT> > result(rows_, cols_);
            for (std::size_t i = 0; i < rows_; ++i)
                for (std::size_t j = 0; j < cols_; ++j)
                    result(i, j) = (*this)(i, j).template As<UnitT>();
            return result;
        }

        [[nodiscard]] const std::vector<T> &Data() const { return data_; }
        [[nodiscard]] std::vector<T> &Data() { return data_; }

    private:
        template<typename U>
        void CheckSize(const MatrixX<U> &other) const {
            if (rows_ != other.Rows() || cols_ != other.Cols())
                throw std::invalid_argument("MatrixX size mismatch");
        }

        std::size_t rows_, cols_;
        std::vector<T> data_;
    };

    template<typename T>
    struct is_matrix_x : std::false_type {
    };

    template<typename T>
    struct is_matrix_x<MatrixX<T> > : std::true_type {
    };

    template<typename T>
    inline constexpr bool is_matrix_x_v = is_matrix_x<T>::value;

    template<typename T, typename U, typename = std::enable_if_t<!is_matrix_x_v<T>> >
    auto operator*(T scalar, const MatrixX<U> &m) {
        return m * scalar;
    }
}

#endif // PHYTH_MATRIX_X_HPP
