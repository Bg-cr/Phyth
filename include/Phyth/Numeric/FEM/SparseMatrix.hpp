#ifndef PHYTH_SPARSE_MATRIX_HPP
#define PHYTH_SPARSE_MATRIX_HPP

#include "Phyth/Core/Quantity.hpp"
#include "../../Tools/VectorX.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace Phyth {
    template<typename T>
    class SparseMatrix {
        static_assert(is_quantity_v<T>, "SparseMatrix only supports Quantity types");
    public:
        SparseMatrix(const std::size_t rows, const std::size_t cols)
            : rows_(rows), cols_(cols) {}

        void Add(const std::size_t row, const std::size_t col, const T& value) {
            if (row >= rows_)
                throw std::out_of_range(
                    "SparseMatrix::Add: row (witch is " + std::to_string(row) + ") >= rows_ (witch is " + std::to_string(rows_) + ")"
                );
            
            if (col >= cols_)
                throw std::out_of_range(
                    "SparseMatrix::Add: col (witch is " + std::to_string(col) + ") >= cols_ (witch is " + std::to_string(cols_) + ")"
                );
            triplets_.push_back({row, col, value});
        }

        void Compress() {
            row_ptr_.assign(rows_ + 1, 0);
            col_idx_.clear();
            values_.clear();

            for (const auto& t : triplets_)
                ++row_ptr_[t.row + 1];

            for (std::size_t i = 0; i < rows_; ++i)
                row_ptr_[i + 1] += row_ptr_[i];

            col_idx_.resize(triplets_.size());
            values_.resize(triplets_.size());

            std::vector<std::size_t> cursor = row_ptr_;
            for (const auto& t : triplets_) {
                std::size_t pos = cursor[t.row]++;
                col_idx_[pos] = t.col;
                values_[pos] = t.value;
            }

            triplets_.clear();
            triplets_.shrink_to_fit();
        }

        [[nodiscard]] VectorX<T> Multiply(const VectorX<T>& x) const {
            if (x.Size() != cols_)
                throw std::invalid_argument("SparseMatrix::Multiply size mismatch");

            VectorX<T> result(rows_, T(0));
            for (std::size_t i = 0; i < rows_; ++i) {
                T sum = T(0);
                for (std::size_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k)
                    sum += values_[k] * x[col_idx_[k]];
                result[i] = sum;
            }
            return result;
        }

        [[nodiscard]] std::size_t Rows() const { return rows_; }
        [[nodiscard]] std::size_t Cols() const { return cols_; }

        [[nodiscard]] std::size_t NonZeros() const {
            return triplets_.empty() ? values_.size() : triplets_.size();
        }

    private:
        struct Triplet {
            std::size_t row, col;
            T value;
        };

        std::size_t rows_, cols_;
        std::vector<Triplet> triplets_;
        std::vector<std::size_t> row_ptr_;
        std::vector<std::size_t> col_idx_;
        std::vector<T> values_;
    };
}

#endif // PHYTH_SPARSE_MATRIX_HPP