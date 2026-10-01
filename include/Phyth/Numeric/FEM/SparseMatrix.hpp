#ifndef PHYTH_SPARSE_MATRIX_HPP
#define PHYTH_SPARSE_MATRIX_HPP

#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Tools/VectorX.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace Phyth {
    template<typename T>
    class SparseMatrix {
        static_assert(is_quantity_v<T>, "SparseMatrix only supports Quantity types");

    public:
        SparseMatrix(std::size_t rows, std::size_t cols)
            : rows_(rows), cols_(cols) {}

        void Add(std::size_t row, std::size_t col, const T& value) {
            if (row >= rows_ || col >= cols_)
                throw std::out_of_range("SparseMatrix::Add");
            triplets_.push_back({row, col, value});
        }

        void Compress() {
            row_ptr_.assign(rows_ + 1, 0);
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

            std::vector<std::size_t> new_row_ptr(rows_ + 1, 0);
            std::vector<std::size_t> new_col_idx;
            std::vector<T> new_values;
            new_col_idx.reserve(col_idx_.size());
            new_values.reserve(values_.size());

            for (std::size_t i = 0; i < rows_; ++i) {
                const std::size_t begin = row_ptr_[i];
                const std::size_t end = row_ptr_[i + 1];

                if (begin == end) {
                    new_row_ptr[i + 1] = new_row_ptr[i];
                    continue;
                }

                std::vector<std::size_t> idx(end - begin);
                std::iota(idx.begin(), idx.end(), begin);
                std::sort(idx.begin(), idx.end(),
                    [&](const std::size_t a, const std::size_t b) {
                        return col_idx_[a] < col_idx_[b];
                    });

                std::size_t prev_col = static_cast<std::size_t>(-1);
                T acc {0};
                bool has_prev = false;

                for (std::size_t k : idx) {
                    if (has_prev && col_idx_[k] == prev_col) {
                        acc += values_[k];
                    } else {
                        if (has_prev) {
                            new_col_idx.push_back(prev_col);
                            new_values.push_back(acc);
                        }
                        prev_col = col_idx_[k];
                        acc = values_[k];
                        has_prev = true;
                    }
                }
                if (has_prev) {
                    new_col_idx.push_back(prev_col);
                    new_values.push_back(acc);
                }

                new_row_ptr[i + 1] = new_col_idx.size();
            }

            row_ptr_ = std::move(new_row_ptr);
            col_idx_ = std::move(new_col_idx);
            values_ = std::move(new_values);
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
        [[nodiscard]] std::size_t NonZeros() const { return values_.size(); }

        [[nodiscard]] const std::vector<std::size_t>& RowPtr() const { return row_ptr_; }
        [[nodiscard]] const std::vector<std::size_t>& ColIdx() const { return col_idx_; }
        [[nodiscard]] const std::vector<T>& Values() const { return values_; }

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