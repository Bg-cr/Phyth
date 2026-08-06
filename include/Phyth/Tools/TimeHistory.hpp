#ifndef PHYTH_TIME_HISTORY_HPP
#define PHYTH_TIME_HISTORY_HPP

#include "Phyth/Core/Units.hpp"
#include "QuantityFuncs.hpp"

#include <vector>
#include <cmath>

namespace Phyth {
    template<typename T>
    class TimeHistory {
    public:
        using ValueType = T;

        void Register(T value) {
            history_.push_back(value);
        }

        /**
         * @brief Get the value at the specified time (linear interpolation)
         *
         * If time happens to fall on the sampling point, return the value of that point directly;
         * Otherwise, perform linear interpolation between adjacent sampling points
         *
         * @throws std::out_of_range If the time exceeds the historical data range
         */
        [[nodiscard]] ValueType GetValueByTime(const Quantity<Second> time) const {
            if (dt_ <= 0_s) {
                throw std::runtime_error("Invalid dt_ in TimeHistory");
            }

            const Scalar tick = time / dt_;

            const auto index = Utils::floor(tick).to<long long>();

            if (index < 0 || index >= history_.size()) {
                throw std::out_of_range("Time out of range");
            }

            ValueType result = history_[index];
    
            // 计算小数部分
            const Scalar frac = tick - Scalar(static_cast<double>(index));
    
            // 如果需要插值，且下一个点存在
            if (frac > Config::epsilon) {
                if (static_cast<size_t>(index + 1) < history_.size()) {
                    // 正常插值
                    result += (history_[index + 1] - result) * frac;
                }
                // 如果下一个点不存在，直接返回当前点（不插值）
                // 这相当于线性外推被截断为常数
            }
            return result;
        }

        /**
         * @brief Get the value at a certain point in time before the current time
         *
         * @param offset The length of time to push forward
         * @throws std::out_of_range If offset is invalid or exceeds history range
         */
        [[nodiscard]] ValueType GetValueByOffset(const Quantity<Second> offset) const {
            if (history_.empty()) {
                throw std::out_of_range("TimeHistory is empty");
            }

            if (offset < 0_s) {
                throw std::out_of_range("Offset cannot be negative");
            }

            const Quantity<Second> time = static_cast<double>(history_.size()) * dt_ - offset;
            if (time < 0_s) {
                throw std::out_of_range("Offset exceeds history range");
            }

            return GetValueByTime(time);
        }

        void SetDeltaTime(const Quantity<Second> dt) {
            if (dt <= 0_s) {
                throw std::invalid_argument("dt must be positive");
            }
            dt_ = dt;
        }

        [[nodiscard]] Quantity<Second> GetDeltaTime() const { return dt_; }
        [[nodiscard]] size_t GetSize() const noexcept { return history_.size(); }
        [[nodiscard]] bool IsEmpty() const noexcept { return history_.empty(); }

    protected:
        std::vector<ValueType> history_;
        Quantity<Second> dt_;
    };
}

#endif //PHYTH_TIME_HISTORY_HPP