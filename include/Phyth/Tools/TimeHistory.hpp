#ifndef PHYTH_TIME_HISTORY_HPP
#define PHYTH_TIME_HISTORY_HPP

#include "Phyth/Core/Units.hpp"
#include "QuantityFuncs.hpp"
#include "Phyth/Physical/PhysicalConfig.hpp"

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
            const Scalar tick = time / Config::dt;

            const auto index = Utils::floor(tick).to<long long>();

            if (index < 0 || index >= history_.size()) {
                throw std::out_of_range("Time out of range");
            }

            ValueType result = history_[index];

            const Scalar frac = tick - Scalar(static_cast<double>(index));

            if (frac > Config::epsilon) {
                if (static_cast<size_t>(index + 1) < history_.size()) {
                    result += (history_[index + 1] - result) * frac;
                }
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

            const Quantity<Second> time = static_cast<double>(history_.size()) * Config::dt - offset;
            if (time < 0_s) {
                throw std::out_of_range("Offset exceeds history range");
            }

            return GetValueByTime(time);
        }

        [[nodiscard]] size_t GetSize() const noexcept { return history_.size(); }
        [[nodiscard]] bool IsEmpty() const noexcept { return history_.empty(); }

    protected:
        std::vector<ValueType> history_;
    };
}

#endif //PHYTH_TIME_HISTORY_HPP