#ifndef PHYTH_TIME_HISTORY_HPP
#define PHYTH_TIME_HISTORY_HPP

#include "Phyth/Core/Units.hpp"
#include "QuantityFuncs.hpp"

#include <vector>

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
            const Scalar tick = time / dt_;
            if (tick < 0 || tick >= history_.size()) {
                throw std::out_of_range("Time out of range");
            }

            ValueType result = history_[tick.to<long long>()];
            if (Utils::fmod(tick, 1_) > Config::epsilon) {
                if (tick >= history_.size() - 1) {
                    throw std::out_of_range("Time out of range");
                }
                result += (history_[tick.to<long long>() + 1] - result) * Utils::fmod(tick, 1_);
            }
            return result;
        }

        /**
         * @brief Get the value of a certain point in time before the current time
         *
         * @param offset The length of time to push forward
         */
        [[nodiscard]] ValueType GetValueByOffset(const Quantity<Second> offset) const {
            return GetValueByTime(history_.size() * dt_ - offset);
        }

        void SetDeltaTime(const Quantity<Second> dt) { dt_ = dt; }
        [[nodiscard]] Quantity<Second> GetDeltaTime() const { return dt_; }

        [[nodiscard]] size_t GetSize() const noexcept { return history_.size(); }

    protected:
        std::vector<ValueType> history_;
        Quantity<Second> dt_;
    };
}

#endif //PHYTH_TIME_HISTORY_HPP
