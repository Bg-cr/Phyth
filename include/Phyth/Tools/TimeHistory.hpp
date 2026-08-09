#ifndef PHYTH_TIME_HISTORY_HPP
#define PHYTH_TIME_HISTORY_HPP

#include "Phyth/Core/Units.hpp"
#include "QuantityFuncs.hpp"
#include "Phyth/Physical/PhysicalConfig.hpp"

#include <vector>
#include <cmath>

namespace Phyth {
    /**
     * @brief Time-based history with linear interpolation
     *
     * TimeHistory stores a sequence of values sampled at fixed time intervals
     * (Config::dt). It provides methods to retrieve values at arbitrary times
     * using linear interpolation between samples.
     *
     * This is useful for:
     *   - Tracking particle trajectories over time
     *   - Storing force or energy history for analysis
     *   - Implementing delay or memory effects in simulations
     *
     * @tparam T The type of value to store (typically Quantity or Vector3<Quantity>)
     *
     * Example:
     *   TimeHistory<Quantity<Meter>> history;
     *   for (double t = 0; t < 10; t += Config::dt) {
     *       history.Register(particle.GetPosition().x);
     *   }
     *   auto pos_at_2_5s = history.GetValueByTime(2.5_s);
     *   auto pos_1s_ago = history.GetValueByOffset(1.0_s);
     */
    template<typename T>
    class TimeHistory {
    public:
        using ValueType = T;

        /**
         * @brief Append a value to the history
         *
         * @param value The value to store at the current time step
         *
         * Each call should correspond to one simulation time step.
         * The time index is implicit: the first call is t=0, second is t=Config::dt, etc.
         */
        void Register(T value) {
            history_.push_back(value);
        }

        /**
         * @brief Get the value at the specified time (linear interpolation)
         *
         * If the time exactly matches a sampling point, returns the stored value.
         * Otherwise, performs linear interpolation between adjacent samples.
         *
         * @param time The time to query (must be >= 0 and within recorded range)
         * @return The interpolated value at the requested time
         * @throws std::out_of_range If time exceeds the historical data range
         *
         * Example:
         *   // If history has samples at t=0, 0.001, 0.002, ...
         *   auto val = history.GetValueByTime(0.0015_s);  // interpolates between t=0.001 and t=0.002
         */
        [[nodiscard]] ValueType GetValueByTime(const Quantity<Second> time) const {
            const Scalar tick = time / Config::dt;

            const auto index = Utils::floor(tick).to<long long>();

            if (index < 0 || index >= static_cast<long long>(history_.size())) {
                throw std::out_of_range("TimeHistory::GetValueByTime: time out of range");
            }

            ValueType result = history_[index];

            const Scalar frac = tick - Scalar(static_cast<double>(index));

            // If there is a fractional part and we have a next sample, interpolate
            if (frac > Config::epsilon) {
                if (static_cast<size_t>(index + 1) < history_.size()) {
                    result += (history_[index + 1] - result) * frac;
                }
            }
            return result;
        }

        /**
         * @brief Get the value at a certain time offset before the current time
         *
         * This is useful for implementing delays or looking up past states.
         * The "current time" is implicitly the last recorded sample point.
         *
         * @param offset The time to look back (must be >= 0 and not exceed history duration)
         * @return The value at (current_time - offset)
         * @throws std::out_of_range If offset is invalid or exceeds history range
         *
         * Example:
         *   // Get the position 0.5 seconds ago
         *   auto past_pos = history.GetValueByOffset(0.5_s);
         *
         * Note: offset = 0_s returns the most recent recorded value (not interpolated
         *       to the exact current time, but the last stored sample point).
         */
        [[nodiscard]] ValueType GetValueByOffset(const Quantity<Second> offset) const {
            if (history_.empty()) {
                throw std::out_of_range("TimeHistory::GetValueByOffset: history is empty");
            }

            if (offset < 0_s) {
                throw std::out_of_range("TimeHistory::GetValueByOffset: offset cannot be negative");
            }

            const Quantity<Second> time = static_cast<double>(history_.size()) * Config::dt - offset;
            if (time < 0_s) {
                throw std::out_of_range("TimeHistory::GetValueByOffset: offset exceeds history range");
            }

            return GetValueByTime(time);
        }

        /** @brief Get the number of stored samples */
        [[nodiscard]] size_t GetSize() const noexcept { return history_.size(); }

        /** @brief Check if the history is empty */
        [[nodiscard]] bool IsEmpty() const noexcept { return history_.empty(); }

    protected:
        std::vector<ValueType> history_;
    };
}

#endif // PHYTH_TIME_HISTORY_HPP
