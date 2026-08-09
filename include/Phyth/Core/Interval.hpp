#ifndef PHYTH_INTERVAL_HPP
#define PHYTH_INTERVAL_HPP

#include <stdexcept>

#include "Quantity.hpp"
#include "Unit.hpp"

namespace Phyth {
    /**
     * @brief Closed interval [min, max] of physical quantities
     *
     * @tparam UnitT The unit type of the bounds
     *
     * Example:
     *   Interval<Meter> temp(200_m, 400_m);
     *   if (temp.Contains(300_m)) { ... }
     */
    template<typename UnitT>
    struct Interval {
        static_assert(is_unit_v<UnitT>, "Interval requires a Unit type as template parameter");

        /**
         * @brief Construct a zero-length interval containing a single value
         *
         * @param value The only value in the interval
         */
        Interval(Quantity<UnitT> value) noexcept
            : min_(value), max_(value) {
        }

        /**
         * @brief Construct a closed interval [min, max]
         *
         * @param min Lower bound (inclusive)
         * @param max Upper bound (inclusive)
         * @throws std::invalid_argument if min > max
         */
        Interval(Quantity<UnitT> min, Quantity<UnitT> max)
            : min_(min), max_(max) {
            if (min_ > max_) {
                throw std::invalid_argument("Interval: min must be less than max");
            }
        }

        /**
         * @brief Check if a value lies within the interval
         *
         * @param x The value to test
         * @return true if min <= x <= max, false otherwise
         */
        bool Contains(Quantity<UnitT> x) const noexcept {
            return min_ <= x && x <= max_;
        }

        /**
         * @brief Check if another interval is fully contained in this one
         *
         * @param interval The interval to test
         * @return true if interval.max <= this.max && interval.min >= this.min
         */
        bool Contains(const Interval &interval) const noexcept {
            return interval.max_ <= max_ && interval.min_ >= min_;
        }

        /** @return Lower bound of the interval (inclusive) */
        Quantity<UnitT> GetMinimum() const noexcept { return min_; }

        /** @return Upper bound of the interval (inclusive) */
        Quantity<UnitT> GetMaximum() const noexcept { return max_; }

    private:
        Quantity<UnitT> min_;
        Quantity<UnitT> max_;
    };
}

#endif // PHYTH_INTERVAL_HPP
