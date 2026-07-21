#ifndef PHYTH_INTERVAL_HPP
#define PHYTH_INTERVAL_HPP

#include <stdexcept>

#include "Quantity.hpp"
#include "Unit.hpp"

namespace Phyth {
    template <typename UnitT>
    struct Interval {
        static_assert(is_unit_v<UnitT>, "Unit type cannot be a non-unit");

        explicit Interval(Quantity<UnitT> value) noexcept
            : min_(value), max_(value) {}

        Interval(Quantity<UnitT> min, Quantity<UnitT> max)
            : min_(min), max_(max) {
            if (min_ > max_) {
                throw std::invalid_argument("Interval: min must be less than max");
            }
        }

        bool Contains(Quantity<UnitT> x)  const noexcept {
            return min_ <= x && x <= max_;
        }

        bool Contains(const Interval &interval) const noexcept {
            return interval.max_ <= max_ && interval.min_ >= min_;
        }

        Quantity<UnitT> GetMinimum() const noexcept { return min_; }
        Quantity<UnitT> GetMaximum() const noexcept { return max_; }

    private:
        Quantity<UnitT> min_;
        Quantity<UnitT> max_;
    };

}

#endif //PHYTH_INTERVAL_HPP
