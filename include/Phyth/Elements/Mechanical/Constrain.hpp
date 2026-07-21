#ifndef PHYTH_CONSTRAINER_HPP
#define PHYTH_CONSTRAINER_HPP

#include <memory>
#include <utility>

#include "Particle.hpp"
#include "Phyth/Core/Interval.hpp"
#include "Phyth/Core/Unit.hpp"

namespace Phyth::Mechanical {

    class DistanceConstrainer {
    public:
        DistanceConstrainer(const Vector3<Quantity<Meter>> &ref_point, std::shared_ptr<Particle> particle)
            : ref_(ref_point), particle_(std::move(particle)), interval_((ref_ - particle_->GetPosition()).Length()) {
            IntervalCheck();
        }

        DistanceConstrainer(const Vector3<Quantity<Meter>> &ref_point, std::shared_ptr<Particle> particle, const Interval<Meter> &interval)
            : ref_(ref_point), particle_(std::move(particle)), interval_(interval) {
            IntervalCheck();
        }

        void Correct(const Quantity<Second> dt) const {
            const Vector3<Quantity<Meter>> difference = particle_->GetPosition() - ref_;
            const Quantity<Meter> length = difference.Length();
            if (interval_.Contains(length))
                return;

            const Vector3<Quantity<Meter>> target_offset = difference.Normalized() *
                (length > interval_.GetMaximum() ? interval_.GetMaximum() : interval_.GetMinimum());

            particle_->SetVelocity(particle_->GetVelocity() - (difference - target_offset) / dt);
        }

    private:
        Vector3<Quantity<Meter>> ref_;
        std::shared_ptr<Particle> particle_;

        Interval<Meter> interval_;

        void IntervalCheck() const {
            if (interval_.GetMinimum() < Quantity<Meter>(0)) {
                throw std::invalid_argument("The minimum of interval must be non-negative");
            }
        }
    };
}

#endif //PHYTH_CONSTRAINER_HPP
