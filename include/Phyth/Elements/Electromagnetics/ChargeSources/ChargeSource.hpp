#ifndef PHYTH_CHARGE_SOURCE_HPP
#define PHYTH_CHARGE_SOURCE_HPP

#include "Phyth/Tools/Vector3.hpp"
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Core/Units.hpp"
#include "../ChargeParticle.hpp"

namespace Phyth::Electromagnetics {
    class ChargeSource {
    public:
        virtual ~ChargeSource() = default;

        [[nodiscard]] virtual Vector3<Quantity<NewtonPerCoulomb>>
        GetElectricFieldAt(const Vector3<Quantity<Meter>>& point) const = 0;

        [[nodiscard]] virtual Quantity<Volt>
        GetElectricPotentialAt(const Vector3<Quantity<Meter>>& point) const = 0;

        [[nodiscard]] virtual Quantity<Coulomb> GetTotalChargeValue() const = 0;

        [[nodiscard]] virtual Quantity<CoulombPerCubicMeter>
        GetChargeDensityAt(const Vector3<Quantity<Meter>>& point) const {
            return 0_Cpm3;
        }

        [[nodiscard]] virtual Quantity<Meter>
        GetMinimumDistanceTo(const Vector3<Quantity<Meter>>& point) const = 0;
    };


    template <typename T>
    inline constexpr bool is_charge_source_v = std::is_base_of_v<ChargeSource, T>;
}

#endif //PHYTH_CHARGE_SOURCE_HPP
