#ifndef PHYTH_UNITS_H
#define PHYTH_UNITS_H

#include "Phyth/Core/Dimension.hpp"
#include "Phyth/System/UnitsRegistry.hpp"
#include "Unit.hpp"
#include "Phyth/Physical/Constants.hpp"

namespace Phyth {
    #define PHYTH_CONCAT_IMPL(a, b) a##b
    #define PHYTH_CONCAT(a, b) PHYTH_CONCAT_IMPL(a, b)

    #define PHYTH_DEFINE_UNIT(dim_type, unit_type, suffix, scale_num, scale_den) \
        using unit_type = Unit<dim_type, std::ratio<scale_num, scale_den>>; \
        \
        constexpr auto operator""_##suffix(long double x) { \
            return Quantity<unit_type>(static_cast<double>(x) * unit_type::scale); \
        } \
        \
        constexpr auto operator""_##suffix(unsigned long long x) { \
            return Quantity<unit_type>(static_cast<double>(x) * unit_type::scale); \
        } \
        \
        static const Phyth::UnitRegistrar<unit_type> PHYTH_CONCAT(registrar_, unit_type)(#unit_type, #suffix);

    #define PHYTH_DEFINE_UNIT_WITH_NAME(dim_type, unit_type, suffix, scale_num, scale_den, symbol) \
        using unit_type = Unit<dim_type, std::ratio<scale_num, scale_den>>; \
        \
        constexpr auto operator""_##suffix(long double x) { \
            return Quantity<unit_type>(static_cast<double>(x) * unit_type::scale); \
        } \
        \
        constexpr auto operator""_##suffix(unsigned long long x) { \
            return Quantity<unit_type>(static_cast<double>(x) * unit_type::scale); \
        } \
        \
        static const Phyth::UnitRegistrar<unit_type> PHYTH_CONCAT(registrar_, unit_type)(#unit_type, symbol);

    // ScalarUnit & Radian, Dimensionless * 1.0

    PHYTH_DEFINE_UNIT(Dimensionless, ScalarUnit, , 1, 1)

    struct RadianTag {};
    using Radian = Unit<Dimensionless, std::ratio<1>, RadianTag>;
    constexpr auto operator""_rad(long double x) {
        return Quantity<Radian>(static_cast<double>(x));
    }
    constexpr auto operator""_rad(unsigned long long x) {
        return Quantity<Radian>(static_cast<double>(x));
    }
    static const UnitRegistrar<Radian> registrar_Radian("radian", "rad");

    struct DegreeTag {};
    using Degree = Unit<Dimensionless, std::ratio<1>, DegreeTag>;
    constexpr auto operator""_deg(long double x) {
        return Quantity<Degree>(static_cast<double>(x) * Consts::PI / 180.0);
    }
    constexpr auto operator""_deg(unsigned long long x) {
        return Quantity<Degree>(static_cast<double>(x) * Consts::PI / 180.0);
    }
    static const UnitRegistrar<Degree> registrar_Degree("degree", "deg");

    PHYTH_DEFINE_UNIT(Length, Meter, m, 1, 1)
    PHYTH_DEFINE_UNIT(Length, Kilometer, km, 1000, 1)
    PHYTH_DEFINE_UNIT(Length, Centimeter, cm, 1, 100)
    PHYTH_DEFINE_UNIT(Length, Millimeter, mm, 1, 1000)
    PHYTH_DEFINE_UNIT(Length, Foot, ft, 3048, 10000)
    PHYTH_DEFINE_UNIT(Length, Inch, in, 254, 10000)
    PHYTH_DEFINE_UNIT(Length, Mile, mi, 1609344, 1000)

    PHYTH_DEFINE_UNIT(Mass, Kilogram, kg, 1, 1)
    PHYTH_DEFINE_UNIT(Mass, Gram, g, 1, 1000)

    PHYTH_DEFINE_UNIT(Time, Second, s, 1, 1)
    PHYTH_DEFINE_UNIT(Time, Millisecond, ms, 1, 1000)
    PHYTH_DEFINE_UNIT(Time, Minute, min, 60, 1)
    PHYTH_DEFINE_UNIT(Time, Hour, h, 3600, 1)

    PHYTH_DEFINE_UNIT(ElectricCurrent, Ampere, A, 1, 1)
    PHYTH_DEFINE_UNIT(ThermodynamicTemperature, Kelvin, K, 1, 1)
    PHYTH_DEFINE_UNIT(AmountOfSubstance, Mole, mol, 1, 1)
    PHYTH_DEFINE_UNIT(LuminousIntensity, Candela, cd, 1, 1)

    PHYTH_DEFINE_UNIT(Force, Newton, N, 1, 1)
    PHYTH_DEFINE_UNIT(Force, Kilonewton, kN, 1000, 1)

    // Joule & NewtonMeter, Energy * 1.0

    PHYTH_DEFINE_UNIT(Energy, Joule, J, 1, 1)

    struct NewtonMeterTag {};
    using NewtonMeter = Unit<Energy, std::ratio<1>, NewtonMeterTag>;
    constexpr auto operator""_Nm(long double x) {
        return Quantity<NewtonMeter>(static_cast<double>(x));
    }
    constexpr auto operator""_Nm(unsigned long long x) {
        return Quantity<NewtonMeter>(static_cast<double>(x));
    }
    static const UnitRegistrar<NewtonMeter> registrar_NewtonMeter("NewtonMeter", "N*m");

    PHYTH_DEFINE_UNIT(Energy, Kilojoule, kJ, 1000, 1)
    PHYTH_DEFINE_UNIT(Energy, ElectronVolt, eV, 160217663, 1000000000000000000LL)

    PHYTH_DEFINE_UNIT(Power, Watt, W, 1, 1)
    PHYTH_DEFINE_UNIT(Power, Kilowatt, kW, 1000, 1)

    PHYTH_DEFINE_UNIT(Pressure, Pascal, Pa, 1, 1)
    PHYTH_DEFINE_UNIT(Pressure, Kilopascal, kPa, 1000, 1)
    PHYTH_DEFINE_UNIT(Pressure, Bar, bar, 100000, 1)

    PHYTH_DEFINE_UNIT(Frequency, Hertz, Hz, 1, 1)
    PHYTH_DEFINE_UNIT(Frequency, Kilohertz, kHz, 1000, 1)

    PHYTH_DEFINE_UNIT(Charge, Coulomb, C, 1, 1)

    PHYTH_DEFINE_UNIT(Voltage, Volt, V, 1, 1)
    PHYTH_DEFINE_UNIT(Voltage, Kilovolt, kV, 1000, 1)

    PHYTH_DEFINE_UNIT(Resistance, Ohm, Ohm, 1, 1)
    PHYTH_DEFINE_UNIT(Resistance, Kiloohm, kOhm, 1000, 1)

    PHYTH_DEFINE_UNIT(Capacitance, Farad, F, 1, 1)
    PHYTH_DEFINE_UNIT(Capacitance, Microfarad, uF, 1, 1000000)

    PHYTH_DEFINE_UNIT(Inductance, Henry, H, 1, 1)
    PHYTH_DEFINE_UNIT(Inductance, Millihenry, mH, 1, 1000)

    PHYTH_DEFINE_UNIT(MagneticFlux, Weber, Wb, 1, 1)
    PHYTH_DEFINE_UNIT(MagneticFluxDensity, Tesla, T, 1, 1)

    PHYTH_DEFINE_UNIT_WITH_NAME(Velocity, MeterPerSecond, mps, 1, 1, "m/s")
    PHYTH_DEFINE_UNIT_WITH_NAME(Velocity, KilometerPerHour, kmph, 1, 36, "km/h")

    PHYTH_DEFINE_UNIT_WITH_NAME(Acceleration, MeterPerSecondSquared, mps2, 1, 1, "m/s^2")

    PHYTH_DEFINE_UNIT_WITH_NAME(Area, MeterSquared, m2, 1, 1, "m^2")
    PHYTH_DEFINE_UNIT_WITH_NAME(Volume, MeterCubed, m3, 1, 1, "m^3")
    PHYTH_DEFINE_UNIT_WITH_NAME(TranslationalStiffness, NewtonPerMeter, Npm, 1, 1, "N/m")
    PHYTH_DEFINE_UNIT_WITH_NAME(DampingCoefficient, KilogramPerSecond, kgps, 1, 1, "kg/s")
    PHYTH_DEFINE_UNIT_WITH_NAME(MassDensity, KilogramPerMeter, kgpm, 1, 1, "kg/m")

    #undef PHYTH_DEFINE_UNIT
    #undef PHYTH_DEFINE_UNIT_WITH_NAME

}

#endif  //PHYTH_UNITS_H