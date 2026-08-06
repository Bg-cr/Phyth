#ifndef PHYTH_UNITS_HPP
#define PHYTH_UNITS_HPP

#include "Phyth/Core/Dimension.hpp"
#include "Phyth/System/UnitsRegistry.hpp"
#include "Unit.hpp"
#include "Phyth/Core/Quantity.hpp"

namespace Phyth {
#define PHYTH_CONCAT_IMPL(a, b) a##b
#define PHYTH_CONCAT(a, b) PHYTH_CONCAT_IMPL(a, b)

    /**
     * @def PHYTH_DEFINE_UNIT
     * @brief Define a unit, literal operator, and registry entry
     * @param dim_type Dimension type
     * @param unit_type Unit type name
     * @param suffix Literal suffix
     * @param scale_num Scale numerator
     * @param scale_den Scale denominator
     */
#define PHYTH_DEFINE_UNIT(dim_type, unit_type, suffix, scale_num, scale_den) \
    using unit_type = Unit<dim_type, std::ratio<scale_num, scale_den>>; \
    \
    constexpr auto operator""_##suffix(long double x) { \
        return Quantity<unit_type>(x); \
    } \
    \
    constexpr auto operator""_##suffix(unsigned long long x) { \
        return Quantity<unit_type>(x); \
    } \
    \
    static const Phyth::UnitRegistrar<unit_type> PHYTH_CONCAT(registrar_, unit_type)(#unit_type, #suffix);

    /**
     * @def PHYTH_DEFINE_UNIT_WITH_NAME
     * @brief Same as DEFINE_UNIT with custom display symbol
     * @param dim_type Dimension type
     * @param unit_type Unit type name
     * @param suffix Literal suffix
     * @param scale_num Scale numerator
     * @param scale_den Scale denominator
     * @param symbol Display symbol
     */
#define PHYTH_DEFINE_UNIT_WITH_NAME(dim_type, unit_type, suffix, scale_num, scale_den, symbol) \
    using unit_type = Unit<dim_type, std::ratio<scale_num, scale_den>>; \
    \
    constexpr auto operator""_##suffix(long double x) { \
        return Quantity<unit_type>(x); \
    } \
    \
    constexpr auto operator""_##suffix(unsigned long long x) { \
        return Quantity<unit_type>(x); \
    } \
    \
    static const Phyth::UnitRegistrar<unit_type> PHYTH_CONCAT(registrar_, unit_type)(#unit_type, symbol);

    /**
     * @def PHYTH_DEFINE_UNIT_WITH_TAG
     * @brief Define a unit with a tag type for type-safe differentiation
     * @param dim_type Dimension type
     * @param unit_type Unit type name
     * @param suffix Literal suffix
     * @param scale_num Scale numerator
     * @param scale_den Scale denominator
     * @param symbol Display symbol
     * @param tag Tag type
     */
#define PHYTH_DEFINE_UNIT_WITH_TAG(dim_type, unit_type, suffix, scale_num, scale_den, symbol, tag) \
        using unit_type = Unit<dim_type, std::ratio<scale_num, scale_den>, tag>; \
        \
        constexpr auto operator""_##suffix(long double x) { \
            return Quantity<unit_type>(x); \
        } \
        \
        constexpr auto operator""_##suffix(unsigned long long x) { \
            return Quantity<unit_type>(x); \
        } \
        \
        static const Phyth::UnitRegistrar<unit_type> PHYTH_CONCAT(registrar_, unit_type)(#unit_type, symbol);

    PHYTH_DEFINE_UNIT(Dimensionless, ScalarUnit, , 1, 1)

    struct RadianTag {
    };

    PHYTH_DEFINE_UNIT_WITH_TAG(Dimensionless, Radian, rad, 1, 1, "rad", RadianTag)

    struct DegreeTag {
    };

    PHYTH_DEFINE_UNIT_WITH_TAG(Dimensionless, Degree, deg, 3141592653589793, 18000000000000000, "deg", DegreeTag)

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

    PHYTH_DEFINE_UNIT(Energy, Joule, J, 1, 1)

    struct NewtonMeterTag {
    };

    PHYTH_DEFINE_UNIT_WITH_TAG(Energy, NewtonMeter, Nm, 1, 1, "N*m", NewtonMeterTag)

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
    PHYTH_DEFINE_UNIT_WITH_NAME(MassDensity, KilogramPerMeter, kgpm, 1, 1, "kg/m^3")

    PHYTH_DEFINE_UNIT_WITH_NAME(ElectricFieldIntensity, NewtonPerCoulomb, NpC, 1, 1, "N/C")

    struct VoltagePerMeterTag {
    };

    PHYTH_DEFINE_UNIT_WITH_TAG(ElectricFieldIntensity, VoltagePerMeter, Vpm, 1, 1, "V/m", VoltagePerMeterTag)

    PHYTH_DEFINE_UNIT_WITH_NAME(LinearChargeDensity, CoulombPerMeter, Cpm, 1, 1, "C/m")
    PHYTH_DEFINE_UNIT_WITH_NAME(SurfaceChargeDensity, CoulombPerSquareMeter, Cpm2, 1, 1, "C/m^2")
    PHYTH_DEFINE_UNIT_WITH_NAME(BulkChargeDensity, CoulombPerCubicMeter, Cpm3, 1, 1, "C/m^3")

    PHYTH_DEFINE_UNIT_WITH_NAME(ElectricPotential, FaradPerMeter, Fpm, 1, 1, "F/m")
    PHYTH_DEFINE_UNIT_WITH_NAME(ElectricDipoleMoment, CoulombMeter, Cm, 1, 1, "C*m")
    PHYTH_DEFINE_UNIT_WITH_NAME(MagneticDipoleMoment, AmpereMeterSquared, Am2, 1, 1, "A*m^2")
    PHYTH_DEFINE_UNIT_WITH_NAME(MomentOfInertia, KilogramMeterSquared, kgm2, 1, 1, "kg*m^2")

    PHYTH_DEFINE_UNIT_WITH_NAME(AngularVelocity, RadianPerSecond, radps, 1, 1, "rad/s");

    struct DegreePerSecondTag {
    };
    PHYTH_DEFINE_UNIT_WITH_TAG(AngularVelocity, DegreePerSecond, degps, 3141592653589793, 18000000000000000, "deg/m", DegreePerSecondTag)

#undef PHYTH_DEFINE_UNIT
#undef PHYTH_DEFINE_UNIT_WITH_NAME
}

#endif  //PHYTH_UNITS_HPP
