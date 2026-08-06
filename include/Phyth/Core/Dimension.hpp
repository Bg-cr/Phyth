#ifndef PHYTH_DIMENSION_HPP
#define PHYTH_DIMENSION_HPP

#include <ratio>
#include <type_traits>
#include <string>
#include <tuple>
#include <cmath>

namespace Phyth {
    /**
     * @brief Compile-time physical dimension with 7 SI base exponents
     * @tparam Ratios std::ratio sequence for [m, kg, s, A, K, mol, cd]
     */
    template<typename... Ratios>
    struct Dimension {
        static constexpr int size = sizeof...(Ratios);
        using ratio_tuple = std::tuple<Ratios...>;
    };

    /**
     * @brief Extracts the N-th ratio type from a Dimension
     * @tparam N Index of the ratio to extract (0-based)
     * @tparam D The Dimension type
     */
    template<int N, typename D>
    struct DimComponent;

    template<int N, typename... Ratios>
    struct DimComponent<N, Dimension<Ratios...>> {
        static_assert(N < sizeof...(Ratios), "Index out of range");
        using type = std::tuple_element_t<N, std::tuple<Ratios...>>;
    };

    template<int N, typename D>
    using DimComponentT = typename DimComponent<N, D>::type;

    template<typename D1, typename D2>
    struct DimMul;

    template<typename... R1, typename... R2>
    struct DimMul<Dimension<R1...>, Dimension<R2...> > {
        static_assert(sizeof...(R1) == sizeof...(R2), "Dimension count mismatch");
        using type = Dimension<std::ratio_add<R1, R2>...>;
    };

    template<typename D1, typename D2>
    using DimMulT = typename DimMul<D1, D2>::type;

    template<typename D1, typename D2>
    struct DimDiv;

    template<typename... R1, typename... R2>
    struct DimDiv<Dimension<R1...>, Dimension<R2...> > {
        static_assert(sizeof...(R1) == sizeof...(R2), "Dimension count mismatch");
        using type = Dimension<std::ratio_subtract<R1, R2>...>;
    };

    template<typename D1, typename D2>
    using DimDivT = typename DimDiv<D1, D2>::type;

    template<int N, typename D>
    struct DimRoot;

    template<int N, typename... Ratios>
    struct DimRoot<N, Dimension<Ratios...> > {
        static_assert(N > 0, "Root index must be positive");
        using type = Dimension<std::ratio_divide<Ratios, std::ratio<N> >...>;
    };

    template<int N, typename D>
    using DimRootT = typename DimRoot<N, D>::type;

    template<typename D>
    struct DimRcp;

    template<typename... Ratios>
    struct DimRcp<Dimension<Ratios...> > {
        using type = Dimension<std::ratio_multiply<Ratios, std::ratio<-1> >...>;
    };

    template<typename D>
    using DimRcpT = typename DimRcp<D>::type;

    template<typename D1, typename D2>
    struct DimEqual;

    template<typename... R1, typename... R2>
    struct DimEqual<Dimension<R1...>, Dimension<R2...> > {
        static constexpr bool value = (std::is_same_v<R1, R2> && ...);
    };

    template<typename D1, typename D2>
    inline constexpr bool dim_equal_v = DimEqual<D1, D2>::value;

    template<typename>
    struct is_dim : std::false_type {
    };

    template<typename... Ratios>
    struct is_dim<Dimension<Ratios...> > : std::true_type {
    };

    template<typename T>
    inline constexpr bool is_dim_v = is_dim<T>::value;

    inline std::string formatPower(const double power) {
        if (std::floor(power + 0.5) == power) {
            // floor(power + 0.5) avoids floating point artifacts like 1.000000 -> 1
            return std::to_string(static_cast<int>(power));
        }
        return std::to_string(power);
    }

    template<typename D>
    struct DimToString;

    template<typename... Ratios>
    struct DimToString<Dimension<Ratios...> > {
        static std::string value() {
            const std::string names[] = {"m", "kg", "s", "A", "K", "mol", "cd"};

            std::string numerator, denominator;
            bool has_num = false, has_den = false;
            int idx = 0;

            ((void) [&] {
                const double power = static_cast<double>(Ratios::num) / Ratios::den;
                if (power == 0) {
                    ++idx;
                    return;
                }

                const std::string unit = names[idx];
                std::string power_str;
                if (power != 1 && power != -1) {
                    power_str = formatPower(power < 0 ? -power : power);
                }

                if (power > 0) {
                    if (has_num) numerator += '*';
                    numerator += unit;
                    if (power != 1) {
                        numerator += "^" + power_str;
                    }
                    has_num = true;
                } else {
                    if (has_den) denominator += '*';
                    denominator += unit;
                    if (power != -1) {
                        denominator += "^" + power_str;
                    }
                    has_den = true;
                }
                ++idx;
            }(), ...);

            std::string result;
            if (has_num) {
                result += numerator;
            } else if (has_den) {
                result += '1';
            }

            if (has_den) {
                result += '/';
                if (denominator.find('*') != std::string::npos) {
                    result += "(" + denominator + ")";
                } else {
                    result += denominator;
                }
            }

            return result;
        }
    };

    using Dimensionless = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<0>,
        std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<0>
    >;

    using Length = Dimension<
        std::ratio<1>, std::ratio<0>, std::ratio<0>,
        std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<0>
    >;

    using Mass = Dimension<
        std::ratio<0>, std::ratio<1>, std::ratio<0>,
        std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<0>
    >;

    using Time = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<1>,
        std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<0>
    >;

    using ElectricCurrent = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<0>,
        std::ratio<1>, std::ratio<0>, std::ratio<0>, std::ratio<0>
    >;

    using ThermodynamicTemperature = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<0>,
        std::ratio<0>, std::ratio<1>, std::ratio<0>, std::ratio<0>
    >;

    using AmountOfSubstance = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<0>,
        std::ratio<0>, std::ratio<0>, std::ratio<1>, std::ratio<0>
    >;

    using LuminousIntensity = Dimension<
        std::ratio<0>, std::ratio<0>, std::ratio<0>,
        std::ratio<0>, std::ratio<0>, std::ratio<0>, std::ratio<1>
    >;

    using Velocity = DimDivT<Length, Time>;
    using Acceleration = DimDivT<Velocity, Time>;
    using Force = DimMulT<Mass, Acceleration>;
    using Energy = DimMulT<Force, Length>;
    using Power = DimDivT<Energy, Time>;
    using Pressure = DimDivT<Force, DimMulT<Length, Length> >;
    using Frequency = DimDivT<Dimensionless, Time>;
    using Charge = DimMulT<Time, ElectricCurrent>;
    using ElectricFieldIntensity = DimDivT<Force, Charge>;
    using Voltage = DimDivT<Power, ElectricCurrent>;
    using Resistance = DimDivT<Voltage, ElectricCurrent>;
    using Capacitance = DimDivT<Charge, Voltage>;
    using Inductance = DimDivT<Voltage, DimDivT<ElectricCurrent, Time> >;
    using MagneticFlux = DimMulT<Voltage, Time>;
    using MagneticFluxDensity = DimDivT<MagneticFlux, DimMulT<Length, Length> >;
    using Area = DimMulT<Length, Length>;
    using Volume = DimMulT<Length, Area>;
    using TranslationalStiffness = DimDivT<Force, Length>;
    using DampingCoefficient = DimDivT<Mass, Time>;
    using MassDensity = DimDivT<Mass, Volume>;
    using ElectricPotential = DimDivT<Capacitance, Length>;
    using EnergyDensity = DimDivT<Energy, Volume>;
    using LinearChargeDensity = DimDivT<Charge, Length>;
    using SurfaceChargeDensity = DimDivT<Charge, Area>;
    using BulkChargeDensity = DimDivT<Charge, Volume>;
    using ElectricDipoleMoment = DimMulT<Charge, Length>;
    using MagneticDipoleMoment = DimMulT<ElectricCurrent, Area>;
    using MomentOfInertia = DimMulT<Mass, Area>;
    using AngularVelocity = DimDivT<Dimensionless, Time>;

    template<typename UnitT>
    struct is_dimensionless {
        static constexpr bool value = std::is_same_v<typename UnitT::DimensionT, Dimensionless>;
    };

    template<typename UnitT>
    inline constexpr bool is_dimensionless_v = is_dimensionless<UnitT>::value;

    template<int N, typename D>
    struct DimPower;

    template<int N, typename... Ratios>
    struct DimPower<N, Dimension<Ratios...> > {
        static constexpr double value = static_cast<double>(
                                            std::tuple_element_t<N, std::tuple<Ratios...> >::num
                                        ) / std::tuple_element_t<N, std::tuple<Ratios...> >::den;
    };

    template<int N, typename D>
    inline constexpr double dim_power_v = DimPower<N, D>::value;

    template<int P, typename D>
    struct DimPowerHelper {
        using type = DimMulT<D, typename DimPowerHelper<P - 1, D>::type>;
    };

    template<typename D>
    struct DimPowerHelper<0, D> {
        using type = Dimensionless;
    };

    template<int Power, typename D>
    struct DimPowerType {
        using type = typename DimPowerHelper<Power, D>::type;
    };

    template<int Power, typename D>
    using DimPowerTypeT = typename DimPowerType<Power, D>::type;
}

#endif //PHYTH_DIMENSION_HPP
