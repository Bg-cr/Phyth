#ifndef PHYTH_PHYSICAL_CONSTS_HPP
#define PHYTH_PHYSICAL_CONSTS_HPP

#include "../Core/Units.hpp"
#include "Constants.hpp"

namespace Phyth::Consts {
    constexpr double ALPHA = 0.00729735256;

    constexpr auto g = 9.80665_mps2;
    constexpr auto G = 6.67430e-11_N * 1_m2 / 1_kg / 1_kg;
    constexpr auto c = 299'792'458_mps;
    constexpr auto h = 6.62607015e-34_J * 1_s;
    constexpr auto hbar = h / TAU;
    constexpr auto k_B = 1.380649e-23_J / 1_K;
    constexpr auto MU_0 = VARPI * 1e-7_N / 1_A / 1_A;
    constexpr auto EPSILON_0 = 8.8541878e-12_Fpm;
    constexpr auto k_E = 1 / VARPI / EPSILON_0;
}

#endif //PHYTH_PHYSICAL_CONSTS_HPP
