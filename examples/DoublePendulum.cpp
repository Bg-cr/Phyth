#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanical;

constexpr Quantity<Second> DELTA_TIME = 1_s / 300;

int main() {
    const auto anchor_point = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor_point->SetFixed(true);
    const auto p1 = std::make_shared<Particle>(1_kg, Vector3{1_m, 1_m, 0_m});
    const DistanceConstrainer c1 {anchor_point, p1};
    const auto p2 = std::make_shared<Particle>(1_kg, Vector3{2_m, 0_m, 0_m});
    const DistanceConstrainer c2 {p1, p2};
    const DistanceConstrainer c3 {p2, p1};

    constexpr auto total_time = 10000_s;
    constexpr int steps = static_cast<int>((total_time / DELTA_TIME / 50).value);

    for (int i = 0; i < total_time / DELTA_TIME; ++i) {
        p1->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        p2->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        c1.Correct(DELTA_TIME);
        c2.Correct(DELTA_TIME);
        c3.Correct(DELTA_TIME);
        p1->Integrate(DELTA_TIME);
        p2->Integrate(DELTA_TIME);

        if (i % steps == 0) {
            // If we record the simulated data every moment, the final analysis will result in image 'DoublePendulumAnalysis.png'
            const auto t = i * DELTA_TIME;

            const auto pos1 = p1->GetPosition();
            const auto vel1 = p1->GetVelocity();
            const auto ke1 = 0.5_ * 1_kg * Utils::pow<2>(vel1.Length());
            const auto pe1 = 1_kg * Consts::g * (pos1.y + 2_m);

            const auto pos2 = p2->GetPosition();
            const auto vel2 = p2->GetVelocity();
            const auto ke2 = 0.5_ * 1_kg * Utils::pow<2>(vel2.Length());
            const auto pe2 = 1_kg * Consts::g * (pos2.y + 2_m);

            const auto ke_total = ke1 + ke2;
            const auto pe_total = pe1 + pe2;
            const auto te_total = ke_total + pe_total;

            std::cout << "Time: " << t << "\n"
                      << "  P1 Position:" << pos1 << "\n"
                      << "  P1 Velocity: " << vel1 << "\n"
                      << "  P2 Position:" << pos2 << "\n"
                      << "  P2 Velocity: " << vel2 << "\n"
                      << "  KE: " << ke_total << "\n"
                      << "  PE: " << pe_total << "\n"
                      << "  TE: " << te_total << "\n"
                      << std::endl;
        }
    }
}
