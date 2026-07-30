#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanics;

constexpr Quantity<Second> DELTA_TIME = 1_s / 1000;

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

            std::cout << "Time: " << i * DELTA_TIME << "\n"
                      << "  P1 Position:" << p1->GetPosition() << "\n"
                      << "  P1 Velocity: " << p1->GetVelocity() << "\n"
                      << "  P2 Position:" << p2->GetPosition() << "\n"
                      << "  P2 Velocity: " << p2->GetVelocity() << "\n"
                      << "  KE: " <<  p1->GetKineticEnergy() + p2->GetKineticEnergy() << "\n"
                      << "  PE: " << p1->GetPotentialEnergy(-2_m) + p2->GetPotentialEnergy(-2_m) << "\n"
                      << "  TE: " << p1->GetTotalEnergy() + p2->GetTotalEnergy() << "\n"
                      << std::endl;
        }
    }
}
