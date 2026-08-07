#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanics;


int main() {
    Config::dt = 1_s / 1000;

    const auto anchor_point = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor_point->SetFixed(true);
    const auto p1 = std::make_shared<Particle>(1_kg, Vector3{1_m, 1_m, 0_m});
    const DistanceConstrainer c1 {anchor_point, p1};
    const auto p2 = std::make_shared<Particle>(1_kg, Vector3{2_m, 0_m, 0_m});
    const DistanceConstrainer c2 {p1, p2};
    const DistanceConstrainer c3 {p2, p1};

    p1->SetComputeForcesFunction(
        [](Particle * p) {
            p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        }
    );

    p2->SetComputeForcesFunction(
        [](Particle *p) {
            p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        }
    );

    constexpr auto total_time = 10000_s;
    const int steps = (total_time / Config::dt / 50).to<int>();

    for (int i = 0; i < total_time / Config::dt; ++i) {
        c1.Correct();
        c2.Correct();
        c3.Correct();

        p1->Integrate();
        p2->Integrate();

        if (i % steps == 0) {
            // If we record the simulated data every moment, the final analysis will result in image 'DoublePendulumAnalysis.png'
            std::cout << "Time: " << i * Config::dt << "\n"
                      << "  P1 Position:" << p1->GetPosition() << "\n"
                      << "  P1 Velocity: " << p1->GetVelocity() << "\n"
                      << "  P2 Position:" << p2->GetPosition() << "\n"
                      << "  P2 Velocity: " << p2->GetVelocity() << "\n"
                      << "  KE: " <<  p1->GetKineticEnergy() + p2->GetKineticEnergy() << "\n"
                      << "  PE: " << p1->GetPotentialEnergy(-2.828_m) + p2->GetPotentialEnergy(-2.828_m) << "\n"
                      << "  TE: " << p1->GetTotalEnergy(-2.828_m) + p2->GetTotalEnergy(-2.828_m) << "\n"
                      << std::endl;
        }
    }
}
