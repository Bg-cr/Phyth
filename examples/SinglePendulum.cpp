#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanics;

constexpr Quantity<Second> DELTA_TIME = 1_s / 300;

int main() {
    const auto anchor_point = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor_point->SetFixed(true);
    const auto particle = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
    const DistanceConstrainer c1 {anchor_point, particle};

    constexpr auto total_time = 2.73_s;
    constexpr int steps = static_cast<int>((total_time / DELTA_TIME / 20).value);

    for (int i = 0; i < total_time / DELTA_TIME; ++i) {
        particle->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        c1.Correct(DELTA_TIME);
        particle->Integrate(DELTA_TIME);

        if (i % steps == 0) {
            // If we record the simulated data every moment, the final analysis will result in image 'SinglePendulumAnalysis.png'

            // If we record 5000 simulated data on average during a 200'000 second simulation,
            // the final analysis will generate the image 'SinglePendulumAnalysis2.png'

            std::cout << "Time: " << i * DELTA_TIME << "\n"
                      << "  Position:" << particle->GetPosition() << "\n"
                      << "  Velocity: " << particle->GetVelocity() << "\n"
                      << "  KE: " << particle->GetKineticEnergy() << "\n"
                      << "  PE: " << particle->GetPotentialEnergy(-1_m) << "\n"
                      << "  TE: " << particle->GetTotalEnergy(-1_m) << "\n"
                      << std::endl;
        }
    }
}
