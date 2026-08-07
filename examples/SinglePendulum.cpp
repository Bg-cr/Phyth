#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanics;


int main() {
    Config::dt = 1_s / 1000;
    
    const auto anchor_point = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor_point->SetFixed(true);
    const auto particle = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
    const DistanceConstrainer c1 {anchor_point, particle};

    particle->SetComputeForcesFunction([](Particle *p){p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});});

    constexpr auto total_time = 2.73_s;
    const int steps = (total_time / Config::dt / 20).to<int>();

    for (int i = 0; i < total_time / Config::dt; ++i) {
        particle->Integrate();
        c1.Correct();

        if (i % steps == 0) {
            // If we record the simulated data every moment, the final analysis will result in image 'SinglePendulumAnalysis.png'

            // If we record 5000 simulated data on average during a 200'000 second simulation,
            // the final analysis will generate the image 'SinglePendulumAnalysis2.png'

            std::cout << "Time: " << i * Config::dt << "\n"
                      << "  Position:" << particle->GetPosition() << "\n"
                      << "  Velocity: " << particle->GetVelocity() << "\n"
                      << "  KE: " << particle->GetKineticEnergy() << "\n"
                      << "  PE: " << particle->GetPotentialEnergy(-1_m) << "\n"
                      << "  TE: " << particle->GetTotalEnergy(-1_m) << "\n"
                      << std::endl;
        }
    }
}
