#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanical;

constexpr Quantity<Second> DELTA_TIME = 1_s / 300;

int main() {
    const auto anchor_point = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor_point->SetFixed(true);
    const auto particle = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
    const DistanceConstrainer c1 {anchor_point, particle};

    for (int i = 0; i < 2.37_s / DELTA_TIME; ++i) {
        particle->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
        c1.Correct(DELTA_TIME);
        particle->Integrate(DELTA_TIME);

        if (i % static_cast<int>((0.2_s / DELTA_TIME).value) == 0) {
            // If we record the simulated data every moment, the final analysis will result in image 'SinglePendulumAnalysis.png'
            const auto t = i * DELTA_TIME;
            const auto pos = particle->GetPosition();
            const auto vel = particle->GetVelocity();

            const auto ke = 0.5_ * 1_kg * Utils::pow<2>(vel.Length());
            const auto pe = 1_kg * Consts::g * (pos.y + 1_m);
            const auto te = ke + pe;

            std::cout << "Time: " << t << "\n"
                      << "  Position:" << pos << "\n"
                      << "  Velocity: " << vel << "\n"
                      << "  KE: " << ke << "\n"
                      << "  PE: " << pe << "\n"
                      << "  TE: " << te << "\n"
                      << std::endl;
        }
    }
}
