# Phyth
> A dimensionally safe C++ header-only mechanics and electromagnetics simulation library.

[![CI](https://github.com/Bg-cr/Phyth/actions/workflows/ci.yml/badge.svg)](https://github.com/Bg-cr/Phyth/actions)
![License](https://img.shields.io/github/license/Bg-cr/Phyth)
![C++](https://img.shields.io/badge/C++-17-blue.svg)
![Header-only](https://img.shields.io/badge/header--only-yes-2ea44f.svg)
[![Code Size](https://img.shields.io/github/languages/code-size/Bg-cr/Phyth)](https://github.com/Bg-cr/Phyth)

## Table of Contents

- [Description](#description)
- [Features](#features)
- [Positioning](#positioning)
- [Installation](#installation)
- [Quick Start](#quick-start)

## Description
Phyth is a C++ simulation library designed for numerical method researchers. 
While maintaining lightweight and safety, it breaks down the simulation process into replaceable independent components,
allowing you to change integrators or constraint solvers like building blocks, 
without being troubled by the complex architecture of large engines.

## Features
- **Type-safe physical quantities**: Literals like `1_kg`, `1_m` make code concise and readable, 
while the type system catches dimension mismatches at compile time—assigning meters to seconds, or using force as mass, 
fails before runtime, eliminating unit-confusion bugs.

- **Semantic unit conversion**: Units of the same dimension (e.g., joules and newton-meters) retain distinct semantic labels. 
Conversions preserve meaning, and output formatting automatically distinguishes `J` from `N*m`—no redundant re-declaration, no ambiguity.

- **Straightforward simulation workflow**: Create particles → apply constraints → define forces → loop over integration and correction. 
The simulation loop follows the same sequence as a physics textbook, making the learning curve remarkably shallow.

- **Interchangeable components**: Constraint solvers (e.g., `DistanceConstraint`) and force models (e.g., `GravityFormula`) are designed as pluggable, 
independent components—swap integrators or implement custom algorithms without restructuring the entire framework.

- **Zero dependencies, lightweight core**: The core library has no third-party dependencies and compiles to a small binary, 
suitable for embedded environments, teaching demonstrations, and research scenarios that demand deep customization.

- **Built-in validation hooks**: Energy monitoring interfaces (e.g., `GetTotalEnergy()`) are available out of the box, 
making it easy to verify physical conservation laws during simulation and simplifying debugging and experiment validation.

## Positioning

Existing physics engines sit at two extremes:

|                   | Industrial Engines            | Scientific Engines                |
|-------------------|-------------------------------|-----------------------------------|
| **Speed**         | ✅ Fast                       | ⚠️ Slow                           |
| **Robustness**    | ✅ Battle-tested              | ⚠️ Academic                       |
| **Transparency**  | ❌ Black box (stages fused)   | ✅ Exposed internals              |
| **Modifiability** | ❌ Fixed solver, hard to swap | ⚠️ Heavy codebase, hard to modify |
| **Dependencies**  | ✅ Minimal                    | ❌ Heavy dependency graph         |

**Phyth sits right in the middle:**

> **Transparent like a scientific engine, lightweight like an industrial one.**

- **More transparent** than industrial engines — simulation stages are exposed as replaceable components, 
so you can swap integrators or constraint solvers without fighting the framework.
- **More lightweight** than scientific engines — zero dependencies, minimal codebase, 
easy to build, modify, and experiment with.

**In short**, Phyth gives researchers the freedom to develop new numerical methods, 
without the friction of production-level complexity or academic-level baggage.

## Installation

### Prerequisites

- C++17 or higher

#### Clone the Repository

```bash
git clone https://github.com/Bg-cr/Phyth.git
```

#### Include in Your Project

Add the `include/` directory to your compiler's header search path:

```cpp
#include <Phyth/Phyth.hpp>
```

Done. No compilation required.

## Quick Start

The following code demonstrates **single pendulum motion under gravity**.
The output results consist of the particle and total mechanical energy for each time step,
which are used to verify numerical stability.

```cpp
#include "Phyth/Phyth.hpp"
#include <iostream>

using namespace Phyth;
using namespace Phyth::Mechanics;

constexpr Quantity<Second> DELTA_TIME = 1_s / 1000;
constexpr Quantity<Second> SIMULATION_TIME = 2.73_s;

int main() {
    const auto anchor = std::make_shared<Particle>(1_kg, Vector3{0_m, 0_m, 0_m});
    anchor->SetFixed(true);
    
    const auto pendulum = std::make_shared<Particle>(1_kg, Vector3{1_m, 0_m, 0_m});
    
    const DistanceConstrainer constraint{anchor, pendulum};

    pendulum->SetComputeForcesFunction([](Particle *p) {
        p->ApplyForce({0_N, -1_kg * Consts::g, 0_N});
    });

    for (Quantity<Second> t = 0_s; t < SIMULATION_TIME; t += DELTA_TIME) {
        pendulum->Integrate(DELTA_TIME);
        constraint.Correct(DELTA_TIME);

        std::cout << t << " " << *pendulum << " " << pendulum->GetTotalEnergy(-1_m) << std::endl;
    }

    return 0;
}
```

**Expected output example**
```
0 s Particle{pos=(1 m, -4.90332e-06 m, 0 m), mass=1 kg} 9.80665 J
0.001 s Particle{pos=(1 m, -1.96133e-05 m, 0 m), mass=1 kg} 9.80665 J
0.002 s Particle{pos=(1 m, -4.41299e-05 m, 0 m), mass=1 kg} 9.80665 J
0.003 s Particle{pos=(1 m, -7.84532e-05 m, 0 m), mass=1 kg} 9.80665 J
0.004 s Particle{pos=(1 m, -0.000122583 m, 0 m), mass=1 kg} 9.80665 J
0.005 s Particle{pos=(1 m, -0.00017652 m, 0 m), mass=1 kg} 9.80665 J
0.006 s Particle{pos=(1 m, -0.000240263 m, 0 m), mass=1 kg} 9.80665 J
0.007 s Particle{pos=(1 m, -0.000313813 m, 0 m), mass=1 kg} 9.80665 J
0.008 s Particle{pos=(1 m, -0.000397169 m, 0 m), mass=1 kg} 9.80665 J
0.009 s Particle{pos=(1 m, -0.000490332 m, 0 m), mass=1 kg} 9.80665 J
0.01 s Particle{pos=(1 m, -0.000593302 m, 0 m), mass=1 kg} 9.80665 J
0.011 s Particle{pos=(1 m, -0.000706079 m, 0 m), mass=1 kg} 9.80665 J
0.012 s Particle{pos=(1 m, -0.000828662 m, 0 m), mass=1 kg} 9.80665 J
0.013 s Particle{pos=(1 m, -0.000961052 m, 0 m), mass=1 kg} 9.80665 J
0.014 s Particle{pos=(0.999999 m, -0.00110325 m, 0 m), mass=1 kg} 9.80665 J
0.015 s Particle{pos=(0.999999 m, -0.00125525 m, 0 m), mass=1 kg} 9.80665 J
...
```

**Verification method**
If the above code is run, the total energy at each moment should be approximately `9.80665 J`.

In fact, after running and analysis the data, it can be concluded that the relative fluctuation within `2.73s` is approximately `1.32e-5`, 
which confirms that the integrator and constraint solver are working properly.


