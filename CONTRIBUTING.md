# Contributing to Phyth

Thanks for your interest in contributing to Phyth!
Whether you're a researcher, student, or developer, all contributions are welcome.

---

## Table of Contents

- [Ways to Contribute](#ways-to-contribute)
- [Code Style](#code-style)
- [How to Add a New Feature](#how-to-add-a-new-feature)
- [Validation Requirements](#validation-requirements)
- [Commit Messages](#commit-messages)
- [Pull Request Process](#pull-request-process)
- [Getting Help](#getting-help)

---

## Ways to Contribute

- **Report a bug** → Use the [Bug Report template](https://github.com/Bg-cr/Phyth/issues/new/choose)
- **Propose a feature** → Use the [Feature Request template](https://github.com/Bg-cr/Phyth/issues/new/choose)
- **Add a new physics module** → See [How to Add a New Feature](#how-to-add-a-new-feature)
- **Improve documentation** → Fix typos, add examples, clarify comments
- **Share your use case** → Start a [Discussion](https://github.com/Bg-cr/Phyth/discussions)

---

## Code Style

Phyth uses a **lightweight** set of conventions (no need to memorize a 100-page guide):

### Naming

| Type                    | Style               | Example                                            |
|-------------------------|---------------------|----------------------------------------------------|
| **Classes / Structs**   | `PascalCase`        | `Particle`, `DistanceConstrainer`, `ElectricField` |
| **Functions / Methods** | `PascalCase`        | `GetTotalEnergy()`, `ApplyForce()`, `Integrate()`  |
| **Variables**           | `snake_case`        | `delta_time`, `anchor`, `pendulum_mass`            |
| **Constants**           | `k` + `PascalCase`  | `kSpeedOfLight`, `kGravitationalConstant`          |
| **Files**               | `PascalCase.hpp`    | `Quantity.hpp`, `Spring.hpp`, `ElectricField.hpp`  |
| **Namespaces**          | `lowercase`         | `Phyth::Mechanics`, `Phyth::Electromagnetics`      |

### Headers

- Use header guard(```#ifdef #define ... #endif```) and not to use `#pragma once`.
- Order includes: **project headers → standard library headers**
- Keep headers minimal — only include what you need

```cpp
// Good
#include "Phyth/Core/Quantity.hpp"
#include "Phyth/Tools/Vector3.hpp"
#include <memory>
#include <vector>
```
```cpp
// Avoid
#include <bits/stdc++.h>  // No!
```

### Comments
- Use `//` for line comments

- Document public APIs with `///` or `/** ... */`

- Focus on why, not what (the code already says what it does)

```cpp
/// Computes gravitational force between two particles.
/// Uses Newton's law of universal gravitation: F = G * m1 * m2 / r^2
Force ComputeGravity(const Particle& a, const Particle& b);
```

---

## How to Add a New Feature
1. Open an issue first to discuss the feature (avoid wasted work)

2. Fork the repository

3. Create a branch from `main`: `git checkout -b feature/your-feature-name`

4. Implement the feature

5. Add validation tests (see below)

6. Update documentation (if needed)

7. Push and open a Pull Request

---

## Validation Requirements
**All new physics modules must include validation tests.**

This is non-negotiable — Phyth is a research-oriented library, and users need to trust the results.

### What to include
Compare against analytical solutions (preferred) or known benchmarks

Add test code to `examples/`

Document the test case and results in the PR description

### Example format
```cpp
// examples/SpringOscillation.cpp
// Validates simple harmonic motion: T = 2π * sqrt(m/k)
// Expected: T = 2.0s for m=1kg, k=π² N/m
// Measured: T = 2.00001s → relative error < 1e-5
```

### Commit Messages
Write commit messages in imperative mood:

```text
Add Coulomb's law validation test
Fix energy drift in spring constraint
Update documentation for Quantity class
```

**Format**:

- First line: ≤ 72 characters, summary of what changed

- Blank line

- (Optional) Detailed description and reasoning

```text
Add adaptive time-step integration

Current fixed-step integrator fails for stiff systems.
This adds Runge-Kutta-Fehlberg (RKF45) with automatic step-size control.
Fixes #100
```

---

## Pull Request Process
1. **PR title**: Brief summary of the change

2. **PR description**: What does this do? Why? Any breaking changes?

3. **CI**: All checks must pass (build + tests)

4. **Review**: At least one maintainer reviews the code

5. **Merge**: Squash and merge (or rebase)

### PR Checklist
□ Code follows style guidelines
□ Added validation tests (if new physics)
□ Updated documentation (if API changed)
□ All CI checks pass

---

## Getting Help
Open a [Discussion](https://github.com/Bg-cr/Phyth/discussions) for questions

Tag `@Bg-cr` in issues for maintainer attention

For urgent issues, email: `[your-email]`

---

## License
By contributing, you agree that your contributions will be licensed under the [MIT License](LICENSE).

Thank you for helping make Phyth better!