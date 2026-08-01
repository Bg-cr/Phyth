#ifndef PHYTH_MECHANICS_FORMULAS_HPP
#define PHYTH_MECHANICS_FORMULAS_HPP

#include "Phyth/Tools/FormulaSolver.hpp"
#include "Phyth/System/Config.hpp"

namespace Phyth {
    enum class GravityFormulaVariables {
        APP_POS,
        APP_MASS,
        REF_POS,
        REF_MASS,

        GRAVITY
    };

    class GravityFormula : public FormulaSolver<GravityFormulaVariables> {
    public:
        GravityFormula() {
            Set(GravityFormulaVariables::APP_POS,
                Vector3{0_m, 0_m, 0_m})
            .Set(GravityFormulaVariables::APP_MASS,
                1_kg)
            .Set(GravityFormulaVariables::REF_POS,
                Vector3{1_m, 0_m, 0_m})
            .Set(GravityFormulaVariables::REF_MASS,
                1_kg)
            .Set(GravityFormulaVariables::GRAVITY,
                1_N);
        }
    protected:
        GravityFormula &SolveImpl() override {
            const auto app_pos = Get<Vector3<Quantity<Meter> > >(GravityFormulaVariables::APP_POS);
            const auto app_mass = Get<Quantity<Kilogram> >(GravityFormulaVariables::APP_MASS);
            const auto ref_pos = Get<Vector3<Quantity<Meter> > >(GravityFormulaVariables::REF_POS);
            const auto ref_mass = Get<Quantity<Kilogram> >(GravityFormulaVariables::REF_MASS);
            const auto gravity = Get<Vector3<Quantity<Newton> > >(GravityFormulaVariables::GRAVITY);

            const auto r = ref_pos - app_pos;

            if (gravity.Length() < Quantity<Newton>(Config::epsilon)) {
                throw std::invalid_argument("Gravity is zero vector.");
            }
            if (app_mass < Quantity<Kilogram>(Config::epsilon)) {
                throw std::invalid_argument("The mass of application particle is a non-positive number.");
            }
            if (ref_mass < Quantity<Kilogram>(Config::epsilon)) {
                throw std::invalid_argument("The mass of reference particle is a non-positive number.");
            }
            if (r.Length() < Quantity<Meter>(Config::epsilon)) {
                throw std::invalid_argument("Application particle coincides with reference particle.");
            }

            const auto solved_r = gravity.Normalized() * Utils::sqrt(Consts::G * app_mass * ref_mass / gravity.Length());
            const auto solved_mass_product = r.LengthSquared() * gravity.Length() / Consts::G;

            Set(GravityFormulaVariables::APP_POS,
                ref_pos - solved_r);
            Set(GravityFormulaVariables::APP_MASS,
                solved_mass_product / ref_mass);
            Set(GravityFormulaVariables::REF_POS,
                solved_r + app_pos);
            Set(GravityFormulaVariables::REF_MASS,
                solved_mass_product / app_mass);
            Set(GravityFormulaVariables::GRAVITY,
                -Consts::G * app_mass * ref_mass / r.LengthSquared() * r.Normalized());

            return *this;
        }
    };

    enum class SpringFormulaVariables {
        APP_POS,
        REF_POS,
        STIFFNESS,
        DAMPING,
        VELOCITY,

        FORCE
    };

    class SpringFormula : public FormulaSolver<SpringFormulaVariables> {
    public:
        SpringFormula() {
            Set(SpringFormulaVariables::APP_POS, Vector3{0_m, 0_m, 0_m})
            .Set(SpringFormulaVariables::REF_POS, Vector3{1_m, 0_m, 0_m})
            .Set(SpringFormulaVariables::STIFFNESS, 1_Npm)
            .Set(SpringFormulaVariables::DAMPING, 1_kgps)
            .Set(SpringFormulaVariables::VELOCITY, Vector3{1_mps, 0_mps, 0_mps})
            .Set(SpringFormulaVariables::FORCE, Vector3{0_N, 0_N, 0_N});
        }
    protected:
        SpringFormula &SolveImpl() override {
            const auto app_pos = Get<Vector3<Quantity<Meter>>>(SpringFormulaVariables::APP_POS);
            const auto ref_pos = Get<Vector3<Quantity<Meter>>>(SpringFormulaVariables::REF_POS);
            const auto stiffness = Get<Quantity<NewtonPerMeter>>(SpringFormulaVariables::STIFFNESS);
            const auto damping = Get<Quantity<KilogramPerSecond>>(SpringFormulaVariables::DAMPING);
            const auto velocity = Get<Vector3<Quantity<MeterPerSecond>>>(SpringFormulaVariables::VELOCITY);
            const auto force = Get<Vector3<Quantity<Newton>>>(SpringFormulaVariables::FORCE);

            const auto x = app_pos - ref_pos;
            const auto solved_x = -(force + damping * velocity) / stiffness;

            Set(SpringFormulaVariables::APP_POS,
                ref_pos + solved_x);
            Set(SpringFormulaVariables::REF_POS,
                app_pos - solved_x);

            Set(SpringFormulaVariables::STIFFNESS,
                -(force + damping * velocity).Dot(x) / x.LengthSquared());
            Set(SpringFormulaVariables::DAMPING,
                -(force + stiffness * x).Dot(velocity) / velocity.LengthSquared());

            Set(SpringFormulaVariables::VELOCITY,
                -(force + stiffness * x) / damping);

            Set(SpringFormulaVariables::FORCE,
                -stiffness * x - damping * velocity);

            return *this;
        }
    };
}

#endif //PHYTH_MECHANICS_FORMULAS_HPP
