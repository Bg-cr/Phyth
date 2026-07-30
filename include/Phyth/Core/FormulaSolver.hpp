#ifndef PHYTH_FORMULA_SOLVER_HPP
#define PHYTH_FORMULA_SOLVER_HPP

#include <any>
#include <map>

namespace Phyth {
    /**
     * @brief Given all variables, treat each variable as an unknown number at one time, and treat the rest as known numbers, solve and overwrite
     *
     * @tparam VariableEnum Variable enum class
     */
    template<typename VariableEnum>
    class FormulaSolver {
    public:
        explicit FormulaSolver() = default;
        virtual ~FormulaSolver() = default;

        template <typename T>
        [[nodiscard]] T Get(const VariableEnum index) const {
            auto it = variables_.find(index);
            if (it == variables_.end()) {
                throw std::runtime_error("Variable not set: " + std::to_string(static_cast<int>(index)));
            }
            try {
                return std::any_cast<T>(it->second);
            } catch (const std::bad_any_cast& e) {
                throw std::runtime_error("Type mismatch for variable: " + std::to_string(static_cast<int>(index)));
            }
        }

        template <typename T>
        FormulaSolver& Set(const VariableEnum index, T value) {
            variables_[index] = value;
            return *this;
        }

        FormulaSolver& ReloadSnapshot() {
            variables_ = prev_snapshot_;
            return *this;
        }

        FormulaSolver& DeleteVariable(const VariableEnum index) {
            variables_.erase(index);
            return *this;
        }

        FormulaSolver& ClearVariables() {
            variables_.clear();
            return *this;
        }

        FormulaSolver& Solve() {
            SaveSnapshot();
            return SolveImpl();
        }

    protected:
        std::map<VariableEnum, std::any> variables_;
        std::map<VariableEnum, std::any> prev_snapshot_;

        virtual FormulaSolver& SolveImpl() = 0;

        void SaveSnapshot() {
            prev_snapshot_ = variables_;
        }
    };
}

#endif //PHYTH_FORMULA_SOLVER_HPP
