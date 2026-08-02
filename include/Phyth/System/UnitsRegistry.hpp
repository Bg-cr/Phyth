#ifndef PHYTH_UNIT_REGISTRY_HPP
#define PHYTH_UNIT_REGISTRY_HPP

#include <string>
#include <unordered_map>
#include <typeindex>
#include <mutex>

namespace Phyth {

    struct UnitInfo {
        std::string name;
        std::string symbol;
    };

    /**
     * @brief Singleton registry for unit names and symbols
     *
     * Stores display information for unit types. Used by QuantityToString()
     * to print human-readable unit names instead of dimension strings.
     */
    class UnitRegistry {
    public:
        static UnitRegistry& GetInstance() {
            static UnitRegistry registry;
            return registry;
        }

        template<typename UnitT>
        void RegisterUnit(const std::string& name, const std::string& symbol) {
            std::lock_guard lock(mutex_);
            registry_[std::type_index(typeid(UnitT))] = {.name = name, .symbol = symbol};
        }

        template<typename UnitT>
        std::string GetName() const {
            if (const auto it = registry_.find(std::type_index(typeid(UnitT))); it != registry_.end()) {
                return it->second.name;
            }
            return "";
        }

        template<typename UnitT>
        std::string GetSymbol() const {
            if (const auto it = registry_.find(std::type_index(typeid(UnitT))); it != registry_.end()) {
                return it->second.symbol;
            }
            return "";
        }

    private:
        UnitRegistry() = default;
        std::unordered_map<std::type_index, UnitInfo> registry_;
        mutable std::mutex mutex_;
    };

    /**
     * @brief Registers a unit with the registry at static initialization time
     *
     * @tparam UnitT The unit type to register
     */
    template<typename UnitT>
    struct UnitRegistrar {
        UnitRegistrar(const std::string& name, const std::string& symbol) {
            UnitRegistry::GetInstance().RegisterUnit<UnitT>(name, symbol);
        }
    };
}

#endif  //PHYTH_UNIT_REGISTRY_HPP