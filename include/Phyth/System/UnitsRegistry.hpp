#ifndef PHYTH_UNIT_REGISTRY_H
#define PHYTH_UNIT_REGISTRY_H

#include <string>
#include <unordered_map>
#include <typeindex>
#include <mutex>

namespace Phyth {

    struct UnitInfo {
        std::string name;
        std::string symbol;
    };

    class UnitRegistry {
    public:
        static UnitRegistry& instance() {
            static UnitRegistry registry;
            return registry;
        }

        template<typename UnitT>
        void registerUnit(const std::string& name, const std::string& symbol) {
            std::lock_guard lock(mutex_);
            registry_[std::type_index(typeid(UnitT))] = {.name = name, .symbol = symbol};
        }

        template<typename UnitT>
        std::string name() const {
            if (const auto it = registry_.find(std::type_index(typeid(UnitT))); it != registry_.end()) {
                return it->second.name;
            }
            return "";
        }

        template<typename UnitT>
        std::string symbol() const {
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

    template<typename UnitT>
    struct UnitRegistrar {
        UnitRegistrar(const std::string& name, const std::string& symbol) {
            UnitRegistry::instance().registerUnit<UnitT>(name, symbol);
        }
    };

}

#endif  //PHYTH_UNIT_REGISTRY_H