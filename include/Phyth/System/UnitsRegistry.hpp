#ifndef PHYTH_UNIT_REGISTRY_HPP
#define PHYTH_UNIT_REGISTRY_HPP

#include <string>
#include <unordered_map>
#include <typeindex>
#include <mutex>

namespace Phyth {
    /**
     * @brief Display information for a unit type
     *
     * Contains human-readable strings used by QuantityToString().
     */
    struct UnitInfo {
        std::string name; /**< Full unit name, e.g., "Meter", "Kilogram" */
        std::string symbol; /**< Unit symbol, e.g., "m", "kg", "N" */
    };

    /**
     * @brief Singleton registry for unit names and symbols
     *
     * UnitRegistry stores display information for unit types. It is used by
     * QuantityToString() to print human-readable unit names instead of raw
     * dimension strings (e.g., "J" instead of "m^2*kg*s^-2").
     *
     * The registry is thread-safe via a mutex.
     *
     * Usage:
     *   // Register a unit manually
     *   UnitRegistry::GetInstance().RegisterUnit<Meter>("Meter", "m");
     *
     *   // Register a unit at static initialization time using UnitRegistrar
     *   static UnitRegistrar<Meter> meterRegistrar("Meter", "m");
     *
     *   // Retrieve registered info
     *   auto name = UnitRegistry::GetInstance().GetName<Meter>();   // "Meter"
     *   auto symbol = UnitRegistry::GetInstance().GetSymbol<Meter>(); // "m"
     */
    class UnitRegistry {
    public:
        /**
         * @brief Get the singleton instance
         *
         * @return Reference to the single UnitRegistry instance
         */
        static UnitRegistry &GetInstance() {
            static UnitRegistry registry;
            return registry;
        }

        /**
         * @brief Register a unit type with its display information
         *
         * @tparam UnitT The unit type to register
         * @param name Full unit name (e.g., "Meter")
         * @param symbol Unit symbol (e.g., "m")
         *
         * This function is thread-safe.
         */
        template<typename UnitT>
        void RegisterUnit(const std::string &name, const std::string &symbol) {
            std::lock_guard lock(mutex_);
            registry_[std::type_index(typeid(UnitT))] = {.name = name, .symbol = symbol};
        }

        /**
         * @brief Get the full name of a unit type
         *
         * @tparam UnitT The unit type to look up
         * @return The unit name, or empty string if not registered
         *
         * This function is thread-safe (const method with mutable mutex).
         */
        template<typename UnitT>
        std::string GetName() const {
            if (const auto it = registry_.find(std::type_index(typeid(UnitT))); it != registry_.end()) {
                return it->second.name;
            }
            return "";
        }

        /**
         * @brief Get the symbol of a unit type
         *
         * @tparam UnitT The unit type to look up
         * @return The unit symbol, or empty string if not registered
         *
         * This function is thread-safe (const method with mutable mutex).
         */
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
        mutable std::mutex mutex_; // mutable to allow locking in const methods
    };

    /**
     * @brief Static registrar for unit types
     *
     * UnitRegistrar registers a unit at static initialization time, before main()
     * is called. This allows units to be registered automatically without requiring
     * explicit registration calls.
     *
     * @tparam UnitT The unit type to register
     *
     * Usage:
     *   // In a .cpp file or header (with inline):
     *   static UnitRegistrar<Meter> meterRegistrar("Meter", "m");
     *   static UnitRegistrar<Second> secondRegistrar("Second", "s");
     *   static UnitRegistrar<Newton> newtonRegistrar("Newton", "N");
     */
    template<typename UnitT>
    struct UnitRegistrar {
        /**
         * @brief Register the unit with the given name and symbol
         *
         * @param name Full unit name
         * @param symbol Unit symbol
         */
        UnitRegistrar(const std::string &name, const std::string &symbol) {
            UnitRegistry::GetInstance().RegisterUnit<UnitT>(name, symbol);
        }
    };
}

#endif // PHYTH_UNIT_REGISTRY_HPP
