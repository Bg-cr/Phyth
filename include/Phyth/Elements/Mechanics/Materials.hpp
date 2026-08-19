#ifndef PHYTH_MATERIALS_HPP
#define PHYTH_MATERIALS_HPP

#include "Phyth/Core/Units.hpp"
#include "IsotropicMaterial.hpp"

#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>

namespace Phyth::Mechanics::Material {
    namespace detail {
        /**
         * @brief Parser for material property files in INI format
         */
        class MaterialINIParser {
        public:
            /**
             * @brief Material property data structure
             *
             * Contains the physical properties required to define
             * an isotropic linear elastic material.
             */
            struct MaterialData {
                Quantity<Pascal> youngs_modulus;    ///< Young's modulus E (Pa)
                Scalar poisson_ratio;               ///< Poisson's ratio nu (dimensionless)
            };

            /**
             * @brief Parse an INI file and extract material definitions
             *
             * @param filename Path to the INI file containing material definitions
             * @return std::unordered_map<std::string, MaterialData> Map of material names to their properties
             *
             * @throws std::runtime_error If the file cannot be opened
             * @throws std::runtime_error If a section header is malformed
             * @throws std::runtime_error If a numeric value cannot be parsed
             */
            static std::unordered_map<std::string, MaterialData> Parse(const std::string &filename) {
                std::ifstream file(filename);
                if (!file.is_open()) {
                    throw std::runtime_error("Cannot open material file: " + filename);
                }

                std::unordered_map<std::string, MaterialData> result;
                std::string line, current_section;
                MaterialData current_data{.youngs_modulus = 0_Pa, .poisson_ratio = 0_};

                while (std::getline(file, line)) {
                    trim(line);

                    if (line.empty() || line[0] == '#' || line[0] == ';') {
                        continue;
                    }

                    if (line[0] == '[') {
                        if (!current_section.empty()) {
                            result[current_section] = current_data;
                        }

                        auto end = line.find(']');
                        if (end == std::string::npos) {
                            throw std::runtime_error("Invalid section header: " + line);
                        }
                        current_section = line.substr(1, end - 1);
                        trim(current_section);
                        current_data = MaterialData{.youngs_modulus = 0.0, .poisson_ratio = 0.0};
                        continue;
                    }

                    if (auto eq_pos = line.find('='); eq_pos != std::string::npos) {
                        std::string key = line.substr(0, eq_pos);
                        std::string value = line.substr(eq_pos + 1);
                        trim(key);
                        trim(value);

                        if (auto comment_pos = value.find_first_of("#;"); comment_pos != std::string::npos) {
                            value = value.substr(0, comment_pos);
                            trim(value);
                        }

                        if (key == "youngs_modulus" || key == "E") {
                            current_data.youngs_modulus = parseNumber(value) * 1_Pa;
                        } else if (key == "poisson_ratio" || key == "nu") {
                            current_data.poisson_ratio = parseNumber(value) * 1_;
                        }
                    }
                }

                if (!current_section.empty()) {
                    result[current_section] = current_data;
                }

                return result;
            }

        private:
            /**
             * @brief Remove leading and trailing whitespace characters
             *
             * @param str String to trim (modified in-place)
             * @note Trims spaces, tabs, carriage returns, and newlines
             */
            static void trim(std::string &str) {
                const auto start = str.find_first_not_of(" \t\r\n");
                if (start == std::string::npos) {
                    str.clear();
                    return;
                }
                const auto end = str.find_last_not_of(" \t\r\n");
                str = str.substr(start, end - start + 1);
            }

            /**
             * @brief Parse a string as a double-precision floating point number
             *
             * @param str String containing a numeric value
             * @return double Parsed numeric value
             *
             * @throws std::runtime_error If the string cannot be parsed as a number
             *
             * @note Supports scientific notation (e.g., "2.1e-3")
             */
            static double parseNumber(const std::string &str) {
                try {
                    // Scientific notation
                    return std::stod(str);
                } catch (const std::exception &) {
                    throw std::runtime_error("Invalid number format: " + str);
                }
            }
        };
    }

    /**
     * @brief Database for managing and retrieving material properties
     *
     * Provides a central registry for isotropic materials, supporting
     * registration, lookup, and loading from configuration files.
     *
     * The database stores materials by name and returns IsotropicMaterial
     * objects for use in structural analysis computations.
     */
    class MaterialDatabase {
    public:
        /**
         * @brief Retrieve a material by its registered name
         *
         * @param name Name of the material to retrieve
         * @return IsotropicMaterial Material properties
         *
         * @throws std::runtime_error If material name is not found in database
         *
         * @note Case-sensitive matching
         * @warning Throws exception if material does not exist
         *
         * Example:
         *   auto material = db.Get("Aluminum");
         *   auto nu = material.GetPoissonRatio();
         */
        [[nodiscard]] IsotropicMaterial Get(const std::string &name) const {
            const auto it = materials_.find(name);
            if (it == materials_.end()) {
                throw std::runtime_error("Material not found: " + name);
            }
            return it->second;
        }

        /**
         * @brief Get a list of all registered material names
         *
         * @return std::vector<std::string> Vector containing all material names
         *
         * @note Order is not guaranteed (unordered_map iteration order)
         *
         * Example:
         *   auto names = db.ListMaterials();
         *   for (const auto& name : names) {
         *       std::cout << "Available: " << name << std::endl;
         *   }
         */
        [[nodiscard]] std::vector<std::string> ListMaterials() const {
            std::vector<std::string> names;
            for (const auto &[name, _]: materials_) {
                names.push_back(name);
            }
            return names;
        }

        /**
         * @brief Register a new material in the database
         *
         * @param name Unique identifier for the material
         * @param material IsotropicMaterial object containing properties
         *
         * @note If a material with the same name exists, it will be overwritten
         * @warning Overwriting existing materials may cause unexpected behavior
         *
         * Example:
         *   db.RegisterMaterial("CustomSteel",
         *       IsotropicMaterial(210e9_Pa, 0.3_));
         */
        void RegisterMaterial(const std::string &name, const IsotropicMaterial &material) {
            materials_[name] = material;
        }

        /**
         * @brief Create a MaterialDatabase from an INI configuration file
         *
         * Parses the specified INI file and registers all materials found.
         *
         * @param path Path to the INI configuration file
         * @return MaterialDatabase& Reference to the constructed database
         * @throws std::runtime_error If the file cannot be opened or parsed
         *
         * @par INI File Format:
         * @code
         *   [MaterialName]
         *   youngs_modulus = 210e9     # E in Pa
         *   poisson_ratio = 0.3        # nu (dimensionless)
         * @endcode
         * or
         * @code
         *   [MaterialName]
         *   E = 210e9     # E in Pa
         *   nu = 0.3        # nu (dimensionless)
         * @endcode
         *
         * Supported keys: `youngs_modulus`/`E` and `poisson_ratio`/`nu`
         * Comments start with `#` or `;`
         *
         * @par Example:
         * @code
         *   auto& db = MaterialDatabase::FromINI("standard.ini");
         *   auto steel = db.Get("Steel");
         * @endcode
         */
        [[nodiscard]] static MaterialDatabase &FromINI(const std::string &path) {
            MaterialDatabase db;
            db.materials_ = {};
            const auto materials = detail::MaterialINIParser::Parse(path);
            for (const auto &[category, material_data] : materials) {
                db.RegisterMaterial(category, IsotropicMaterial(material_data.youngs_modulus, material_data.poisson_ratio));
            }
            return db;
        }

    private:
        std::unordered_map<std::string, IsotropicMaterial> materials_;  ///< Internal storage mapping names to materials
    };

    /**
     * @brief Get the global singleton material database instance
     *
     * Returns a reference to the static MaterialDatabase instance.
     * Use this for application-wide material access and registration.
     *
     * @return MaterialDatabase& Reference to the global database
     *
     * Example:
     *   auto& db = Phyth::Mechanics::Material::GetMaterialDatabase();
     *   db.RegisterMaterial("Titanium", IsotropicMaterial(110e9_Pa, 0.34_));
     *   auto titanium = db.Get("Titanium");
     */
    inline MaterialDatabase &GetMaterialDatabase() {
        static MaterialDatabase db;
        return db;
    }
}

#endif //PHYTH_MATERIALS_HPP