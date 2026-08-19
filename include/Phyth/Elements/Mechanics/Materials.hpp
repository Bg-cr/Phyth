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
        class MaterialINIParser {
        public:
            struct MaterialData {
                Quantity<Pascal> youngs_modulus;
                Scalar poisson_ratio;
            };

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
            static void trim(std::string &str) {
                const auto start = str.find_first_not_of(" \t\r\n");
                if (start == std::string::npos) {
                    str.clear();
                    return;
                }
                const auto end = str.find_last_not_of(" \t\r\n");
                str = str.substr(start, end - start + 1);
            }

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

    class MaterialDatabase {
    public:
        [[nodiscard]] IsotropicMaterial Get(const std::string &name) const {
            const auto it = materials_.find(name);
            if (it == materials_.end()) {
                throw std::runtime_error("Material not found: " + name);
            }
            return it->second;
        }

        [[nodiscard]] std::vector<std::string> ListMaterials() const {
            std::vector<std::string> names;
            for (const auto &[name, _]: materials_) {
                names.push_back(name);
            }
            return names;
        }

        void RegisterMaterial(const std::string &name, const IsotropicMaterial &material) {
            materials_[name] = material;
        }

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
        std::unordered_map<std::string, IsotropicMaterial> materials_;
    };

    inline MaterialDatabase &GetMaterialDatabase() {
        static MaterialDatabase db;
        return db;
    }
}

#endif //PHYTH_MATERIALS_HPP
