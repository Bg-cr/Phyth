#ifndef PHYTH_MATERIALS_HPP
#define PHYTH_MATERIALS_HPP

#include "Phyth/Core/Units.hpp"
#include "IsotropicMaterial.hpp"

#include <unordered_map>
#include <string>
#include <vector>

namespace Phyth::Mechanics::Material {
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

    private:
        std::unordered_map<std::string, IsotropicMaterial> materials_ = {
            {"Steel", IsotropicMaterial(205e9_Pa, 0.305_)},
            {"SteelStainless304", IsotropicMaterial(193e9_Pa, 0.285_)},
            {"SteelStainless316", IsotropicMaterial(193e9_Pa, 0.30_)},
            {"SteelTool", IsotropicMaterial(210e9_Pa, 0.292_)},
            {"Aluminum6061", IsotropicMaterial(69.5e9_Pa, 0.33_)},
            {"Aluminum7075", IsotropicMaterial(71.7e9_Pa, 0.33_)},
            {"Aluminum2024", IsotropicMaterial(73.1e9_Pa, 0.33_)},
            {"AluminumCast", IsotropicMaterial(70e9_Pa, 0.33_)},
            {"Copper", IsotropicMaterial(120e9_Pa, 0.34_)},
            {"CopperBeryllium", IsotropicMaterial(130e9_Pa, 0.30_)},
            {"TitaniumGrade2", IsotropicMaterial(105e9_Pa, 0.32_)},
            {"TitaniumGrade5", IsotropicMaterial(113e9_Pa, 0.325_)},
            {"TitaniumGrade9", IsotropicMaterial(105e9_Pa, 0.33_)},
            {"Magnesium", IsotropicMaterial(45e9_Pa, 0.30_)},
            {"MagnesiumAZ31", IsotropicMaterial(45e9_Pa, 0.35_)},
            {"Nickel", IsotropicMaterial(200e9_Pa, 0.31_)},
            {"Inconel718", IsotropicMaterial(207e9_Pa, 0.297_)},
            {"Inconel625", IsotropicMaterial(208e9_Pa, 0.31_)},
            {"HastelloyC276", IsotropicMaterial(205e9_Pa, 0.30_)},
            {"Monel400", IsotropicMaterial(179e9_Pa, 0.32_)},
            {"Brass", IsotropicMaterial(100e9_Pa, 0.35_)},
            {"Bronze", IsotropicMaterial(120e9_Pa, 0.34_)},
            {"PhosphorBronze", IsotropicMaterial(110e9_Pa, 0.35_)},
            {"Beryllium", IsotropicMaterial(300e9_Pa, 0.08_)},
            {"Tungsten", IsotropicMaterial(400e9_Pa, 0.28_)},
            {"Molybdenum", IsotropicMaterial(325e9_Pa, 0.30_)},
            {"Zinc", IsotropicMaterial(70e9_Pa, 0.25_)},
            {"Lead", IsotropicMaterial(16e9_Pa, 0.44_)},
            {"Tin", IsotropicMaterial(50e9_Pa, 0.36_)},
            {"Gold", IsotropicMaterial(78e9_Pa, 0.42_)},
            {"Silver", IsotropicMaterial(76e9_Pa, 0.37_)},
            {"Platinum", IsotropicMaterial(168e9_Pa, 0.38_)},
            {"Palladium", IsotropicMaterial(121e9_Pa, 0.39_)},
            {"PZT5H", IsotropicMaterial(63e9_Pa, 0.315_)},
            {"PZT5A", IsotropicMaterial(61e9_Pa, 0.325_)},
            {"PZT4", IsotropicMaterial(79e9_Pa, 0.325_)},
            {"PZT8", IsotropicMaterial(78e9_Pa, 0.31_)},
            {"PZT5J", IsotropicMaterial(50e9_Pa, 0.31_)},
            {"BaTiO3", IsotropicMaterial(223e9_Pa, 0.22_)},
            {"PMNPT", IsotropicMaterial(72.5e9_Pa, 0.32_)},
            {"PZNPT", IsotropicMaterial(75e9_Pa, 0.31_)},
            {"PVDF", IsotropicMaterial(1.8e9_Pa, 0.325_)},
            {"Alumina", IsotropicMaterial(355e9_Pa, 0.215_)},
            {"Alumina96", IsotropicMaterial(326e9_Pa, 0.22_)},
            {"Alumina99", IsotropicMaterial(367e9_Pa, 0.245_)},
            {"Zirconia", IsotropicMaterial(205e9_Pa, 0.25_)},
            {"ZirconiaTZP", IsotropicMaterial(214e9_Pa, 0.275_)},
            {"SiliconCarbide", IsotropicMaterial(450e9_Pa, 0.16_)},
            {"SiliconNitride", IsotropicMaterial(309e9_Pa, 0.275_)},
            {"BoronCarbide", IsotropicMaterial(450e9_Pa, 0.17_)},
            {"BoronNitride", IsotropicMaterial(132e9_Pa, 0.205_)},
            {"TitaniumCarbide", IsotropicMaterial(450e9_Pa, 0.185_)},
            {"TungstenCarbide", IsotropicMaterial(505e9_Pa, 0.22_)},
            {"SodaLimeGlass", IsotropicMaterial(71e9_Pa, 0.22_)},
            {"BorosilicateGlass", IsotropicMaterial(63.5e9_Pa, 0.20_)},
            {"FusedSilica", IsotropicMaterial(73e9_Pa, 0.16_)},
            {"Diamond", IsotropicMaterial(1060e9_Pa, 0.14_)},
            {"Silicon", IsotropicMaterial(149e9_Pa, 0.285_)},
            {"Germanium", IsotropicMaterial(112e9_Pa, 0.25_)},
            {"GalliumArsenide", IsotropicMaterial(96e9_Pa, 0.275_)},
            {"IndiumPhosphide", IsotropicMaterial(75e9_Pa, 0.26_)},
            {"TitaniumDioxide", IsotropicMaterial(219e9_Pa, 0.285_)},
            {"YAG", IsotropicMaterial(305e9_Pa, 0.275_)},
            {"Sapphire", IsotropicMaterial(390e9_Pa, 0.27_)},
            {"Spinels", IsotropicMaterial(200e9_Pa, 0.255_)},
            {"Pyrex", IsotropicMaterial(62.8e9_Pa, 0.20_)},
            {"VYCOR", IsotropicMaterial(67.6e9_Pa, 0.19_)},
            {"Rubber", IsotropicMaterial(0.016e9_Pa, 0.49_)},
            {"Silicone", IsotropicMaterial(0.091e9_Pa, 0.445_)},
            {"PDMS", IsotropicMaterial(0.00137e9_Pa, 0.475_)},
            {"Neoprene", IsotropicMaterial(0.08e9_Pa, 0.47_)},
            {"Viton", IsotropicMaterial(0.0066e9_Pa, 0.49_)},
            {"BunaN", IsotropicMaterial(0.007e9_Pa, 0.49_)},
            {"Butyl", IsotropicMaterial(0.2e9_Pa, 0.45_)},
            {"EPDM", IsotropicMaterial(0.006e9_Pa, 0.48_)},
            {"Polyurethane", IsotropicMaterial(1.83e9_Pa, 0.325_)},
            {"PMMA", IsotropicMaterial(3.97e9_Pa, 0.37_)},
            {"PTFE", IsotropicMaterial(0.50e9_Pa, 0.455_)},
            {"Nylon6", IsotropicMaterial(2.8e9_Pa, 0.39_)},
            {"Nylon66", IsotropicMaterial(3.3e9_Pa, 0.41_)},
            {"ABS", IsotropicMaterial(2.2e9_Pa, 0.38_)},
            {"Epoxy", IsotropicMaterial(3.2e9_Pa, 0.375_)},
            {"Polycarbonate", IsotropicMaterial(2.4e9_Pa, 0.37_)},
            {"Polypropylene", IsotropicMaterial(1.85e9_Pa, 0.415_)},
            {"Polyethylene", IsotropicMaterial(1.9e9_Pa, 0.39_)},
            {"Acrylic", IsotropicMaterial(2.95e9_Pa, 0.37_)},
            {"Delrin", IsotropicMaterial(3.0e9_Pa, 0.40_)},
            {"PET", IsotropicMaterial(3.7e9_Pa, 0.405_)},
            {"PEEK", IsotropicMaterial(3.65e9_Pa, 0.39_)},
            {"QuartzIsoApprox", IsotropicMaterial(72e9_Pa, 0.17_)},
            {"LiNbO3IsoApprox", IsotropicMaterial(120e9_Pa, 0.30_)},
            {"LiTaO3IsoApprox", IsotropicMaterial(110e9_Pa, 0.31_)},
            {"KNbO3IsoApprox", IsotropicMaterial(90e9_Pa, 0.30_)},
            {"NaNbO3IsoApprox", IsotropicMaterial(85e9_Pa, 0.31_)},
            {"BiFeO3IsoApprox", IsotropicMaterial(100e9_Pa, 0.28_)},
            {"AlNIsoApprox", IsotropicMaterial(310e9_Pa, 0.24_)},
            {"ZnOIsoApprox", IsotropicMaterial(140e9_Pa, 0.30_)},
            {"BoneIsoApprox", IsotropicMaterial(18e9_Pa, 0.30_)}
        };
    };

    inline MaterialDatabase &GetMaterialDatabase() {
        static MaterialDatabase db;
        return db;
    }
}

#endif //PHYTH_MATERIALS_HPP
