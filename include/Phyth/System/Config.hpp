#ifndef PHYTH_CONFIG_HPP
#define PHYTH_CONFIG_HPP

namespace Phyth {
    enum class OutputMode {
        Auto,
        Raw
    };


    namespace Config {
        inline auto output_mode = OutputMode::Auto;
        inline double epsilon = 1e-6;
        inline int max_iterations = 50;
        inline double tolerance = 1e-12;

        inline void Reset() noexcept {
            output_mode = OutputMode::Auto;
            epsilon = 1e-6;
            max_iterations = 50;
            tolerance = 1e-12;
        }
    }
}
#endif  //PHYTH_CONFIG_HPP