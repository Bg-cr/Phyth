#ifndef PHYTH_CONFIG_HPP
#define PHYTH_CONFIG_HPP

namespace Phyth {
    enum class OutputMode {
        Auto,
        Raw
    };


    namespace Config {
        static inline auto output_mode = OutputMode::Auto;
        static inline double epsilon = 1e-6;

        inline void Reset() noexcept {
            output_mode = OutputMode::Auto;
            epsilon = 1e-6;
        }
    }
}
#endif  //PHYTH_CONFIG_HPP