#ifndef PHYTH_CONFIG_HPP
#define PHYTH_CONFIG_HPP

namespace Phyth {
    enum class OutputMode {
        Auto,
        Raw
    };


    struct Config {
        static void Reset() {
            output_mode = OutputMode::Auto;
            epsilon = 1e-6;
        }

        static OutputMode output_mode;
        static long double epsilon;

    };

    inline OutputMode Config::output_mode = OutputMode::Auto;
    inline long double Config::epsilon = 1e-6;

}

#endif  //PHYTH_CONFIG_HPP