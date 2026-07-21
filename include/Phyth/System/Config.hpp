#ifndef PHYTH_CONFIG_H
#define PHYTH_CONFIG_H

namespace Phyth {
    enum class OutputMode {
        Auto,
        Raw
    };

    class Config {
    public:
        static void setOutputMode(const OutputMode mode) {
            output_mode_ = mode;
        }

        static OutputMode getOutputMode() {
            return output_mode_;
        }

        static void reset() {
            output_mode_ = OutputMode::Auto;
        }

    private:
        static OutputMode output_mode_;
        static bool show_unit_space_;
    };

    inline OutputMode Config::output_mode_ = OutputMode::Auto;
    inline bool Config::show_unit_space_ = true;

    inline void setOutputMode(const OutputMode mode) {
        Config::setOutputMode(mode);
    }

}

#endif  //PHYTH_CONFIG_H