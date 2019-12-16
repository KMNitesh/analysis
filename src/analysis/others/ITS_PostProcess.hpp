#ifndef TINKER_ITS_POSTPROCESS_HPP
#define TINKER_ITS_POSTPROCESS_HPP

#include <string_view>

class ITS_PostProcess {
public:
    [[nodiscard]] static std::string_view title() { return "ITS Postprocess"; }

    static void process();

    static void print(std::ofstream &ofs, bool change_log_base, int step, const std::vector<double> &values);

    inline static const auto log_base = log(10.0);
};


#endif //TINKER_ITS_POSTPROCESS_HPP
