#ifndef TINKER_ITS_REWEIGHT_HPP
#define TINKER_ITS_REWEIGHT_HPP

#include <string_view>

class ITS_Reweight {
public:
    [[nodiscard]] static std::string_view title() { return "ITS Reweight"; }

    static void process();

    [[nodiscard]] static std::pair<bool, std::vector<boost::fusion::vector<int, std::vector<double>>>>
    read_fb(std::istream &is);

    [[nodiscard]] static std::pair<bool, std::vector<std::pair<double, double>>> read_pot(std::istream &is);

};


#endif //TINKER_ITS_REWEIGHT_HPP
