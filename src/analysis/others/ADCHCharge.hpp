#ifndef TINKER_ADCHCHARGE_HPP
#define TINKER_ADCHCHARGE_HPP

#include <boost/fusion/sequence.hpp>
#include <string_view>
#include <boost/optional.hpp>

class ADCHCharge {
public:
    [[nodiscard]] static std::string_view title() {
        return "Atomic dipole corrected Hirshfeld population (ADCH) Charge";
    }

    static void process_interactive(boost::optional<std::string> file);

    static void process();

    static void process(const std::string &file, const std::vector<int> &atoms);

    [[nodiscard]] static boost::optional<boost::fusion::vector<std::string, double, double>> read_charge(
        const std::string &line);
};

#endif  // TINKER_ADCHCHARGE_HPP
