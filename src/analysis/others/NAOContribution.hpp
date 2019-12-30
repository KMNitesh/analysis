#ifndef TINKER_NAOCONTRIBUTION_HPP
#define TINKER_NAOCONTRIBUTION_HPP

#include <string_view>

class NAOContribution {
public:
    [[nodiscard]] static std::string_view title() { return "NAO contribution for atom group"; }

    static void process();

    [[nodiscard]] static boost::optional<boost::fusion::vector<int, double> > parseLine(const std::string &line);

    [[nodiscard]] static std::map<int, double, std::greater<>>
    read_contributions(const std::vector<int> &attrs, std::istream &is, int orbital_number);

    static void
    print_contributions(std::string_view descriptions, std::map<int, double, std::greater<>> &contributions);
};


#endif //TINKER_NAOCONTRIBUTION_HPP
