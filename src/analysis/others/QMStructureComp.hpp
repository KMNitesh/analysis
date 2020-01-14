#ifndef TINKER_QMSTRUCTURECOMP_HPP
#define TINKER_QMSTRUCTURECOMP_HPP

#include <string_view>
#include <boost/fusion/sequence.hpp>

class QMStructureComp {
public:
    [[nodiscard]] static std::string_view title() { return "QM Structure Comparison"; }

    static void process();

    [[nodiscard]] static std::vector<boost::fusion::vector<uint, double, double, double>>
    read_47_file(std::istream &is);

    [[nodiscard]] static std::vector<boost::fusion::vector<std::string, double, double, double>>
    read_log_file(std::istream &is);
};


#endif //TINKER_QMSTRUCTURECOMP_HPP
