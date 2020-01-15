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

    [[nodiscard]] static const std::string &get_element_name(unsigned int i) { return element_table.at(i - 1); }

private:

    inline static const std::vector<std::string> element_table{
            "H", "He",
            "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
            "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
            "Xe",
            "Cs", "Ba",
            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
            "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
            "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
            "Rf"
    };
};


#endif //TINKER_QMSTRUCTURECOMP_HPP
