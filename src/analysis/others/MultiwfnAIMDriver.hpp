#ifndef TINKER_MULTIWFNAIMDRIVER_HPP
#define TINKER_MULTIWFNAIMDRIVER_HPP

#include <tuple>

class MultiwfnAIMDriver {
public:

    [[nodiscard]] static std::string_view title() { return "QTAIM analysis of BCP"; }

    static void process();

    static void process_interactive();

private:

    struct BCP_property {
        double Density_of_all_electrons;
        double Hb;
        double Laplacian_of_electron_density;
        double Ellipticity;
    };


    struct BCP {
        BCP(int i, int j) : atom_pair_belonging{i, j} {}

        int index{};
        std::tuple<double, double, double> coord;
        std::pair<int, int> atom_pair_belonging;

        BCP_property property{};
    };

    static void readBCP(std::istream &is, std::vector<BCP> &bcp_vector);

    static void print(const std::vector<BCP> &bcp_vector);

    [[nodiscard]] static std::tuple<std::string, int, std::vector<int>> inputParameter();

    static void printProperty(const boost::format &fmt, const std::vector<BCP> &bcp_vector,
                              std::string_view name, double BCP_property::* field);

    static void process(const std::string &file, const std::vector<std::pair<int, int>> &bonds);

};


#endif //TINKER_MULTIWFNAIMDRIVER_HPP
