#include <boost/spirit/include/qi.hpp>
#include "ITS_Reweight.hpp"
#include "utils/common.hpp"
#include "dsl/Interpreter.hpp"

void ITS_Reweight::process() {

    std::string fb_filename = choose_file("Enter fb.dat file > ").isExist(true).extension("dat");
    std::ifstream fb_ifs(fb_filename);
    if (!fb_ifs) {
        std::cerr << "ERROR! Cannot open fb file <" << fb_filename << ">\n";
        exit(EXIT_FAILURE);
    }

    std::string pot_filename = choose_file("Enter pot file > ").isExist(true).extension("dat");

    std::ifstream pot_ifs(pot_filename);
    if (!pot_ifs) {
        std::cerr << "ERROR! Cannot open pot file <" << pot_filename << ">\n";
        exit(EXIT_FAILURE);
    }

    std::string factor_filename = choose_file("Enter factor output file > ").isExist(false);

    std::ofstream factor_os(factor_filename);
    if (!factor_os) {
        std::cerr << "ERROR! Cannot open factor output file <" << factor_filename << ">\n";
        exit(EXIT_FAILURE);
    }

    auto[ok, fb_data] = read_fb(fb_ifs);
    if (!ok) {
        std::cerr << "ERROR! the fb file is ill-formed !\n";
        exit(EXIT_FAILURE);
    }

    auto[is_ok, pot_data] = read_pot(pot_ifs);
    if (!is_ok) {
        std::cerr << "ERROR! the pot file is ill-formed !\n";
        exit(EXIT_FAILURE);
    }

    int nb = choose(1, "Enter nb > ");
    double lowT = choose(0.0, "Enter lowT > ");
    double highT = choose(0.0, "Enter highT > ");
    double vshift = choose(0.0, "Enter vshift(kJ/mol) > ");
    double objtemp = choose(0.0, "Enter objtemp > ");
    double integrate_step = choose(0.0, "Enter integrate_step > ");

    std::vector<double> mybeta(nb);

    constexpr auto coff = 1000 / (6.022 * 1.380653);
    const auto beta0 = coff / objtemp;

    const auto dtemp = (highT - lowT) / (nb - 1);

    for (int i = 0; i < nb; ++i) {
        mybeta[i] = coff / (lowT + i * dtemp);
    }

    std::vector<double> gf(nb);

    double gfsum;
    auto it = std::begin(fb_data);
    for (auto[time, v] : pot_data) {
        while (std::abs(boost::fusion::at_c<0>(*it) * integrate_step - time) > std::numeric_limits<double>::epsilon()) {
            ++it;
        }
        auto &fb = boost::fusion::at_c<1>(*it);
        auto vb = v - vshift;
        for (int i = 0; i < nb; ++i) {
            gf[i] = -mybeta[i] * vb + fb[i];
        }

        gfsum = gf[0];
        for (int i = 1; i < nb; ++i) {
            if (gfsum > gf[i])
                gfsum = gfsum + std::log(1.0 + std::exp(gf[i] - gfsum));
            else
                gfsum = gf[i] + std::log(1.0 + std::exp(gfsum - gf[i]));
        }

        auto factor = std::exp(-beta0 * v - gfsum);

        factor_os << std::setw(10) << std::fixed << time
                  << std::setw(10) << v
                  << std::setw(15) << std::scientific << factor
                  << '\n';
    }

}

std::pair<bool, std::vector<boost::fusion::vector<int, std::vector<double>>>> ITS_Reweight::read_fb(std::istream &is) {

    is.unsetf(std::ios::skipws);
    std::vector<boost::fusion::vector<int, std::vector<double>>> fb_data;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;
    auto parser = (int_ >> ':' >> +(double_ >> ',')) % eol >> *eol;

    boost::spirit::istream_iterator begin(is);
    boost::spirit::istream_iterator end;

    auto match = phrase_parse(begin, end, parser, ascii::space - eol, fb_data);

    return {match && begin == end, fb_data};
}

std::pair<bool, std::vector<std::pair<double, double>>> ITS_Reweight::read_pot(std::istream &is) {

    is.unsetf(std::ios::skipws);
    std::vector<std::pair<double, double>> sequences;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;


    qi::rule<boost::spirit::istream_iterator, std::pair<double, double>,
            decltype(SkipperGrammar<boost::spirit::istream_iterator>() - eol)>
            item = (double_ >> double_)[_val = construct<std::pair<double, double>>(_1, _2)];

    auto parser = (item % eol) >> *eol;

    boost::spirit::istream_iterator begin(is);
    boost::spirit::istream_iterator end;

    auto match = phrase_parse(begin, end, parser, SkipperGrammar<boost::spirit::istream_iterator>() - eol, sequences);

    return {match && begin == end, sequences};
}
