#include <boost/spirit/include/qi.hpp>
#include "ITS_PostProcess.hpp"
#include "utils/common.hpp"

void ITS_PostProcess::process() {
    std::string input_file = choose_file("Input dat file > ").isExist(true);
    std::string output_file = choose_file("Output file > ").isExist(false);

    std::ifstream ifs(input_file);
    if (!ifs) {
        std::cerr << "ERROR !! Cannot open input file " << input_file << '\n';
        exit(EXIT_FAILURE);
    }

    std::ofstream ofs(output_file);
    if (!ofs) {
        std::cerr << "ERROR !! Cannot open output file " << output_file << '\n';
        exit(EXIT_FAILURE);
    }

    bool change_log_base = choose_bool("Change log base > ");

    int step;
    std::vector<double> values;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;
    auto parser = int_[ref(step) = _1] >> ':' >> +(double_[push_back(boost::phoenix::ref(values), _1)] >> ',');

    std::string line;

    while (!ifs.eof()) {
        std::getline(ifs, line);
        if (line.empty()) continue;
        auto it = begin(line);
        if (!phrase_parse(it, end(line), parser, ascii::space) || it != end(line)) {
            std::cerr << "Input file content is ill-formed\n";
            exit(EXIT_FAILURE);
        }
        print(ofs, change_log_base, step, values);
        values.clear();
    }
}

void ITS_PostProcess::print(std::ofstream &ofs, bool change_log_base, int step, const std::vector<double> &values) {
    ofs << std::setw(15) << step;
    for (auto &v : values) ofs << std::setw(15) << (change_log_base ? (v / log_base) : v);
    ofs << '\n';
}


