
#include "BondEnergyCompare.hpp"
#include "utils/BondEnergyCalculator.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/range/algorithm.hpp>

#include "utils/common.hpp"
#include <boost/spirit/include/qi.hpp>

namespace {

auto parse(const std::string &prompt) {
    for (;;) {
        boost::fusion::vector<std::vector<boost::fusion::vector<int, boost::optional<int>>>,
                              std::vector<boost::fusion::vector<int, boost::optional<int>>>>
            ast;
        auto line = input(prompt);
        namespace qi = boost::spirit::qi;
        if (auto it = std::begin(line); qi::phrase_parse(it, std::end(line),
                                                         ((qi::int_ >> -('-' >> qi::int_)) % ',') >> ':' >>
                                                             ((qi::int_ >> -('-' >> qi::int_)) % ',') >> qi::eoi,
                                                         qi::ascii::space, ast) &&
                                        it == std::end(line)) {
            return ast;
        }
        std::cerr << "Error re type !!!!\n";
    }
}

bool isIn(std::vector<boost::fusion::vector<int, boost::optional<int>>> &ast, int index) {
    for (const auto &item : ast) {
        const auto &first = boost::fusion::at_c<0>(item);
        const auto &second = boost::fusion::at_c<1>(item);
        if (second.has_value()) {
            if (index >= first and index <= second.get())
                return true;
        } else {
            if (index == first)
                return true;
        }
    }
    return false;
}
} // namespace

void BondEnergyCompare::process() {
    std::string filename_1 = choose_file("dump1 : ").extension("json").can_empty(false).isExist(true);
    std::string filename_2 = choose_file("dump2 : ").extension("json").can_empty(false).isExist(true);

    auto index_map = parse("Index map > ");

    std::ifstream ifs1(filename_1), ifs2(filename_2);

    nlohmann::json j1, j2;
    ifs1 >> j1;
    ifs2 >> j2;

    const auto data1 = j1["data"];
    const auto data2 = j2["data"];
    std::vector<std::pair<std::pair<int, int>, double>> diff;
    std::map<std::pair<int, int>, boost::accumulators::accumulator_set<
                                      double, boost::accumulators::features<boost::accumulators::tag::mean>>>
        diff_map;

    for (std::size_t i = 0; i < std::min(data1.size(), data2.size()); ++i) {
        const auto &field1 = data1[i];
        const auto &field2 = data2[i];
        auto it1 = begin(field1);
        auto it2 = begin(field2);
        for (; it1 != end(field1) and it2 != end(field2);) {
            const auto &v1 = it1->get<std::pair<int, BondEnergyCalculator::Term>>();
            const auto &v2 = it2->get<std::pair<int, BondEnergyCalculator::Term>>();
            if (!isIn(boost::fusion::at_c<0>(index_map), v1.first)) {
                ++it1;
                continue;
            }
            if (!isIn(boost::fusion::at_c<1>(index_map), v2.first)) {
                ++it2;
                continue;
            }
            diff_map[{v1.first, v2.first}](v1.second.total() - v2.second.total());
            ++it1;
            ++it2;
        }
    }
    diff.reserve(diff_map.size());
    for (const auto &[key, value] : diff_map)
        diff.emplace_back(key, boost::accumulators::mean(value));

    boost::sort(diff, [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    std::string output_file = choose_file("output file > ");

    std::ofstream ofs(output_file);

    for (const auto &[key, value] : diff) {
        ofs << boost::format("%10d %10d %20.6f\n") % key.first % key.second % value;
    }
}