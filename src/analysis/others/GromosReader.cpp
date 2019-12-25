#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/xpressive/xpressive_static.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/irange.hpp>
#include <boost/spirit/include/karma.hpp>
#include <pugixml.hpp>
#include "config.h"
#include "GromosReader.hpp"
#include "utils/common.hpp"

namespace bf = boost::fusion;

namespace {
    BOOST_PHOENIX_ADAPT_FUNCTION(void, trim, boost::trim, 1)

    struct Validator : boost::static_visitor<bool> {
        Validator(std::size_t maxTotalEnergyGroupNum, std::size_t maxEnergyGroupNum)
                : max_total_energy_group_num(maxTotalEnergyGroupNum), max_energy_group_num(maxEnergyGroupNum) {}

        bool operator()(const bf::vector<char, std::pair<uint, uint>> &i) const {
            auto g1 = bf::at_c<1>(i).first;
            auto g2 = bf::at_c<1>(i).second;

            return g1 >= 1 and g1 <= max_energy_group_num and g2 >= 1 and g2 <= max_energy_group_num;
        }

        bool operator()(const uint &i) const {
            return i >= 1 and i < max_total_energy_group_num;
        }

    private:
        std::size_t max_total_energy_group_num;
        std::size_t max_energy_group_num;
    };

    struct TitlePrinter : boost::static_visitor<std::string> {

        TitlePrinter(const std::vector<std::string> &menuStrings,
                     const std::array<std::string, 3> &energyNames,
                     const std::vector<std::string> &groupNames) :
                menuStrings(menuStrings),
                energy_names(energyNames),
                group_names(groupNames) {}

        std::string operator()(const bf::vector<char, std::pair<uint, uint>> &i) const {
            auto g1 = bf::at_c<1>(i).first;
            auto g2 = bf::at_c<1>(i).second;

            return energy_names[bf::at_c<0>(i) - 'a'] + "(" + group_names[g1 - 1] + "," + group_names[g2 - 1] + ")";
        }

        std::string operator()(const uint &i) const {
            return menuStrings[i - 1];
        }

    private:
        const std::vector<std::string> &menuStrings;
        const std::array<std::string, 3> &energy_names;
        const std::vector<std::string> &group_names;
    };

    struct EnergyPrinter : boost::static_visitor<double> {

        EnergyPrinter(const GromosReader::Energy &e, size_t groupNums) : e(e), group_nums(groupNums) {}

        double operator()(const bf::vector<char, std::pair<uint, uint>> &i) const {
            auto g1 = bf::at_c<1>(i).first;
            auto g2 = bf::at_c<1>(i).second;

            auto index = ((group_nums + 1) * group_nums - (group_nums + 2 - g1) * (group_nums + 1 - g1)) / 2 + g2 - 1;
            return e.intergroups[bf::at_c<0>(i) - 'a'][index];
        }

        double operator()(const uint &i) const {
            return e.energies[i - 1];
        }

    private:
        const GromosReader::Energy &e;
        std::size_t group_nums;
    };
}

void GromosReader::process() {

    std::string filename = choose_file("GROMOS omd/xml file > ").isExist(true);

    Bundle bundle;

    auto extension = ext_filename(filename);
    if (extension == "xml") {
        bundle = readXml(filename);
    } else if (extension == "omd") {
        bundle = readOmd(filename);
    } else {
        throw std::runtime_error("Error file extension [must be one of xml/omd]");
    }

    auto&[energies, menuStrings, energy_names, group_names] = bundle;

    if (extension == "omd" and choose_bool("Save to xml ? ", Default(false))) {
        std::ofstream ofstream;
        for (;;) {
            std::string filename = choose_file("xml filename > ").isExist(false).extension("xml");
            ofstream.open(filename);
            if (!ofstream) {
                std::cerr << "Can not open file <" << filename << ">\n";
                continue;
            }
            break;
        }
        save2Xml(energies, menuStrings, energy_names, group_names, ofstream);
        std::exit(EXIT_SUCCESS);
    }
    printEnergies(energies, menuStrings, energy_names, group_names);

}

void GromosReader::printMenu(const std::vector<std::string> &menuStrings, std::size_t width) {
    std::cout << ">>>  MENU  <<<<\n" << std::left;
    for (const auto &element : menuStrings | boost::adaptors::indexed(1)) {
        std::cout << '(' << element.index() << ") " << std::setw(width + 2) << element.value();
        if (element.index() % 4 == 0) std::cout << '\n';
    }
    if (menuStrings.size() % 4 != 0) std::cout << '\n';
}

void GromosReader::save2Xml(const std::vector<Energy> &energies,
                            const std::vector<std::string> &menuStrings,
                            const std::array<std::string, 3> &energy_names,
                            const std::vector<std::string> &group_names,
                            std::ostream &os) {

    pugi::xml_document doc;

    auto decl = doc.prepend_child(pugi::node_declaration);
    decl.append_attribute("version") = "1.0";
    decl.append_attribute("encoding") = "UTF-8";

    doc.append_child(pugi::node_comment).set_value("Generated by CAC-ANA [" CACANA_AUTHOR "]");

    auto finish_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    auto finish_time_format = std::put_time(std::localtime(&finish_time), "%F %T");

    std::stringstream ss;
    ss << "Data Time : " << finish_time_format;
    doc.append_child(pugi::node_comment).set_value(ss.str().c_str());

    auto meta_data_node = doc.append_child("meta_data");

    auto attribute_node = meta_data_node.append_child("attribute");
    attribute_node.append_child("step");
    attribute_node.append_child("time").append_attribute("unit") = "ps";

    auto total_energy_node = meta_data_node.append_child("total_energy");

    for (auto &energy_label : menuStrings) {
        total_energy_node.append_child("energy").append_attribute("name") = energy_label.c_str();
    }

    auto energy_group_node = meta_data_node.append_child("energy_group");

    for (auto &name : energy_names) {
        auto node = energy_group_node.append_child("energy");
        node.append_attribute("name") = name.c_str();
        for (auto it1 = group_names.begin(); it1 != group_names.end(); ++it1) {
            for (auto it2 = it1; it2 != group_names.end(); ++it2) {
                auto component = node.append_child("component");
                component.append_attribute("group1") = it1->c_str();
                component.append_attribute("group2") = it2->c_str();
            }
        }
    }

    auto data_node = doc.append_child("data");

    using namespace boost::spirit::karma;
    auto serializer = double_ % ',';

    for (auto &e : energies) {
        auto frame_node = data_node.append_child("frame");
        frame_node.append_attribute("step") = e.step;
        frame_node.append_attribute("time") = e.time;
        std::string generated;
        std::back_insert_iterator<std::string> sink{generated};
        generate(sink, serializer, e.energies);
        for (auto &eg : e.intergroups) {
            sink = ',';
            generate(sink, serializer, eg);
        }
        frame_node.append_child(pugi::node_pcdata).set_value(generated.c_str());
    }

    doc.save(os);
}

GromosReader::Bundle GromosReader::readXml(std::istream &is) {
    pugi::xml_document doc;
    if (auto result = doc.load(is); !result) {
        throw std::runtime_error("Can not open xml stream");
    }

    std::vector<Energy> energies;
    std::vector<std::string> menuStrings;
    std::array<std::string, 3> energy_names;
    std::vector<std::string> group_names;

    auto meta_data_node = doc.child("meta_data");

    auto total_energy_node = meta_data_node.child("total_energy");

    for (auto &node : total_energy_node.children("energy")) {
        menuStrings.emplace_back(node.attribute("name").value());
    }
    auto energy_group_node = meta_data_node.child("energy_group");

    for (const auto &element :  energy_group_node.children("energy") | boost::adaptors::indexed()) {
        auto index = element.index();
        auto &node = element.value();

        energy_names[index] = node.attribute("name").value();

        for (auto &component : node.children("component")) {
            std::string group1_name = component.attribute("group1").value();

            if (boost::range::find(group_names, group1_name) == std::end(group_names)) {
                group_names.emplace_back(std::move(group1_name));
            }
        }
    }

    auto data_node = doc.child("data");

    for (auto &node : data_node.children("frame")) {
        Energy e;
        e.step = node.attribute("step").as_uint();
        e.time = node.attribute("time").as_double();

        std::string value_string = node.text().as_string();
        std::vector<std::string> values = split(value_string, ",");

        auto it = std::begin(values);
        for ([[maybe_unused]] auto i: boost::irange(menuStrings.size())) {
            e.energies.push_back(std::stod(*it));
            ++it;
        }

        auto left = (values.size() - menuStrings.size()) / 3;

        for (auto &array : e.intergroups) {
            for ([[maybe_unused]] auto i: boost::irange(left)) {
                array.push_back(std::stod(*it));
                ++it;
            }
        }
        energies.emplace_back(std::move(e));
    }

    return {energies, menuStrings, energy_names, group_names};
}

GromosReader::Bundle GromosReader::readOmd(const std::string &filename) {

    // mmap all content from file
    boost::iostreams::mapped_file_source file;
    file.open(filename);
    if (!file) {
        throw std::runtime_error("ERROR! GROMOS omd file cannot open");
    }

    std::vector<Energy> energies;

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    auto timestep_parser = "TIMESTEP" >> eol >>
                                      uint_ >> double_ >> eol
                                      >> "END" >> eol;

    auto energies_parser = "ENERGIES" >> eol
                                      >> +(as_string[lexeme[+(char_ - ':')]][trim(_1)]
                                              >> ':' >> double_ >> eol)
                                      >> eol;

    auto parser = timestep_parser >> energies_parser;

    auto energy_group_parser =
            repeat(3)[as_string[lexeme[+(char_ - char_("0-9"))]][trim(_1)]
                    >> +as_string[lexeme[+char_("0-9") >> char_('-') >> +char_("0-9")]]
                    >> eol >> +(omit[lexeme[+char_("0-9") >> char_('-') >> +char_("0-9")]] >> +double_ >> eol)
                    >> eol];

    using namespace boost::xpressive;
    auto rex = (s1 = "TIMESTEP" >> -*_ >> _n >> _n)
            >> (s2 = -*_ >> _n >> _n)
            >> (s3 = boost::xpressive::repeat<3, 3>(-*_ >> _n >> _n));

    bool bFilled = false;

    std::vector<std::string> menuStrings;


    std::array<std::string, 3> energy_names;
    std::vector<std::string> group_names;

    for (cregex_iterator pos(file.begin(), file.end(), rex), end; pos != end; ++pos) {

        bf::vector<uint, double, std::vector<bf::vector<std::string, double>>> attribute;
        if (auto it = (*pos)[1].first;
                !(phrase_parse(it, (*pos)[1].second, parser, ascii::space - boost::spirit::eol, attribute) &&
                  it == (*pos)[1].second)) {
            throw std::runtime_error("omd file is ill-formed !");
        }

        std::vector<bf::vector<std::string, std::vector<std::string>, std::vector<double>>> energy_group_attribute;
        if (auto it = (*pos)[3].first;
                !(phrase_parse(it, (*pos)[3].second, energy_group_parser, ascii::space - boost::spirit::eol,
                               energy_group_attribute) && it == (*pos)[3].second)) {
            throw std::runtime_error("omd file is ill-formed !");
        }

        Energy e;
        e.step = bf::at_c<0>(attribute);
        e.time = bf::at_c<1>(attribute);

        for (auto &bpv : bf::at_c<2>(attribute)) {
            e.energies.push_back(bf::at_c<1>(bpv));
            if (!bFilled) {
                menuStrings.push_back(bf::at_c<0>(bpv));
            }
        }
        for (const auto &element : energy_group_attribute | boost::adaptors::indexed()) {
            if (!bFilled) {
                energy_names[element.index()] = bf::at_c<0>(element.value());
                group_names = bf::at_c<1>(element.value());
            }
            e.intergroups[element.index()] = bf::at_c<2>(element.value());
        }

        bFilled = true;
        energies.push_back(std::move(e));
    }
    return {energies, menuStrings, energy_names, group_names};
}

void GromosReader::printEnergies(const std::vector<Energy> &energies,
                                 const std::vector<std::string> &menuStrings,
                                 const std::array<std::string, 3> &energy_names,
                                 const std::vector<std::string> &group_names) {

    std::size_t max_length = boost::accumulate(menuStrings, std::size_t(0),
                                               [](auto init, const auto &s) { return std::max(init, s.size()); });

    printMenu(menuStrings, max_length);

    std::cout << " >>>--- Energy Group ---<<<< \n";
    std::cout << "Component : (a) " << energy_names[0] << "  (b) " << energy_names[1] << "  (c) " << energy_names[2]
              << '\n';
    std::cout << "Group : " << std::left;
    for (const auto &element : group_names | boost::adaptors::indexed(1)) {
        std::cout << '[' << element.index() << "] " << element.value();
        if (static_cast<std::size_t>(element.index()) != group_names.size()) std::cout << "  ";
    }
    std::cout << '\n';

    using namespace boost::spirit::qi;
    using namespace boost::phoenix;

    boost::spirit::qi::rule<std::string::iterator, std::pair<uint, uint>(), ascii::space_type>
            pair_parser_ = (uint_ >> ',' >> uint_)[_val = construct<std::pair<uint, uint>>(_1, _2)];

    boost::spirit::qi::rule<std::string::iterator,
            bf::vector<char, std::pair<uint, uint>>(), ascii::space_type>
            sel_paser_ = char_("a-c") >> '(' >> pair_parser_ >> ')';

    boost::spirit::qi::rule<std::string::iterator,
            boost::variant<uint, bf::vector<char, std::pair<uint, uint>>>(), ascii::space_type>
            item_parser_ = uint_ | sel_paser_;

    boost::spirit::qi::rule<std::string::iterator,
            std::vector<boost::variant<uint, bf::vector<char, std::pair<uint, uint>>>>(), ascii::space_type>
            input_parser = +item_parser_;

    for (;;) {

        std::cout << "Select Items > ";
        std::string line;
        std::getline(std::cin, line);

        std::vector<boost::variant<uint, bf::vector<char, std::pair<uint, uint>>>> attribute;

        if (auto it = begin(line);
                !(phrase_parse(it, end(line), input_parser, ascii::space, attribute) && it == end(line))) {
            std::cerr << "Paser Error !\n" << line << '\n';
            for (auto iter = line.begin(); iter != it; ++iter) std::cout << " ";
            std::cout << "^\n";
            continue;
        }

        Validator validator(menuStrings.size(), group_names.size());

        if (!boost::accumulate(attribute, true, [&validator](bool init, auto &v) {
            return init && boost::apply_visitor(validator, v);
        })) {
            std::cerr << "number out of range\n";
            continue;
        }

        std::ofstream ofs;
        for (;;) {
            std::string outfile = choose_file("Output xvg file > ").isExist(false);
            ofs.open(outfile);
            if (ofs) break;
            std::cerr << "Cannot open file <" << outfile << ">\n";
        }

        ofs << '#' << std::setw(9) << "Time";
        TitlePrinter titlePrinter(menuStrings, energy_names, group_names);
        boost::for_each(attribute, [&](auto &i) { ofs << std::setw(15) << boost::apply_visitor(titlePrinter, i); });
        ofs << '\n';

        for (auto &e : energies) {
            ofs << std::setw(10) << e.time;
            EnergyPrinter energyPrinter(e, group_names.size());
            boost::for_each(attribute,
                            [&](auto &i) { ofs << std::setw(15) << boost::apply_visitor(energyPrinter, i); });
            ofs << '\n';
        }
        break;
    }

}

GromosReader::Bundle GromosReader::readXml(const std::string &filename) {
    std::ifstream ifstream(filename);
    return readXml(ifstream);
}
