#include <boost/fusion/include/at_c.hpp>
#include "File47CoordindateFormat.hpp"
#include "QMStructureComp.hpp"
#include "utils/common.hpp"

void File47CoordindateFormat::process() {

    std::string _47file = choose_file("Enter .47 file > ").isExist(true).extension("47");
    std::ifstream ifs_47file{_47file};

    auto coord = QMStructureComp::read_47_file(ifs_47file);

    namespace bf = boost::fusion;
    std::cout << std::fixed << std::setprecision(6);
    for (const auto &line : coord) {
        std::cout << std::setw(2) << QMStructureComp::get_element_name(bf::at_c<0>(line))
                  << std::setw(12) << bf::at_c<1>(line)
                  << std::setw(12) << bf::at_c<2>(line)
                  << std::setw(12) << bf::at_c<3>(line) << '\n';
    }

}
