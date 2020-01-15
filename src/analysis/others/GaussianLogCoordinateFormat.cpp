#include <boost/fusion/include/at_c.hpp>
#include "GaussianLogCoordinateFormat.hpp"
#include "QMStructureComp.hpp"
#include "utils/common.hpp"

void GaussianLogCoordinateFormat::process() {

    std::string log_filename = choose_file("Enter gaussian log file > ").isExist(true).extension("log");
    std::ifstream ifs_logfile{log_filename};

    auto coord = QMStructureComp::read_log_file(ifs_logfile);

    namespace bf = boost::fusion;
    std::cout << std::fixed << std::setprecision(10);
    for (const auto &line : coord) {
        std::cout << std::setw(2) << bf::at_c<0>(line)
                  << std::setw(15) << bf::at_c<1>(line)
                  << std::setw(15) << bf::at_c<2>(line)
                  << std::setw(15) << bf::at_c<3>(line) << '\n';
    }
}
