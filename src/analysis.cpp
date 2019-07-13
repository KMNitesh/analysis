#include <iostream>
#include <vector>
#include <list>
#include <memory>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "common.hpp"

#include "mainUtils.hpp"

using namespace std;

/*
 *  This is the main menu the user select when the program starts
 *  evergy function of option may has its own submenu, by using different handling models
 */

int mainMenu() {
    std::cout << "Main Menu\n";
    std::cout << "(0) Trajectory Analysis\n";
    std::cout << "(1) Print Topology\n";
    return choose<int>(0, 1, "select :");
};


int main(int argc, char *argv[]) {

    std::cout << "Build DateTime : " << __DATE__ << " " << __TIME__ << endl;

    std::cout << "current work dir : " << boost::filesystem::current_path() << std::endl;

    auto desc = make_program_options();
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    } catch (std::exception &e) {
        std::cerr << e.what() << '\n';
        std::cout << desc;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help")) {
        cout << desc;
        exit(EXIT_SUCCESS);
    }

    /*
     * examine all input files are already exist
     */
    std::vector<std::string> xyzfiles;
    if (vm.count("file")) {
        xyzfiles = vm["file"].as<std::vector<string>>();
        for (auto &xyzfile : xyzfiles) {
            if (!boost::filesystem::exists(xyzfile)) {
                cerr << "The file " << xyzfile << " is bad !" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    /*
     *  when --target or -x option given,  fast trajectory convert mode is active
     *  Convert trajctory format for convenience
     */
    if (vm.count("target")) {
        fastTrajectoryConvert(vm, xyzfiles);
        exit(EXIT_SUCCESS);
    }


    /*
     * Enter script execution, for non-interactive mode
     */
    if (vm.count("script")) {
        executeScript(desc, vm, xyzfiles, argc, argv);
        return EXIT_SUCCESS;
    }


    /*
     *  examine the tpr file of gromcas and show the content
     */
    if (mainMenu() == 1) {
        printTopolgy(vm);
        return EXIT_SUCCESS;
    }


    /*
     *  At last, handle normal trajctories
     */
    processTrajectory(desc, vm, xyzfiles, argc, argv);

    return EXIT_SUCCESS;

}



