#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "ana_module/IRSpectrum.hpp"
#include "ana_module/IRSpectrumDeltaDipole.hpp"
#include "config.h"
#include "mainUtils.hpp"
#include "others/ADCHCharge.hpp"
#include "others/Averager.hpp"
#include "others/CrossCorrelation.hpp"
#include "others/DelocalizationIndex.hpp"
#include "others/File47CoordindateFormat.hpp"
#include "others/GQuadruplexPdb2gmx.hpp"
#include "others/GaussianLogCoordinateFormat.hpp"
#include "others/GmxTopologyPrinter.hpp"
#include "others/GroRenumber.hpp"
#include "others/GromosReader.hpp"
#include "others/HOOH_Calculator.hpp"
#include "others/ITS_PostProcess.hpp"
#include "others/ITS_Reweight.hpp"
#include "others/MultiwfnAIMDriver.hpp"
#include "others/NBOOrbitalComposition.hpp"
#include "others/NBOSpin.hpp"
#include "others/QMStructureComp.hpp"
#include "others/RamanSpectrum.hpp"
#include "others/TrajConverter.hpp"
#include "utils/ProgramConfiguration.hpp"
#include "utils/common.hpp"

void printDSLDetails() {
    std::cout << "\nAuthor : " << CACANA_AUTHOR << "\n\n";
    std::cout << SCRIPT_SYNTAX_SUMMARY;
}

int main(int argc, char *argv[]) {
    std::cout << "Build DateTime : " << __DATE__ << " " << __TIME__ << '\n';
    std::cout << "current work dir : " << boost::filesystem::current_path() << '\n';

#ifndef NDEBUG
    std::cout << " Program Arugments < \n";
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "] = " << argv[i] << '\n';
    }
    std::cout << " > Program Arugments  \n";
    std::cout << "Envirment Variables :\n";
    std::cout << "ANALYSIS_VECTOR_RESERVE = " << getDefaultVectorReserve() << '\n';
#endif

    auto desc = make_program_options();
    boost::program_options::variables_map vm;
    try {
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).run(), vm);
    } catch (std::exception &e) {
        std::cerr << e.what() << '\n';
        std::cout << desc;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help")) {
        std::cout << desc;
        printDSLDetails();
        exit(EXIT_SUCCESS);
    }

    verbose_message = vm.contains("verbose");

    program_configuration = vm.contains("config")
                                ? std::make_unique<ProgramConfiguration>(vm["config"].as<std::string>())
                                : std::make_unique<ProgramConfiguration>();

    bool keep_silent = vm["silent"].as<bool>();

    /*
     * examine all input files are already exist
     */
    std::vector<std::string> xyzfiles;
    if (vm.count("file")) {
        xyzfiles = vm["file"].as<std::vector<std::string>>();
        for (auto &xyzfile : xyzfiles) {
            if (!boost::filesystem::exists(xyzfile)) {
                std::cerr << "The file " << xyzfile << " is bad !\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    if (vm.count("target") and vm.count("script")) {
        std::cerr << "target and script option cannot coexist\n";
        exit(EXIT_FAILURE);
    }
    if (vm.count("target") and vm.count("script-file")) {
        std::cerr << "target and script-file option cannot coexist\n";
        exit(EXIT_FAILURE);
    }
    if (vm.count("script") and vm.count("script-file")) {
        std::cerr << "script and script-file option cannot coexist\n";
        exit(EXIT_FAILURE);
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
    if (vm.count("script") || vm.count("script-file")) {
        executeScript(desc, vm, xyzfiles, argc, argv);
        return EXIT_SUCCESS;
    }

    if (vm.count("aim")) {
        MultiwfnAIMDriver::process(vm["aim"].as<std::string>());
        return EXIT_SUCCESS;
    }

    if (vm.count("di")) {
        DelocalizationIndex::process();
        return EXIT_SUCCESS;
    }

    if (vm.count("adch")) {
        ADCHCharge::process();
        return EXIT_SUCCESS;
    }

    /*
     *  This is the main menu the user select when the program starts
     *  evergy function of option may has its own submenu, by using different handling models
     */

    std::vector<std::function<void()>> actions{
        [&] { processTrajectory(desc, vm, xyzfiles, argc, argv); },
        [&] { printTopolgy(vm); },
        [&] { IRSpectrum::calculateSpectrum(getOutputFilename(vm)); },
        [&] { IRSpectrumDeltaDipole::calculateSpectrum(getOutputFilename(vm)); },
        [&] { RamanSpectrum::calculateSpectrum(getOutputFilename(vm)); },
        [&] { CrossCorrelation::calculate(getOutputFilename(vm)); },
        [&] { GmxTopologyPrinter::print(getTopologyFilename(vm), getPrmFilename(vm), getOutputFilename(vm)); },
        [] { GQuadruplexPdb2gmx::convert(); },
        [] { GQuadruplexPdb2gmx::superpose_and_move(); },
        [] { NBOSpin::process(); },
        [] { GQuadruplexPdb2gmx::renumberAtomAndResidueNum(); },
        [] { Averager::process(); },
        [] { ITS_PostProcess::process(); },
        [] { ITS_Reweight::process(); },
        [] { GromosReader::process(); },
        [] { MultiwfnAIMDriver::process_interactive(); },
        [] { NBOOrbitalComposition::process(); },
        [] { DelocalizationIndex::process_interactive(); },
        [&] {
            ADCHCharge::process_interactive(vm.count("fchk") ? vm["fchk"].as<std::string>()
                                                             : boost::optional<std::string>{});
        },
        [] { TrajConverter::process(); },
        [] { QMStructureComp::process(); },
        [] { File47CoordindateFormat::process(); },
        [] { GaussianLogCoordinateFormat::process(); },
        [] { HOOH_Calculator::process(); },
        [] { GroRenumber::process(); }};

    auto mainMenu = [&] {
        if (!keep_silent) {
            std::cout << "Main Menu\n";
            std::cout << " (0) Trajectory Analysis\n";
            std::cout << " (1) Print Topology\n";
            std::cout << " (2) Infrared radiation (IR) Spectrum\n";
            std::cout << " (3) Infrared radiation (IR) Spectrum from DeltaDipole\n";
            std::cout << " (4) " << RamanSpectrum::title() << '\n';
            std::cout << " (5) " << CrossCorrelation::title() << '\n';
            std::cout << " (6) " << GmxTopologyPrinter::title() << '\n';
            std::cout << " (7) " << GQuadruplexPdb2gmx::title() << '\n';
            std::cout << " (8) "
                      << "Superpose and move for Residues" << '\n';
            std::cout << " (9) " << NBOSpin::title() << '\n';
            std::cout << "(10) Renumber atom and residue num\n";
            std::cout << "(11) " << Averager::title() << '\n';
            std::cout << "(12) " << ITS_PostProcess::title() << '\n';
            std::cout << "(13) " << ITS_Reweight::title() << '\n';
            std::cout << "(14) " << GromosReader::title() << '\n';
            std::cout << "(15) " << MultiwfnAIMDriver::title() << '\n';
            std::cout << "(16) " << NBOOrbitalComposition::title() << '\n';
            std::cout << "(17) " << DelocalizationIndex::title() << '\n';
            std::cout << "(18) " << ADCHCharge::title() << '\n';
            std::cout << "(19) " << TrajConverter::title() << '\n';
            std::cout << "(20) " << QMStructureComp::title() << '\n';
            std::cout << "(21) " << File47CoordindateFormat::title() << '\n';
            std::cout << "(22) " << GaussianLogCoordinateFormat::title() << '\n';
            std::cout << "(23) " << HOOH_Calculator::title() << '\n';
            std::cout << "(24) " << GroRenumber::title() << '\n';
        }
        return choose<int>(0, actions.size() - 1, "select : ");
    };

    actions.at(mainMenu())();

    return EXIT_SUCCESS;
}
