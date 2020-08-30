#include "config.h"

#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/type_index.hpp>

#include "ana_module/IRSpectrum.hpp"
#include "ana_module/IRSpectrumDeltaDipole.hpp"
#include "mainUtils.hpp"
#include "others/ADCHCharge.hpp"
#include "others/Averager.hpp"
#include "others/BondEnergyCompare.hpp"
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
#include "others/MMPBSA.hpp"
#include "others/MultiwfnAIMDriver.hpp"
#include "others/NBOOrbitalComposition.hpp"
#include "others/NBOSpin.hpp"
#include "others/QMStructureComp.hpp"
#include "others/RamanSpectrum.hpp"
#include "others/TrajConverter.hpp"
#include "others/US_Pull_Wins.hpp"
#include "utils/ProgramConfiguration.hpp"
#include "utils/common.hpp"

void printDSLDetails() {
    std::cout << "\nAuthor : " << CACANA_AUTHOR << "\n\n";
    std::cout << SCRIPT_SYNTAX_SUMMARY;
}

int main(int argc, char *argv[]) {
    try {
        std::cout << "Build DateTime : " << __DATE__ << " " << __TIME__ << '\n';
        std::cout << "current work dir : " << boost::filesystem::current_path() << '\n';

#ifndef NDEBUG
        std::cout << " Program Arguments < \n";
        for (int i = 0; i < argc; i++) {
            std::cout << "argv[" << i << "] = " << argv[i] << '\n';
        }
        std::cout << " > Program Arguments  \n";
        std::cout << "Environment Variables :\n";
        std::cout << "ANALYSIS_VECTOR_RESERVE = " << getDefaultVectorReserve() << '\n';
#endif

        auto desc = make_program_options();
        boost::program_options::variables_map vm;
        try {
            boost::program_options::store(
                boost::program_options::command_line_parser(argc, argv).options(desc).run(), vm);
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
        debug_mode = vm.contains("debug");

        program_configuration =
            vm.contains("config")
                ? std::make_unique<ProgramConfiguration>(vm["config"].as<std::string>())
                : std::make_unique<ProgramConfiguration>();

        bool keep_silent = vm["silent"].as<bool>();

        /*
         * examine all input files are already exist
         */
        std::vector<std::string> xyzfiles;
        if (vm.contains("file")) {
            xyzfiles = vm["file"].as<std::vector<std::string>>();
            for (auto &xyzfile : xyzfiles) {
                if (!boost::filesystem::exists(xyzfile)) {
                    std::cerr << "The file " << xyzfile << " is not exist!\n";
                    exit(EXIT_FAILURE);
                }
            }
        }

        if (vm.contains("target") and vm.contains("script")) {
            std::cerr << "target and script option cannot coexist\n";
            exit(EXIT_FAILURE);
        }
        if (vm.contains("target") and vm.contains("script-file")) {
            std::cerr << "target and script-file option cannot coexist\n";
            exit(EXIT_FAILURE);
        }
        if (vm.contains("script") and vm.contains("script-file")) {
            std::cerr << "script and script-file option cannot coexist\n";
            exit(EXIT_FAILURE);
        }
        /*
         *  when --target or -x option given,  fast trajectory convert mode is active
         *  Convert trajctory format for convenience
         */
        if (vm.contains("target")) {
            fastTrajectoryConvert(vm, xyzfiles);
            exit(EXIT_SUCCESS);
        }

        /*
         * Enter script execution, for non-interactive mode
         */
        if (vm.contains("script") || vm.contains("script-file")) {
            executeScript(desc, vm, xyzfiles, argc, argv);
            return EXIT_SUCCESS;
        }

        if (vm.contains("aim")) {
            MultiwfnAIMDriver::process(vm["aim"].as<std::string>());
            return EXIT_SUCCESS;
        }

        if (vm.contains("di")) {
            DelocalizationIndex::process();
            return EXIT_SUCCESS;
        }

        if (vm.contains("adch")) {
            ADCHCharge::process();
            return EXIT_SUCCESS;
        }

        /*
         *  This is the main menu the user select when the program starts
         *  evergy function of option may has its own submenu, by using different handling models
         */

        struct Action {
            Action(std::string_view name, std::function<void()> action)
                : name(name), action(std::move(action)) {}

            std::string name;
            std::function<void()> action;
        };
        Action actions[]{
            {"Trajectory Analysis", [&] { processTrajectory(desc, vm, xyzfiles, argc, argv); }},
            {"Print Topology", [&] { printTopolgy(vm); }},
            {"Infrared radiation (IR) Spectrum",
             [&] { IRSpectrum::calculateSpectrum(getOutputFilename(vm)); }},
            {"Infrared radiation (IR) Spectrum from DeltaDipole",
             [&] { IRSpectrumDeltaDipole::calculateSpectrum(getOutputFilename(vm)); }},
            {RamanSpectrum::title(),
             [&] { RamanSpectrum::calculateSpectrum(getOutputFilename(vm)); }},
            {CrossCorrelation::title(),
             [&] { CrossCorrelation::calculate(getOutputFilename(vm)); }},
            {GmxTopologyPrinter::title(),
             [&] {
                 GmxTopologyPrinter::print(getTopologyFilename(vm), getPrmFilename(vm),
                                           getOutputFilename(vm));
             }},
            {GQuadruplexPdb2gmx::title(), [] { GQuadruplexPdb2gmx::convert(); }},
            {"Superpose and move for Residues", [] { GQuadruplexPdb2gmx::superpose_and_move(); }},
            {NBOSpin::title(), [] { NBOSpin::process(); }},
            {"Renumber atom and residue num",
             [] { GQuadruplexPdb2gmx::renumberAtomAndResidueNum(); }},
            {Averager::title(), [] { Averager::process(); }},
            {ITS_PostProcess::title(), [] { ITS_PostProcess::process(); }},
            {ITS_Reweight::title(), [] { ITS_Reweight::process(); }},
            {GromosReader::title(), [] { GromosReader::process(); }},
            {MultiwfnAIMDriver::title(), [] { MultiwfnAIMDriver::process_interactive(); }},
            {NBOOrbitalComposition::title(), [] { NBOOrbitalComposition::process(); }},
            {DelocalizationIndex::title(), [] { DelocalizationIndex::process_interactive(); }},
            {ADCHCharge::title(),
             [&] {
                 ADCHCharge::process_interactive(vm.count("fchk") ? vm["fchk"].as<std::string>()
                                                                  : boost::optional<std::string>{});
             }},
            {TrajConverter::title(), [] { TrajConverter::process(); }},
            {QMStructureComp::title(), [] { QMStructureComp::process(); }},
            {File47CoordindateFormat::title(), [] { File47CoordindateFormat::process(); }},
            {GaussianLogCoordinateFormat::title(), [] { GaussianLogCoordinateFormat::process(); }},
            {HOOH_Calculator::title(), [] { HOOH_Calculator::process(); }},
            {GroRenumber::title(), [] { GroRenumber::process(); }},
            {BondEnergyCompare::title(), [] { BondEnergyCompare::process(); }},
            {MMPBSA::title(), [&] { MMPBSA::process(getTopologyFilename(vm)); }},
            {US_Pull_Wins::title(), [] { US_Pull_Wins::process(); }}};

        auto mainMenu = [&] {
            if (!keep_silent) {
                std::cout << "Main Menu\n";
                for (const auto &element : actions | boost::adaptors::indexed())
                    std::cout << "(" << element.index() << ") " << element.value().name << '\n';
            }
            return choose<int>(0, std::extent_v<decltype(actions)> - 1, "select : ");
        };

        actions[mainMenu()].action();

    } catch (std::exception &e) {
        std::cerr << "Exception(" << boost::typeindex::type_id_runtime(e).pretty_name()
                  << ") : " << e.what() << std::endl;
    } catch (...) {
        boost::typeindex::stl_type_index sti = *std::current_exception().__cxa_exception_type();
        std::cerr << "Exception(" << sti.pretty_name() << ")\n";
    }

    return EXIT_SUCCESS;
}
