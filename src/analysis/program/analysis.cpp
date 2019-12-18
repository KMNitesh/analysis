#include "config.h"
#include <iostream>
#include <vector>
#include <list>
#include <memory>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "utils/common.hpp"

#include "mainUtils.hpp"
#include "ana_module/IRSpectrum.hpp"
#include "ana_module/IRSpectrumDeltaDipole.hpp"
#include "others/RamanSpectrum.hpp"
#include "others/CrossCorrelation.hpp"
#include "others/GmxTopologyPrinter.hpp"
#include "others/GQuadruplexPdb2gmx.hpp"
#include "others/NBOSpin.hpp"
#include "others/Averager.hpp"
#include "others/ITS_PostProcess.hpp"
#include "others/ITS_Reweight.hpp"
#include "others/GromosReader.hpp"

using namespace std;


void printDSLDetails() {

    std::cout << "Author : " << CACANA_AUTHOR << "\n\n";

    cout <<
         "\n"
         "                         <Summary Syntax for Control Script>\n"
         "\n"
         " <-- Trajectory Analysis Commands -->"
         "\n"
         " Radical Distribution Function :\n"
         "   rdf ( M : AmberMask, L : AmberMask, max_dist : double = 10.0,  width : double = 0.01,\n"
         "                                                   intramol : bool = false, out : string )\n"
         "\n"
         " Rotational Autocorrelation Function :\n"
         "   rotacf ( vector : VectorSelector, P : int, time_increment_ps : double = 0.1,\n"
         "                                       max_time_grap_ps : double, out : string )\n"
         "\n"
         " Rotational Autocorrelation Function within Solvation Shell :\n"
         "   rotacfCutoff ( M : AmberMask, L : AmberMask, vector : VectorSelector, P : int,\n"
         "                               cutoff : double, time_increment_ps : double = 0.1,\n"
         "                                         max_time_grap_ps : double, out : string )\n"
         "\n"
         " DemixIndex of Bi-component System :\n"
         "   demix ( component1 : AmberMask, component2 : AmberMask, grid : Grid, out : string )\n"
         "\n"
         " Residence Time of Ligand in Solvation Shell :\n"
         "   residenceTime ( M : AmberMask, L : AmberMask, cutoff : double, time_star : int, out : string )\n"
         "\n"
         " Self-Diffusion Coefficient Calculation based on Einstein Equation :\n"
         "   diffuse ( mask : AmberMask, time_increment_ps : double = 0.1, total_frames : int, out : string )\n"
         "\n"
         " Self-Diffusion Coefficient Calculation based on Einstein Equation within Solvation Shell:\n"
         "   diffuseCutoff ( M : AmberMask, L : AmberMask, cutoff : double,\n"
         "                                       time_increment_ps : double = 0.1, out : string )\n"
         "\n"
         "\n"
         " <-- Auxiliary Functions -->\n"
         "     VectorSelector DipoleVector  ( mask  : AmberMask )\n"
         "     VectorSelector TwoAtomVector ( mask1 : AmberMask, mask2 : AmberMask )\n"
         "     VectorSelector NormalVector  ( mask1 : AmberMask, mask2 : AmberMask, mask3 : AmberMask )\n"
         "     Grid  grid ( x : int, y : int, z : int )\n"
         "\n"
         "     go ( start : int = 1, end : int = 0, step : int = 1, nthreads = 0 )\n"
         "     readTop ( file : string )\n"
         "     trajin  ( file : string )\n"
         "     readFF  ( file : string )\n"
         "\n"
         "\n"
         " <-- Extended Backus-Naur Form for Script -->\n"
         "\n"
         " <literal> ::= <int> | <double> | <bool> | <AmberMask> | <string>\n"
         " <for_stmt> ::= \"for\" \"(\" <expr> \";\" <condtion_expr> \";\" <expr> \")\" \"{\" <stmts> \"}\"\n"
         " <while_stmt> ::= \"while\" \"(\" <condtion_expr> \")\" \"{\" <stmts> \"}\"\n"
         " <if_else_stmt> ::= \"if\" \"(\" <condtion_expr> \")\" \"{\" <stmts> \"}\" [ \"else\" \"{\" <stmts> \"}\" ]\n"
         " <do_while_stmt> ::= \"do\" \"{\" <stmts> \"}\" \"while\" \"(\" <condtion_expr> \")\" \";\"\n"
         "\n"
         "  <- Implemented Operator Precedence ->\n"
         " Associativity ( L = Left-to-right, R = Right-to-left )\n"
         " L          a()             Function call\n"
         "---------------------------------------------------------------------------------\n"
         " L          a++  a--        Suffix/postfix increament and decrement\n"
         "---------------------------------------------------------------------------------\n"
         " R          ++a  --a        Prefix increament and decrement\n"
         " R           +a   -a        Unray plus and minus\n"
         " R            !             Logical NOT\n"
         "---------------------------------------------------------------------------------\n"
         " L          a*b  a/b  a%b   Multiplication, division, and remainder\n"
         "---------------------------------------------------------------------------------\n"
         " L          a+b  a-b        Addition and subtraction\n"
         "---------------------------------------------------------------------------------\n"
         " L           <   <=         Less, less equal relational operaters\n"
         " L           >   >=         Great, great equal relational operaters\n"
         "---------------------------------------------------------------------------------\n"
         " L          ==   !=         Equal, not equal releational operaters\n"
         "---------------------------------------------------------------------------------\n"
         " L           &              Bitwise AND\n"
         "---------------------------------------------------------------------------------\n"
         " L           |              Bitwise OR(inclusive or)\n"
         "---------------------------------------------------------------------------------\n"
         " L          &&              Logical AND\n"
         "---------------------------------------------------------------------------------\n"
         " L          ||              Logical OR\n"
         "---------------------------------------------------------------------------------\n"
         " R           =              Direct assignment\n"
         " R          +=   -=         Compound assignment by sum and difference\n"
         " R          *=   /=   %=    Compound assignment by product, quotient, and remainder\n"
         " R          &=   |=         Compound assignment by bitwise AND and OR\n"
         "\n"
         "  <- Internal Functions ->\n"
         "      int(), double(), bool(), sqrt(), abs(), log(), exp(), pow(), print()\n";

}


int main(int argc, char *argv[]) {

    std::cout << "Build DateTime : " << __DATE__ << " " << __TIME__ << endl;

    std::cout << "current work dir : " << boost::filesystem::current_path() << std::endl;

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
        printDSLDetails();
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

    /*
     *  This is the main menu the user select when the program starts
     *  evergy function of option may has its own submenu, by using different handling models
     */
    auto mainMenu = [] {
        std::cout << "Main Menu\n";
        std::cout << "(0) Trajectory Analysis\n";
        std::cout << "(1) Print Topology\n";
        std::cout << "(2) Infrared radiation (IR) Spectrum\n";
        std::cout << "(3) Infrared radiation (IR) Spectrum from DeltaDipole\n";
        std::cout << "(4) " << RamanSpectrum::title() << '\n';
        std::cout << "(5) " << CrossCorrelation::title() << '\n';
        std::cout << "(6) " << GmxTopologyPrinter::title() << '\n';
        std::cout << "(7) " << GQuadruplexPdb2gmx::title() << '\n';
        std::cout << "(8) " << "Superpose and move for Residues" << '\n';
        std::cout << "(9) " << NBOSpin::title() << '\n';
        std::cout << "(10) Renumber atom and residue num\n";
        std::cout << "(11) " << Averager::title() << '\n';
        std::cout << "(12) " << ITS_PostProcess::title() << '\n';
        std::cout << "(13) " << ITS_Reweight::title() << '\n';
        std::cout << "(14) " << GromosReader::title() << '\n';
        return choose<int>(0, 14, "select : ");
    };

    std::vector<std::function<void()>> actions{
            [&] { processTrajectory(desc, vm, xyzfiles, argc, argv); },
            [&] { printTopolgy(vm); },
            [&] { IRSpectrum::calculateSpectrum(getOutputFilename(vm)); },
            [&] { IRSpectrumDeltaDipole::calculateSpectrum(getOutputFilename(vm)); },
            [&] { RamanSpectrum::calculateSpectrum(getOutputFilename(vm)); },
            [&] { CrossCorrelation::calculate(getOutputFilename(vm)); },
            [&] { GmxTopologyPrinter::print(getTopologyFilename(vm), getPrmFilename(vm), getOutputFilename(vm)); },
            [&] { GQuadruplexPdb2gmx::convert(); },
            [&] { GQuadruplexPdb2gmx::superpose_and_move(); },
            [&] { NBOSpin::process(); },
            [&] { GQuadruplexPdb2gmx::renumberAtomAndResidueNum(); },
            [&] { Averager::process(); },
            [&] { ITS_PostProcess::process(); },
            [&] { ITS_Reweight::process(); },
            [&] { GromosReader::process(); }
    };

    actions.at(mainMenu())();

    return EXIT_SUCCESS;
}




