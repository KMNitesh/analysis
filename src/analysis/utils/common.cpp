//
// Created by xiamr on 3/17/19.
//

#include "common.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <regex>
#include <utility>

#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"

bool enable_read_velocity = false;
bool enable_tbb = false;
bool enable_outfile = false;

Forcefield forcefield;
bool enable_forcefield = false;

std::vector<std::string> split(const std::string &str, const std::string &sep) {
    std::vector<std::string> ret_;
    boost::split(ret_, str, boost::is_any_of(sep));
    return ret_;
}

std::vector<std::string> split(const std::string &str) {
    std::vector<std::string> ret_;
    boost::regex e("\\s+");
    std::string s = str;
    boost::regex_split(std::back_inserter(ret_), s, e);
    return ret_;
}

std::vector<std::string> split_quoted(const std::string &str) {
    std::vector<std::string> ret_;
    std::istringstream iss(str);
    std::string s;
    while (iss >> std::quoted(s)) {
        ret_.push_back(s);
    }
    return ret_;
}

std::string input(const std::string &prompt, std::istream &in, std::ostream &out) {
    out << prompt;
    std::string inputline;
    std::getline(in, inputline);
    if (isatty(STDIN_FILENO) == 0) {
        out << inputline << std::endl;
    }
    return inputline;
}

std::string ext_filename(const std::string &filename) {
    auto field = split(filename, ".");
    assert(field.size() > 1);

    auto ext = boost::to_lower_copy(field[field.size() - 1]);

    assert(!ext.empty());
    return ext;
}

FileType getFileType(const std::string &filename) {
    auto extension = boost::filesystem::extension(filename);
    boost::to_lower(extension);

    const std::unordered_map<std::string, FileType> mapping = {
        {".xtc", FileType::XTC},      {".trr", FileType::TRR},   {".nc", FileType::NC},   {".mdcrd", FileType::NC},
        {".xyz", FileType::ARC},      {".arc", FileType::ARC},   {".tpr", FileType::TPR}, {".prmtop", FileType::PRMTOP},
        {".parm7", FileType::PRMTOP}, {".mol2", FileType::MOL2}, {".prm", FileType::PRM}, {".gro", FileType::GRO},
        {".traj", FileType::TRAJ},    {".json", FileType::JSON}};

    auto it = mapping.find(extension);
    return it != mapping.end() ? it->second : FileType::UnKnown;
}

bool choose_bool(const std::string &prompt, Default<bool> defaultValue, std::istream &in, std::ostream &out) {
    while (true) {
        std::string input_line = input(prompt, in, out);
        boost::trim(input_line);
        if (input_line.empty()) {
            if (defaultValue) return defaultValue.getValue();
        }
        boost::to_lower(input_line);

        if (input_line == "y") {
            return true;
        } else if (input_line == "n") {
            return false;
        } else {
            out << "Input Error, must be either y or n !!\n";
        }
    }
}

std::string choose_file(const std::string &prompt, bool exist, std::string ext, bool can_empty, std::istream &in,
                        std::ostream &out) {
    return choose_file(prompt, in, out).isExist(exist).extension(std::move(ext)).can_empty(can_empty);
}

double atom_distance(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2,
                     const std::shared_ptr<Frame> &frame) {
    return std::sqrt(atom_distance2(atom1, atom2, frame));
}

double atom_distance2(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2,
                      const std::shared_ptr<Frame> &frame) {
    auto xr = atom1->x - atom2->x;
    auto yr = atom1->y - atom2->y;
    auto zr = atom1->z - atom2->z;
    frame->image(xr, yr, zr);
    return xr * xr + yr * yr + zr * zr;
}

boost::program_options::options_description make_program_options() {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "show this help message")(
        "topology,p", po::value<std::string>()->value_name("topology-file-name"), "topology file")(
        "file,f", po::value<std::vector<std::string>>()->multitoken()->composing()->value_name("trajectory-file-name"),
        "trajectory file")("output,o", po::value<std::string>()->value_name("output-file-name"), "output file")(
        "prm", po::value<std::string>()->value_name("tinker-prm-file-name"), "force field file")(
        "target,x", po::value<std::string>()->value_name("trajectout-file-name"), "target trajectory file")(
        "script", po::value<std::string>()->value_name("script-content"), "script command for non-interactive use")(
        "script-file", po::value<std::string>()->value_name("script-file-name"), "read command from script file")(
        "aim", po::value<std::string>()->value_name("options"), "QTAIM analysis based on Multiwfn")(
        "di", "Delocalization Index based on Multiwfn")("adch",
                                                        "Atomic dipole corrected Hirshfeld population (ADCH) Charge")(
        "fchk", po::value<std::string>()->value_name("file-name"), "Gaussian Fchk File");

    return desc;
}

std::string print_cmdline(int argc, const char *const argv[]) {
    std::string cmdline;

    const char *const *p = argv;

    while (argc-- > 0) {
        if (p != argv) {
            cmdline += " ";
        }
        cmdline += *p;
        p++;
    }
    return cmdline;
}

std::string getOutputFilename(const boost::program_options::variables_map &vm) {
    return vm.count("output") ? vm["output"].as<std::string>() : choose_file("output file :").isExist(false);
}

std::string getTopologyFilename(const boost::program_options::variables_map &vm) {
    return vm.count("topology") ? vm["topology"].as<std::string>() : choose_file("topology file :").isExist(true);
}

std::string getTrajectoryFilename(const boost::program_options::variables_map &vm) {
    return vm.count("file") ? vm["file"].as<std::string>() : choose_file("trajectory file :").isExist(true);
}

std::string getPrmFilename(const boost::program_options::variables_map &vm) {
    return vm.count("prm") ? vm["prm"].as<std::string>() : choose_file("Tinker prm file :").isExist(true);
}

std::size_t getDefaultVectorReserve() {
    auto p = std::getenv("ANALYSIS_VECTOR_RESERVE");
    return p ? std::stoi(p) : 100000;
}
