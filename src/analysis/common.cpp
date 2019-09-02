//
// Created by xiamr on 3/17/19.
//


#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <regex>
#include <boost/filesystem.hpp>
#include "common.hpp"
#include "frame.hpp"
#include "forcefield.hpp"


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
            {".xtc",   FileType::XTC},
            {".trr",   FileType::TRR},
            {".nc",    FileType::NC},
            {".mdcrd", FileType::NC},
            {".xyz",   FileType::ARC},
            {".arc",   FileType::ARC},
            {".tpr",   FileType::TPR},
            {".mol2",  FileType::MOL2},
            {".prm",   FileType::PRM},
            {".gro",   FileType::GRO},
            {".traj",  FileType::TRAJ}
    };

    auto it = mapping.find(extension);
    return it != mapping.end() ? it->second : FileType::UnKnown;

}

std::string choose_file(const std::string &prompt, bool exist, std::string ext, bool can_empty,
                        std::istream &in, std::ostream &out) {
    while (true) {
        std::string input_line = input(prompt, in, out);
        boost::trim(input_line);
        if (!input_line.empty()) {
            if (ext.length()) {
                if (ext_filename(input_line) != boost::to_lower_copy(ext)) {
                    out << "wrong file extesion name : must be " << ext << std::endl;
                    continue;
                }
            }
            if (!exist) return input_line;
            std::fstream in(input_line, std::ofstream::in);
            if (in.good()) {
                return input_line;
                break;
            } else {
                out << "The file is bad [retype]" << std::endl;
            }
        }
        if (can_empty) return "";
    }
}

double
atom_distance(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2, std::shared_ptr<Frame> &frame) {
    return std::sqrt(atom_distance2(atom1, atom2, frame));
}

double
atom_distance2(const std::shared_ptr<Atom> &atom1, const std::shared_ptr<Atom> &atom2, std::shared_ptr<Frame> &frame) {
    auto xr = atom1->x - atom2->x;
    auto yr = atom1->y - atom2->y;
    auto zr = atom1->z - atom2->z;
    frame->image(xr, yr, zr);
    return xr * xr + yr * yr + zr * zr;
}

po::options_description make_program_options() {
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "show this help message")
            ("topology,p", po::value<std::string>()->value_name("topology-file-name"), "topology file")
            ("file,f", po::value<std::vector<std::string>>()->multitoken()->composing()
                    ->value_name("trajectory-file-name"), "trajectory file")
            ("output,o", po::value<std::string>()->value_name("output-file-name"), "output file")
            ("prm", po::value<std::string>()->value_name("tinker-prm-file-name"), "force field file")
            ("target,x", po::value<std::string>()->value_name("trajectout-file-name"), "target trajectory file")
            ("script", po::value<std::string>()->value_name("script-content"), "script command for non-interactive use")
            ("script-file", po::value<std::string>()->value_name("script-file-name"), "read command from script file");

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

std::string getOutputFilename(const po::variables_map &vm) {
    return vm.count("output") ? vm["output"].as<std::string>() : choose_file("output file :", false);
}
std::string getTopologyFilename(const po::variables_map &vm) {
    return vm.count("topology") ? vm["topology"].as<std::string>() : choose_file("topology file :", true);
}

std::string getTrajectoryFilename(const po::variables_map &vm) {
    return vm.count("file") ? vm["file"].as<std::string>() : choose_file("trajectory file :", true);
}

std::string getPrmFilename(const po::variables_map &vm) {
    return vm.count("prm") ? vm["prm"].as<std::string>() : choose_file("Tinker prm file :", true);
}

std::size_t getDefaultVectorReserve() {
    auto p = std::getenv("ANALYSIS_VECTOR_RESERVE");
    return p ? std::stoi(p) : 100000;
}