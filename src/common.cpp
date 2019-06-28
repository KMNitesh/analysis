//
// Created by xiamr on 3/17/19.
//

#include "config.h"
#include <tuple>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <regex>
#include "common.hpp"
#include "frame.hpp"
#include "forcefield.hpp"
#include <boost/filesystem.hpp>


bool enable_read_velocity = false;
bool enable_tbb = false;
bool enable_outfile = false;

Forcefield forcefield;
bool enable_forcefield = false;

std::fstream outfile;


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


std::string input(const std::string &prompt) {
    std::cout << prompt;
    std::string inputline;
    std::getline(std::cin, inputline);
    if (isatty(STDIN_FILENO) == 0) {
        std::cout << inputline << std::endl;
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
//    if (!boost::filesystem::is_regular_file(filename)){
//        throw std::runtime_error("ERROR !! not regular file !");
//    }
    auto extension = boost::filesystem::extension(filename);
    boost::to_lower(extension);
    /*
     *      XTC,TRR,NC,ARC,TPR,MOL2,PRM, UnKnown
     */
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
            {".gro",   FileType::GRO}
    };

    auto it = mapping.find(extension);
    return it != mapping.end() ? it->second : FileType::UnKnown;

}

std::string choose_file(const std::string &prompt, bool exist, std::string ext, bool can_empty) {
    while (true) {
        std::string input_line = input(prompt);
        boost::trim(input_line);
        if (!input_line.empty()) {
            if (ext.length()) {
                if (ext_filename(input_line) != boost::to_lower_copy(ext)) {
                    std::cerr << "wrong file extesion name : must be " << ext << std::endl;
                    continue;
                }
            }
            if (!exist) return input_line;
            std::fstream in(input_line, std::ofstream::in);
            if (in.good()) {
                return input_line;
                break;
            } else {
                std::cerr << "The file is bad [retype]" << std::endl;
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
            ("topology,p", po::value<std::string>(), "topology file")
            ("file,f", po::value<std::vector<std::string>>()->multitoken()->composing(), "trajectory file")
            ("output,o", po::value<std::string>(), "output file")
            ("prm", po::value<std::string>(), "force field file")
            ("target,x", po::value<std::string>(), "target trajectory file");

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
