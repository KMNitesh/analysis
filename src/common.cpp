//
// Created by xiamr on 3/17/19.
//

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "common.hpp"


bool file_exist(const std::string &name) {
    std::fstream in(name, std::ofstream::in);
    return in.good();
}


std::vector<std::string> split(const std::string &str, const std::string &sep) {
    std::vector<std::string> ret_;
    boost::split(ret_, str, boost::is_any_of(sep));
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

    return boost::to_lower_copy(field[field.size() - 1]);
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



