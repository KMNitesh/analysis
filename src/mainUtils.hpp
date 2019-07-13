//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_MAINUTILS_HPP
#define TINKER_MAINUTILS_HPP

#include <vector>
#include <string>


namespace boost::program_options {
    class variables_map;
}


void fastTrajectoryConvert(const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles);

void printTopolgy(const boost::program_options::variables_map &vm);

void processTrajectory(const boost::program_options::options_description &desc,
                       const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles,
                       int argc, char *argv[]);

void executeScript(const boost::program_options::options_description &desc,
                   const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles,
                   int argc, char *argv[]);


#endif //TINKER_MAINUTILS_HPP
