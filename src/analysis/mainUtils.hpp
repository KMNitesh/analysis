//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_MAINUTILS_HPP
#define TINKER_MAINUTILS_HPP

#include "std.hpp"
#include <boost/optional.hpp>

class BasicAnalysis;

class TrajectoryReader;

namespace boost::program_options {
    class variables_map;
}


void fastTrajectoryConvert(const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles);

void printTopolgy(const boost::program_options::variables_map &vm);

void processTrajectory(const boost::program_options::options_description &desc,
                       const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles,
                       int argc, char *argv[]);

void executeScript(const boost::program_options::options_description &desc,
                   const boost::program_options::variables_map &vm, std::vector<std::string> &xyzfiles,
                   int argc, char *argv[]);

std::shared_ptr<Frame>
getFrame(std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> &task_list,
         const int start, const int step_size, const int total_frames,
         std::shared_ptr<TrajectoryReader> &reader);

int
executeAnalysis(const std::vector<std::string> &xyzfiles, int argc, char *const *argv, const std::string &scriptContent,
                boost::optional<std::string> &script_file, boost::optional<std::string> &topology,
                boost::optional<std::string> &forcefield_file, const boost::optional<std::string> &output_file,
                std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> &task_list,
                int start,
                int total_frames,
                int step_size,
                int nthreads
);


#endif //TINKER_MAINUTILS_HPP
