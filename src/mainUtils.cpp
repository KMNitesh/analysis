//
// Created by xiamr on 6/14/19.
//


#include <string>
#include <tbb/tbb.h>
#include <cmath>
#include <memory>
#include <list>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>

#include "PrintTopolgy.hpp"
#include "taskMenu.hpp"
#include "BasicAnalysis.hpp"
#include "forcefield.hpp"

#include "trajectoryreader.hpp"

#include "mainUtils.hpp"
#include "common.hpp"
#include "frame.hpp"

#include "Trajconv.hpp"
#include "LanguageGrammar.hpp"
#include "RotAcf.hpp"
#include "RotAcfCutoff.hpp"


using namespace std;

void processOneFrame(shared_ptr<Frame> &frame,
                     shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list) {
    for (auto &task : *task_list) {
        task->process(frame);
    }
}

void processFirstFrame(shared_ptr<Frame> &frame,
                       shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list) {
    for (auto &task : *task_list) {
        task->processFirstFrame(frame);
    }
}

void fastTrajectoryConvert(const boost::program_options::variables_map &vm, const vector<std::string> &xyzfiles) {
    if (!vm.count("file")) {
        cerr << "ERROR !! trajctory file not set !\n";
        exit(EXIT_FAILURE);
    }

    auto reader = make_shared<TrajectoryReader>();
    if (boost::algorithm::one_of_equal<std::initializer_list<FileType>>(
            {FileType::NC, FileType::XTC, FileType::TRR}, getFileType(xyzfiles[0]))) {

        if (!vm.count("topology")) {
            cerr << "ERROR !! topology file not set !\n";
            exit(EXIT_FAILURE);
        }
        string topol = vm["topology"].as<string>();
        if (!boost::filesystem::exists(topol)) {
            cerr << "ERROR !! topology file not exist !\n";
            exit(EXIT_FAILURE);
        }
        reader->add_topology(topol);
    } else {
        if (vm.count("topology")) {
            cerr << "WRANING !!  do not use topolgy file !\bn";
        }
    }

    for (auto &traj : xyzfiles) {
        string ext = ext_filename(traj);
        if (ext == "traj" && xyzfiles.size() != 1) {
            cerr << "traj file can not use multiple files\n";
            exit(EXIT_FAILURE);
        }
        reader->add_filename(traj);
    }

    shared_ptr<Frame> frame;
    int current_frame_num = 0;

    Trajconv writer;
    try {
        writer.fastConvertTo(vm["target"].as<string>());
    } catch (std::runtime_error &e) {
        std::cerr << e.what() << '\n';
        exit(EXIT_FAILURE);
    }

    cout << "Fast Trajectory Convert...\n";
    while ((frame = reader->readOneFrame())) {
        current_frame_num++;
        if (current_frame_num == 1) {
            cout << "Processing Coordinate Frame  " << current_frame_num << "   " << flush;
            writer.processFirstFrame(frame);
        } else if (current_frame_num % 10 == 0) {
            cout << "\rProcessing Coordinate Frame  " << current_frame_num << "   " << flush;
        }
        writer.process(frame);
    }
    writer.CleanUp();
    cout << "\nMission Complete" << endl;
}

void printTopolgy(const boost::program_options::variables_map &vm) {
    PrintTopolgy printer;
    if (vm.count("prm")) {
        auto ff = vm["prm"].as<std::string>();
        if (boost::filesystem::exists(ff)) {
            forcefield.read(ff);
        } else {
            std::cout << "force field file " << ff << " is bad !" << std::endl;
        }
    }
    if (vm.count("topology")) {
        string topol = vm["topology"].as<string>();
        if (boost::filesystem::exists(topol)) {
            printer.action(topol);
            exit(EXIT_SUCCESS);
        }
        cout << "topology file " << topol << " is bad ! please retype !" << endl;
    }
    printer.action(choose_file("input topology file : ", true));
}

void executeScript(const boost::program_options::options_description &desc,
                   const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles,
                   int argc, char *argv[]) {
    if (!vm.count("file")) {
        std::cerr << "input trajectory file is not set !" << std::endl;
        std::cerr << desc;
        exit(EXIT_FAILURE);
    }

    string scriptContent;

    if (vm.count("script")) {
        scriptContent = vm["script"].as<string>();

    } else if (vm.count("script-file")) {

        try {
            ifstream f(vm["script-file"].as<string>());
            string content((std::istreambuf_iterator<char>(f)),
                           std::istreambuf_iterator<char>());
            scriptContent = content;

        } catch (std::exception &e) {
            std::cerr << e.what() << '\n';
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "program option error  \n";
        exit(EXIT_FAILURE);
    }


    LanguageGrammar<std::string::iterator, qi::ascii::space_type> grammar;


    Language ast;
    auto it = scriptContent.begin();
    bool status = qi::phrase_parse(it, scriptContent.end(), grammar, qi::ascii::space, ast);

    if (!(status and it == scriptContent.end())) {
        std::cerr << "Syntax Parse Error\n";
        std::cerr << "error-pos : " << std::endl;
        std::cout << scriptContent << std::endl;
        for (auto iter = scriptContent.begin(); iter != it; ++iter) std::cerr << " ";
        std::cerr << "^" << std::endl;
        exit(EXIT_FAILURE);
    }


    struct ExecutionEngine : boost::static_visitor<shared_ptr<BasicAnalysis>> {

        shared_ptr<BasicAnalysis> operator()(RotAcfNode &ast) {
            auto task = make_shared<RotAcf>();
            try {
                task->readAST(ast);
            } catch (std::exception &e) {
                std::cerr << e.what() << '\n';
                exit(EXIT_FAILURE);
            }
            return task;
        }

        shared_ptr<BasicAnalysis> operator()(RotAcfCutoffNode &ast) {
            auto task = make_shared<RotAcfCutoff>();
            try {
                task->readAST(ast);
            } catch (std::exception &e) {
                std::cerr << e.what() << '\n';
                exit(EXIT_FAILURE);
            }
            return task;
        }

        shared_ptr<BasicAnalysis> operator()(GoNode &ast) {
            if (ast->start <= 0) {
                cerr << "start frame cannot less than 1\n";
                exit(EXIT_FAILURE);
            } else {
                this->start = ast->start;
            }

            if (ast->end <= this->start and ast->end != 0) {
                cerr << "end frame cannot less than start frame\n";
                exit(EXIT_FAILURE);
            } else {
                end = ast->end;
            }

            if (ast->step <= 0) {
                cerr << "frame step cannot less than 1\n";
                exit(EXIT_FAILURE);
            } else {
                this->step = ast->step;
            }

            this->nthreads = ast->nthreads;

            return {};
        }

        int start = 1;
        int end = 0;
        int step = 1;
        int nthreads = 0;

    } engine;


    auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();


    for (auto &stmt : ast) {
        auto task = boost::apply_visitor(engine, *stmt);
        if (task) {
            task_list->push_back(task);
        }
    }

    if (enable_forcefield) {
        if (vm.count("topology") && getFileType(vm["topology"].as<std::string>()) != FileType::ARC) {

        } else if (vm.count("prm")) {
            auto ff = vm["prm"].as<std::string>();
            if (boost::filesystem::exists(ff)) {
                forcefield.read(ff);
            } else {
                std::cerr << "force field file " << ff << " is bad  !\n";
                exit(EXIT_FAILURE);
            }
        } else {
            std::cerr << "force field file not given !\n";
            exit(EXIT_FAILURE);
        }
    }

    int start = engine.start;
    int step_size = engine.step;
    int total_frames = engine.end;

    tbb::task_scheduler_init tbb_init(tbb::task_scheduler_init::deferred);
    if (engine.nthreads > sysconf(_SC_NPROCESSORS_ONLN)) {
        cout << "WARNING !! nthreads larger than max core of system CPU, will use automatic mode";
        engine.nthreads = 0;
    }

    if (enable_tbb) {
        tbb_init.initialize(engine.nthreads == 0 ? tbb::task_scheduler_init::automatic : engine.nthreads);
    }
    int current_frame_num = 0;

    auto reader = std::make_shared<TrajectoryReader>();
    bool b_added_topology = true;
    auto ext = ext_filename(xyzfiles[0]);
    if (ext == "nc" or ext == "xtc" or ext == "trr") {
        b_added_topology = false;
    } else {
        if (vm.count("topology")) {
            std::cerr << "WRANING !!  do not use topolgy file !\bn";
        }
    }
    for (auto &xyzfile : xyzfiles) {
        reader->add_filename(xyzfile);
        std::string ext = ext_filename(xyzfile);
        if (ext == "traj" && xyzfiles.size() != 1) {
            std::cout << "traj file can not use multiple files" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!b_added_topology) {
            if (vm.count("topology")) {
                std::string topol = vm["topology"].as<std::string>();
                if (boost::filesystem::exists(topol)) {
                    reader->add_topology(topol);
                    b_added_topology = true;
                    continue;
                }
                std::cerr << "topology file " << topol << " is bad ! please retype !\n";
                exit(EXIT_FAILURE);
            }
            std::cerr << "topology file  not given !\n";
            exit(EXIT_FAILURE);
        }
    }

    std::fstream outfile;
    if (enable_outfile) {
        if (vm.count("output")) {
            cerr << "Output file option do not need\n";
            exit(EXIT_FAILURE);
        }
    }
    std::shared_ptr<Frame> frame;
    int Clear = 0;
    while ((frame = reader->readOneFrame())) {
        current_frame_num++;
        if (total_frames != 0 and current_frame_num > total_frames)
            break;
        if (current_frame_num % 10 == 0) {
            if (Clear) {
                std::cout << "\r";
            }
            std::cout << "Processing Coordinate Frame  " << current_frame_num << "   " << std::flush;
            Clear = 1;
        }
        if (current_frame_num >= start && (current_frame_num - start) % step_size == 0) {
            if (current_frame_num == start) {
                if (forcefield.isVaild()) {
                    forcefield.assign_forcefield(frame);
                }
                processFirstFrame(frame, task_list);
            }
            processOneFrame(frame, task_list);
        }
    }
    std::cout << std::endl;

    for (auto &task : *task_list) {
        ofstream ofs(task->getOutfileName());
        ofs << "#  workdir > " << boost::filesystem::current_path() << '\n';
        ofs << "#  cmdline > " << print_cmdline(argc, argv) << '\n';
        if (vm.count("script-file")) {
            ofs << "#  script-file > " << vm["script-file"].as<string>() << '\n';
        }
        ofs << "#  script BEGIN>\n" << scriptContent << "\n#  script END>\n";
        task->print(ofs);
    }
    std::cout << "Mission Complete\n";
}

void processTrajectory(const boost::program_options::options_description &desc,
                       const boost::program_options::variables_map &vm, const std::vector<std::string> &xyzfiles,
                       int argc, char *argv[]) {
    if (!vm.count("file")) {
        std::cerr << "input trajectory file is not set !" << std::endl;
        std::cerr << desc;
        exit(EXIT_FAILURE);
    }

    auto task_list = getTasks();

    while (enable_forcefield) {
        if (vm.count("topology") && getFileType(vm["topology"].as<std::string>()) != FileType::ARC) {
            break;
        }
        if (vm.count("prm")) {
            auto ff = vm["prm"].as<std::string>();
            if (boost::filesystem::exists(ff)) {
                forcefield.read(ff);
                break;
            }
            std::cout << "force field file " << ff << " is bad ! please retype !" << std::endl;

        }
        forcefield.read(choose_file("force field filename:", true));
        break;
    }
    int start = choose(1, INT32_MAX, "Enter the start frame[1]:", true, 1);
    int step_size = choose(1, INT32_MAX, "Enter the step size[1]:", true, 1);
    int total_frames = choose(0, INT32_MAX, "How many frames to read [all]:", true);

    tbb::task_scheduler_init tbb_init(tbb::task_scheduler_init::deferred);
    if (enable_tbb) {
        int threads = choose<int>(0, sysconf(_SC_NPROCESSORS_ONLN),
                                  "How many cores to used in parallel[automatic]:", true);
        tbb_init.initialize(threads == 0 ? tbb::task_scheduler_init::automatic : threads);
    }
    int current_frame_num = 0;

    auto reader = std::make_shared<TrajectoryReader>();
    bool b_added_topology = true;
    auto ext = ext_filename(xyzfiles[0]);
    if (ext == "nc" or ext == "xtc" or ext == "trr") {
        b_added_topology = false;
    } else {
        if (vm.count("topology")) {
            std::cerr << "WRANING !!  do not use topolgy file !\bn";
        }
    }
    for (auto &xyzfile : xyzfiles) {
        reader->add_filename(xyzfile);
        std::string ext = ext_filename(xyzfile);
        if (ext == "traj" && xyzfiles.size() != 1) {
            std::cout << "traj file can not use multiple files" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!b_added_topology) {
            if (vm.count("topology")) {
                std::string topol = vm["topology"].as<std::string>();
                if (boost::filesystem::exists(topol)) {
                    reader->add_topology(topol);
                    b_added_topology = true;
                    continue;
                }
                std::cout << "topology file " << topol << " is bad ! please retype !" << std::endl;
            }
            reader->add_topology(choose_file("input topology file : ", true));
            b_added_topology = true;
        }
    }


    auto input_line = input("Do you want to use multiple files [No]:");
    boost::trim(input_line);
    if (!input_line.empty()) {
        if (input_line[0] == 'Y' or input_line[0] == 'y') {
            while (true) {
                input_line = choose_file("next file [Enter for End]:", true, "", true);
                if (input_line.empty()) break;
                if (ext_filename(input_line) == "traj") {
                    std::cout << "traj file can not use multiple files [retype]" << std::endl;
                    continue;
                }
                reader->add_filename(input_line);
            }
        }
    }

    std::fstream outfile;
    if (enable_outfile) {
        outfile.open(vm.count("output") ? vm["output"].as<std::string>() : choose_file("Output file: ", false),
                     std::ios_base::out);
    }
    std::shared_ptr<Frame> frame;
    int Clear = 0;
    while ((frame = reader->readOneFrame())) {
        current_frame_num++;
        if (total_frames != 0 and current_frame_num > total_frames)
            break;
        if (current_frame_num % 10 == 0) {
            if (Clear) {
                std::cout << "\r";
            }
            std::cout << "Processing Coordinate Frame  " << current_frame_num << "   " << std::flush;
            Clear = 1;
        }
        if (current_frame_num >= start && (current_frame_num - start) % step_size == 0) {
            if (current_frame_num == start) {
                if (forcefield.isVaild()) {
                    forcefield.assign_forcefield(frame);
                }
                processFirstFrame(frame, task_list);
            }
            processOneFrame(frame, task_list);
        }
    }
    std::cout << std::endl;


    if (outfile.is_open()) {
        outfile << "#  workdir > " << boost::filesystem::current_path() << '\n';
        outfile << "#  cmdline > " << print_cmdline(argc, argv) << '\n';
    }

    for (auto &task : *task_list) {
        task->print(outfile);
    }
    if (outfile.is_open()) outfile.close();
    std::cout << "Mission Complete" << std::endl;
}

