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
#include <chrono>
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
#include "TypeUtility.hpp"
#include "Trajconv.hpp"
#include "RotAcf.hpp"
#include "RotAcfCutoff.hpp"
#include "Interpreter.hpp"
#include "ThrowAssert.hpp"
#include "DipoleVectorSelector.hpp"
#include "NormalVectorSelector.hpp"
#include "TwoAtomVectorSelector.hpp"
#include "RadicalDistribtuionFunction.hpp"
#include "DemixIndexOfTwoGroup.hpp"
#include "ResidenceTime.hpp"
#include "Diffuse.hpp"
#include "DiffuseCutoff.hpp"
#include "AngleDistributionBetweenTwoVectorWithCutoff.hpp"

using namespace std;

int executeAnalysis(const vector<string> &xyzfiles, int argc, char *const *argv, const string &scriptContent,
                    boost::optional<string> &script_file, boost::optional<string> &topology,
                    boost::optional<string> &forcefield_file, const boost::optional<string> &output_file,
                    shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list,
                    int start,
                    int total_frames,
                    int step_size,
                    int nthreads
);

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
                   const boost::program_options::variables_map &vm, std::vector<std::string> &xyzfiles,
                   int argc, char *argv[]) {
    string scriptContent;

    boost::optional<string> script_file;
    if (vm.count("script")) {
        scriptContent = vm["script"].as<string>();

    } else if (vm.count("script-file")) {

        try {
            script_file = vm["script-file"].as<string>();
            ifstream f(script_file.value());
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

    boost::optional<string> topology;
    if (vm.count("topology")) topology = vm["topology"].as<string>();
    boost::optional<string> forcefield_file;
    if (vm.count("prm")) forcefield_file = vm["prm"].as<string>();
    boost::optional<string> output_file;
    if (vm.count("output")) output_file = vm["output"].as<string>();

    boost::any ast;
    auto it = scriptContent.begin();
    bool status = qi::phrase_parse(it, scriptContent.end(), InterpreterGrammarT(), SkipperT(), ast);

    if (!(status and it == scriptContent.end())) {
        std::cerr << "Syntax Parse Error\n";
        std::cerr << "error-pos : " << std::endl;
        std::cout << scriptContent << std::endl;
        for (auto iter = scriptContent.begin(); iter != it; ++iter) std::cerr << " ";
        std::cerr << "^" << std::endl;
        exit(EXIT_FAILURE);
    }

    Interpreter interpreter;

    auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();

    interpreter.registerFunction(
                    "rdf", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<RadicalDistribtuionFunction>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))),
                                    AutoConvert(get<3>(args.at(5))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function rdf (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("M")
            .addArgument<Atom::Node>("L")
            .addArgument<double, int>("max_dist", 10.0)
            .addArgument<double, int>("width", 0.01)
            .addArgument<bool>("intramol", false)
            .addArgument<string>("out");


    interpreter.registerFunction(
                    "rotacf", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<RotAcf>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function rotacf (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        cout << task->description();
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<shared_ptr<VectorSelector>>("vector")
            .addArgument<int>("P")
            .addArgument<double, int>("time_increment_ps", 0.1)
            .addArgument<double, int>("max_time_grap_ps")
            .addArgument<string>("out");


    interpreter.registerFunction(
                    "rotacfCutoff", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<RotAcfCutoff>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))),
                                    AutoConvert(get<3>(args.at(5))),
                                    AutoConvert(get<3>(args.at(6))),
                                    AutoConvert(get<3>(args.at(7))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function rotacfCutoff (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        cout << task->description();
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("M")
            .addArgument<Atom::Node>("L")
            .addArgument<shared_ptr<VectorSelector>>("vector")
            .addArgument<int>("P")
            .addArgument<double, int>("cutoff")
            .addArgument<double, int>("time_increment_ps", 0.1)
            .addArgument<double, int>("max_time_grap_ps")
            .addArgument<string>("out");


    interpreter.registerFunction(
                    "demix", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<DemixIndexOfTwoGroup>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function demix (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("component1")
            .addArgument<Atom::Node>("component2")
            .addArgument<Grid>("grid")
            .addArgument<string>("out");

    interpreter.registerFunction(
                    "residenceTime", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<ResidenceTime>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function residenceTime (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("M")
            .addArgument<Atom::Node>("L")
            .addArgument<double, int>("cutoff")
            .addArgument<double, int>("time_star")
            .addArgument<string>("out");


    interpreter.registerFunction(
                    "diffuse", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<Diffuse>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function diffuse (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        cout << task->description();
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("mask")
            .addArgument<double, int>("time_increment_ps", 0.1)
            .addArgument<int>("total_frames")
            .addArgument<string>("out");

    interpreter.registerFunction(
                    "diffuseCutoff", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<DiffuseCutoff>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))));
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function diffuseCutoff (" << __FILE__ << ":" << __LINE__ << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        cout << task->description();
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("M")
            .addArgument<Atom::Node>("L")
            .addArgument<double, int>("cutoff")
            .addArgument<double, int>("time_increment_ps", 0.1)
            .addArgument<string>("out");


    interpreter.registerFunction(
                    "crossAngleCutoff", [&task_list](auto &args) -> boost::any {
                        auto task = make_shared<AngleDistributionBetweenTwoVectorWithCutoff>();
                        try {
                            task->setParameters(
                                    AutoConvert(get<3>(args.at(0))),
                                    AutoConvert(get<3>(args.at(1))),
                                    AutoConvert(get<3>(args.at(2))),
                                    AutoConvert(get<3>(args.at(3))),
                                    AutoConvert(get<3>(args.at(4))),
                                    AutoConvert(get<3>(args.at(5))),
                                    AutoConvert(get<3>(args.at(6))),
                                    AutoConvert(get<3>(args.at(7))),
                                    AutoConvert(get<3>(args.at(8)))
                            );
                        } catch (std::exception &e) {
                            cerr << e.what() << " for function crossAngleCutoff (" << __FILE__ << ":" << __LINE__
                                 << ")\n";
                            exit(EXIT_FAILURE);
                        }
                        task_list->emplace_back(task);
                        cout << task->description();
                        return shared_ptr<BasicAnalysis>(task);
                    })
            .addArgument<Atom::Node>("M")
            .addArgument<Atom::Node>("L")
            .addArgument<shared_ptr<VectorSelector>>("vector1")
            .addArgument<shared_ptr<VectorSelector>>("vector2")
            .addArgument<double, int>("angle_max", 180)
            .addArgument<double, int>("angle_width", 0.5)
            .addArgument<double, int>("cutoff1", 0)
            .addArgument<double, int>("cutoff2")
            .addArgument<string>("out");
    //  const Atom::Node &M,
    //        const Atom::Node &L,
    //        std::shared_ptr<VectorSelector> vector1,
    //        std::shared_ptr<VectorSelector> vector2,
    //        double angle_max,
    //        double angle_width,
    //        double cutoff1,
    //        double cutoff2,
    //        const std::string &outfilename

    interpreter.registerFunction(
            "DipoleVector", [](auto &args) -> boost::any {
                auto vector = make_shared<DipoleVectorSelector>();
                try {
                    vector->setParameters(boost::any_cast<Atom::Node>(get<3>(args.at(0))));
                } catch (std::exception &e) {
                    cerr << e.what() << " for function DipoleVector (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
                return shared_ptr<VectorSelector>(vector);
            }).addArgument<Atom::Node>("mask");

    interpreter.registerFunction(
            "NormalVector", [](auto &args) -> boost::any {
                auto vector = make_shared<NormalVectorSelector>();
                try {
                    vector->setParameters(AutoConvert(get<3>(args.at(0))),
                                          AutoConvert(get<3>(args.at(1))),
                                          AutoConvert(get<3>(args.at(2))));
                } catch (std::exception &e) {
                    cerr << e.what() << " for function NormalVector (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
                return shared_ptr<VectorSelector>(vector);
            }).addArgument<Atom::Node>("mask1").addArgument<Atom::Node>("mask2").addArgument<Atom::Node>("mask3");

    interpreter.registerFunction(
            "TwoAtomVector", [](auto &args) -> boost::any {
                auto vector = make_shared<TwoAtomVectorSelector>();
                try {
                    vector->setParameters(AutoConvert(get<3>(args.at(0))), AutoConvert(get<3>(args.at(1))));
                } catch (std::exception &e) {
                    cerr << e.what() << " for function TwoAtomVector (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
                return shared_ptr<VectorSelector>(vector);
            }).addArgument<Atom::Node>("mask1").addArgument<Atom::Node>("mask2");

    interpreter.registerFunction(
            "Grid", [](auto &args) -> boost::any {
                try {
                    return Grid{AutoConvert(get<3>(args.at(0))),
                                AutoConvert(get<3>(args.at(1))),
                                AutoConvert(get<3>(args.at(2)))};
                } catch (std::exception &e) {
                    cerr << e.what() << " for function Grid (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
            }).addArgument<int>("x").addArgument<int>("y").addArgument<int>("z");

    interpreter.registerFunction(
            "readTop", [&topology](auto &args) -> boost::any {
                try {
                    topology = AutoConvert(get<3>(args.at(0)));
                    return topology.value();
                } catch (std::exception &e) {
                    cerr << e.what() << " for function readTop (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
            }).addArgument<string>("file");

    interpreter.registerFunction(
            "trajin", [&xyzfiles](auto &args) -> boost::any {
                try {
                    string xyz = AutoConvert(get<3>(args.at(0)));
                    xyzfiles.emplace_back(xyz);
                    return xyz;
                } catch (std::exception &e) {
                    cerr << e.what() << " for function trajin (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
            }).addArgument<string>("file");

    interpreter.registerFunction(
            "readFF", [&forcefield_file](auto &args) -> boost::any {
                try {
                    forcefield_file = AutoConvert(get<3>(args.at(0)));
                    return forcefield_file.value();
                } catch (std::exception &e) {
                    cerr << e.what() << " for function readFF (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }
            }).addArgument<string>("file");

    interpreter.registerFunction("go", [&](auto &args) -> boost::any {
                int start, total_frames, step_size, nthreads;
                try {
                    start = AutoConvert(get<3>(args.at(0)));
                    total_frames = AutoConvert(get<3>(args.at(1)));
                    step_size = AutoConvert(get<3>(args.at(2)));
                    nthreads = AutoConvert(get<3>(args.at(3)));
                } catch (std::exception &e) {
                    cerr << e.what() << " for function go (" << __FILE__ << ":" << __LINE__ << ")\n";
                    exit(EXIT_FAILURE);
                }

                return executeAnalysis(xyzfiles, argc, argv, scriptContent, script_file, topology, forcefield_file,
                                       output_file, task_list,
                                       start, total_frames, step_size, nthreads);

            }).addArgument<int>("start", 1).addArgument<int>("end", 0)
            .addArgument<int>("step", 1).addArgument<int>("nthreads", 0);

    interpreter.execute(ast);

    if (!task_list->empty()) {
        cout << "Program terminated with " << task_list->size() << " task(s) unperformed !\n";
        cout << "Don't forget to put go function in the end of script\n";
    }
    cout << " < Mission Complete >\n";
}

int executeAnalysis(const vector<string> &xyzfiles, int argc, char *const *argv, const string &scriptContent,
                    boost::optional<string> &script_file, boost::optional<string> &topology,
                    boost::optional<string> &forcefield_file, const boost::optional<string> &output_file,
                    shared_ptr<list<shared_ptr<BasicAnalysis>>> &task_list, int start, int total_frames, int step_size,
                    int nthreads) {
    if (task_list->empty()) {
        cerr << "Empty task in the pending list, skip go function ...\n";
        return 0;
    }

    auto start_time = chrono::steady_clock::now();

    cout << boost::format("Start Process...  start = %d, end = %s, step = %d, nthreads = %s\n")
            % start % (total_frames == 0 ? "all" : to_string(total_frames)) % step_size
            % (nthreads == 0 ? "automatic" : to_string(nthreads));

    if (start <= 0) {
        cerr << "start frame cannot less than 1\n";
        exit(EXIT_FAILURE);
    }
    if (total_frames <= start and total_frames != 0) {
        cerr << format("end(%d) frame cannot less than start(%d) frame\n", total_frames, start);
        exit(EXIT_FAILURE);
    }
    if (step_size <= 0) {
        cerr << "frame step cannot less than 1\n";
        exit(EXIT_FAILURE);
    }
    if (nthreads < 0) {
        cerr << "thread number cannot less than zero\n";
        exit(EXIT_FAILURE);
    }
    if (enable_forcefield) {
        if (topology && getFileType(topology.value()) != FileType::ARC) {

        } else if (forcefield_file) {
            if (boost::filesystem::exists(forcefield_file.value())) {
                forcefield.read(forcefield_file.value());
            } else {
                cerr << "force field file " << forcefield_file.value() << " is bad  !\n";
                exit(EXIT_FAILURE);
            }
        } else {
            cerr << "force field file not given !\n";
            exit(EXIT_FAILURE);
        }
    }

    tbb::task_scheduler_init tbb_init(tbb::task_scheduler_init::deferred);
    if (nthreads > sysconf(_SC_NPROCESSORS_ONLN)) {
        cout << "WARNING !! nthreads larger than max core of system CPU, will use automatic mode";
        nthreads = 0;
    }


    tbb_init.initialize(nthreads == 0 ? tbb::task_scheduler_init::automatic : nthreads);

    int current_frame_num = 0;

    auto reader = make_shared<TrajectoryReader>();
    bool b_added_topology = true;
    if (boost::algorithm::one_of_equal<initializer_list<FileType>>(
            {FileType::NC, FileType::XTC, FileType::TRR}, getFileType(xyzfiles[0]))) {
        b_added_topology = false;
    } else {
        if (topology) {
            cerr << "WRANING !!  do not use topolgy file !\bn";
        }
    }

    if (xyzfiles.empty()) {
        cerr << "trajectory file not set\n";
        exit(EXIT_FAILURE);
    }
    for (auto &xyzfile : xyzfiles) {
        reader->add_filename(xyzfile);
        string ext = ext_filename(xyzfile);
        if (ext == "traj" && xyzfiles.size() != 1) {
            cout << "traj file can not use multiple files" << endl;
            exit(EXIT_FAILURE);
        }
        if (!b_added_topology) {
            if (topology) {
                if (boost::filesystem::exists(topology.value())) {
                    reader->add_topology(topology.value());
                    b_added_topology = true;
                    continue;
                }
                cerr << "topology file " << topology.value() << " is bad ! please retype !\n";
                exit(EXIT_FAILURE);
            }
            cerr << "topology file  not given !\n";
            exit(EXIT_FAILURE);
        }
    }

    if (enable_outfile) {
        if (output_file) {
            cerr << "Output file option do not need\n";
            exit(EXIT_FAILURE);
        }
    }
    shared_ptr<Frame> frame;
    int Clear = 0;
    while ((frame = reader->readOneFrame())) {
        current_frame_num++;
        if (total_frames != 0 and current_frame_num > total_frames)
            break;
        if (current_frame_num % 10 == 0) {
            if (Clear) {
                cout << "\r";
            }
            cout << "Processing Coordinate Frame  " << current_frame_num << "   " << flush;
            Clear = 1;
        }
        if (current_frame_num >= start && (current_frame_num - start) % step_size == 0) {
            if (current_frame_num == start) {
                if (forcefield.isValid()) {
                    forcefield.assign_forcefield(frame);
                }
                processFirstFrame(frame, task_list);
            }
            processOneFrame(frame, task_list);
        }
    }
    cout << '\n';

    tbb::parallel_for_each(*task_list, [&](auto &task) {
        ofstream ofs(task->getOutfileName());
        ofs << "#  workdir > " << boost::filesystem::current_path() << '\n';
        ofs << "#  cmdline > " << print_cmdline(argc, argv) << '\n';
        if (script_file) {
            ofs << "#  script-file > " << script_file.value() << '\n';
        }
        ofs << "#  script BEGIN>\n" << scriptContent << "\n#  script END>\n";
        task->print(ofs);
    });

    int totol_task_count = task_list->size();

    cout << "Complete " << totol_task_count << " task(s)         Run Time "
         << chrono_cast(chrono::steady_clock::now() - start_time) << '\n';

    task_list->clear();
    return totol_task_count;
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

    for (;;) {
        if (vm.count("topology") && getFileType(vm["topology"].as<std::string>()) != FileType::ARC) {
            break;
        }
        if (vm.count("prm")) {
            auto ff = vm["prm"].as<std::string>();
            if (boost::filesystem::exists(ff)) {
                forcefield.read(ff);
                break;
            } else if (enable_forcefield) {
                std::cout << "force field file " << ff << " is bad ! please retype !" << std::endl;
            } else {
                break;
            }
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
    if (boost::algorithm::one_of_equal<std::initializer_list<FileType>>(
            {FileType::NC, FileType::XTC, FileType::TRR}, getFileType(xyzfiles[0]))) {
        b_added_topology = false;
    } else {
        if (vm.count("topology")) {
            std::cerr << "WRANING !!  do not use topolgy file !\bn";
        }
    }
    for (auto &xyzfile : xyzfiles) {
        reader->add_filename(xyzfile);
        if (getFileType(xyzfile) == FileType::TRAJ && xyzfiles.size() != 1) {
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

    std::ofstream outfile;
    if (enable_outfile) {
        outfile.open(getOutputFilename(vm));
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
                if (forcefield.isValid()) {
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
        outfile << "#  start frame  > " << start << '\n';
        outfile << "#  step (frame) > " << step_size << '\n';
        outfile << "#  end   frame  > " << (total_frames == 0 ? "all" : std::to_string(total_frames)) << '\n';
    }

    for (auto &task : *task_list) {
        task->print(outfile);
    }
    if (outfile.is_open()) outfile.close();
    std::cout << "Mission Complete" << std::endl;
}

