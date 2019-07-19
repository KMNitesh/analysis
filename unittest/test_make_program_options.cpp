//
// Created by xiamr on 6/10/19.
//

#include <gmock/gmock.h>
#include <type_traits>
#include <string>
#include "../src/common.hpp"

using namespace testing;


TEST(make_program_options,TargetTrajectory) {
    const char *argv[] = {
            "analysis",
            "-x",
            "x.xtc"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    ASSERT_THAT(vm.count("target"), Eq(1));
}

TEST(make_program_options, InvalidOption) {
    const char *argv[] = {
            "analysis",
            "-a"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    ASSERT_ANY_THROW(po::store(po::command_line_parser(argc, argv).options(desc).run(), vm));

}

TEST(make_program_options, HelpOption) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.tpr",
            "-h"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    ASSERT_THAT(vm.count("help"), Eq(1));
}

TEST(make_program_options, TopologyOptionOnly) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.tpr"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.tpr"));
    ASSERT_THAT(vm.count("file"), Eq(0));
    ASSERT_THAT(vm.count("output"), Eq(0));
    ASSERT_THAT(vm.count("prm"), Eq(0));
    ASSERT_THAT(vm.count("help"), Eq(0));
}

TEST(make_program_options, TopologyWithOutpoutOption) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.tpr",
            "-o",
            "out.xvg"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.tpr"));
    ASSERT_THAT(vm.count("file"), Eq(0));
    ASSERT_THAT(vm.count("output"), Eq(1));
    ASSERT_THAT(vm["output"].as<std::string>(), Eq("out.xvg"));
    ASSERT_THAT(vm.count("prm"), Eq(0));
    ASSERT_THAT(vm.count("help"), Eq(0));
}

TEST(make_program_options, MultiTrajectoryWithCompose) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.mol2",
            "-f",
            "b.xtc",
            "-f",
            "c.nc"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.mol2"));

    ASSERT_THAT(vm.count("file"), Eq(1));

    auto trajs = vm["file"].as<std::vector<std::string>>();

    ASSERT_THAT(trajs.size(), Eq(2));
    ASSERT_THAT(trajs[0], Eq("b.xtc"));
    ASSERT_THAT(trajs[1], Eq("c.nc"));
}

TEST(make_program_options, MultiTrajectoryWithMultitoken) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.mol2",
            "-f",
            "b.xtc",
            "c.nc",
            "-o",
            "out.xvg"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.mol2"));

    ASSERT_THAT(vm.count("file"), Eq(1));

    auto trajs = vm["file"].as<std::vector<std::string>>();

    ASSERT_THAT(trajs.size(), Eq(2));
    ASSERT_THAT(trajs[0], Eq("b.xtc"));
    ASSERT_THAT(trajs[1], Eq("c.nc"));

    ASSERT_THAT(vm.count("output"), Eq(1));
    ASSERT_THAT(vm["output"].as<std::string>(), Eq("out.xvg"));
}

TEST(make_program_options, MultiTrajectoryWithMultitokenWithScript) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.mol2",
            "-f",
            "b.xtc",
            "c.nc",
            "-o",
            "out.xvg",
            "--script",
            "rotacf(vector = normalVector([:1], [@O], [@1500]), P=1, time_increment_ps = 0.1, max_time_grap_ps = 100)"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.mol2"));

    ASSERT_THAT(vm.count("file"), Eq(1));

    auto trajs = vm["file"].as<std::vector<std::string>>();

    ASSERT_THAT(trajs.size(), Eq(2));
    ASSERT_THAT(trajs[0], Eq("b.xtc"));
    ASSERT_THAT(trajs[1], Eq("c.nc"));

    ASSERT_THAT(vm.count("output"), Eq(1));
    ASSERT_THAT(vm["output"].as<std::string>(), Eq("out.xvg"));

    ASSERT_THAT(vm.count("script"), Eq(1));
    ASSERT_THAT(vm["script"].as<std::string>(),
                Eq("rotacf(vector = normalVector([:1], [@O], [@1500]), P=1, time_increment_ps = 0.1, max_time_grap_ps = 100)"));
}

TEST(make_program_options, ScriptFile) {
    const char *argv[] = {
            "analysis",
            "-p",
            "a.mol2",
            "-f",
            "b.xtc",
            "c.nc",
            "-o",
            "out.xvg",
            "--script-file",
            "a.script"
    };

    int argc = std::extent_v<decltype(argv)>;

    auto desc = make_program_options();
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    ASSERT_THAT(vm.count("topology"), Eq(1));
    ASSERT_THAT(vm["topology"].as<std::string>(), Eq("a.mol2"));

    ASSERT_THAT(vm.count("file"), Eq(1));

    auto trajs = vm["file"].as<std::vector<std::string>>();

    ASSERT_THAT(trajs.size(), Eq(2));
    ASSERT_THAT(trajs[0], Eq("b.xtc"));
    ASSERT_THAT(trajs[1], Eq("c.nc"));

    ASSERT_THAT(vm.count("output"), Eq(1));
    ASSERT_THAT(vm["output"].as<std::string>(), Eq("out.xvg"));

    ASSERT_THAT(vm.count("script-file"), Eq(1));
    ASSERT_THAT(vm["script-file"].as<std::string>(), Eq("a.script"));
}