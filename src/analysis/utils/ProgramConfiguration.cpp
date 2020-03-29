
#include "ProgramConfiguration.hpp"
#include "data_structure/atom.hpp"
#include "utils/common.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/dll.hpp>
#include <boost/filesystem.hpp>
#include <pugixml.hpp>

bool ProgramConfiguration::try_load_path(const boost::filesystem::path &path) {
    LOG("try load ", path.native());

    if (boost::filesystem::exists(path) and boost::filesystem::is_regular_file(path)) {
        load_config(path.string());
        LOG(" ...OK\n");
        return true;
    }
    LOG(" ...Fail\n");
    return false;
}

ProgramConfiguration::ProgramConfiguration() {
    const boost::filesystem::path name = ".cac-ana.xml";

    for (auto path = boost::filesystem::current_path();;) {
        if (try_load_path(path / name))
            return;
        auto parent_path = path.parent_path();
        if (parent_path.parent_path() == parent_path)
            break;

        path = parent_path;
    }
    if (try_load_path(std::getenv("HOME") / name))
        return;
    if (try_load_path(boost::dll::program_location().parent_path() / name))
        return;

    LOG("Default configuration is used\n");
}

ProgramConfiguration::ProgramConfiguration(const boost::filesystem::path &config_file) { load_config(config_file); }

void ProgramConfiguration::load_config(const boost::filesystem::path &config_file) {

    pugi::xml_document doc;

    if (auto result = doc.load_file(config_file.c_str()); result) {

        LOG("  load");

        const auto &root = doc.child("cac-ana");
        vector_size = root.child("vector-size").attribute("value").as_ullong();
        nthreads = root.child("nthreads").attribute("value").as_int();

        for (const auto &mask : root.children("macro")) {
            macro_mask.emplace_back(boost::trim_copy(std::string(mask.attribute("name").as_string())),
                                    parse_atoms(mask.attribute("mask").as_string(), true));
        }
    }
}