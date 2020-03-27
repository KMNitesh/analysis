
#include "ProgramConfiguration.hpp"
#include "utils/common.hpp"
#include <boost/dll.hpp>
#include <boost/filesystem.hpp>
#include <pugixml.hpp>

ProgramConfiguration::ProgramConfiguration() {
    if (auto path = boost::filesystem::path(std::getenv("HOME")) / ".cac-ana";
        boost::filesystem::exists(path) and boost::filesystem::is_regular_file(path)) {
        load_config(path.string());
        return;
    }
    if (auto path = boost::dll::program_location().parent_path() / ".cac-ana";
        boost::filesystem::exists(path) and boost::filesystem::is_regular_file(path)) {
        load_config(path.string());
        return;
    }
}

ProgramConfiguration::ProgramConfiguration(const std::string &config_file) { load_config(config_file); }

void ProgramConfiguration::load_config(const std::string &config_file) {

    pugi::xml_document doc;

    if (auto result = doc.load_file(config_file.c_str()); result) {
        if (verbose_message) {
            std::cout << "load config file : " << config_file << '\n';
        }
        const auto &root = doc.child("cac-ana");
        vector_size = root.child("vector-size").attribute("value").as_ullong();
    }
}