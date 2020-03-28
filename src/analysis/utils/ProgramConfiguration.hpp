
#ifndef TINKER_PROGRAMCONFIGURATION_HPP
#define TINKER_PROGRAMCONFIGURATION_HPP

#include "data_structure/atom.hpp"
#include "utils/std.hpp"

namespace boost::filesystem {
class path;
}

class ProgramConfiguration {
public:
    ProgramConfiguration();

    ProgramConfiguration(const boost::filesystem::path &config_file);

    void load_config(const boost::filesystem::path &config_file);

    auto get_vector_size() { return vector_size; }

    auto get_nthreads() { return nthreads; }

    const auto &get_macro_mask() { return macro_mask; }

private:
    bool try_load_path(const boost::filesystem::path &path);

    std::size_t vector_size = 0;
    int nthreads = 0;

    std::vector<std::pair<std::string, AmberMask>> macro_mask;
};

#endif // TINKER_PROGRAMCONFIGURATION_HPP