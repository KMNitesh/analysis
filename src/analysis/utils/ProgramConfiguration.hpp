
#ifndef TINKER_PROGRAMCONFIGURATION_HPP
#define TINKER_PROGRAMCONFIGURATION_HPP

#include <cstddef>
#include <string>

class ProgramConfiguration {
public:
    ProgramConfiguration();
    
    ProgramConfiguration(const std::string &config_file);

    void load_config(const std::string &config_file);

    auto get_vector_size() { return vector_size; }

private:
    std::size_t vector_size = 0;
};

#endif // TINKER_PROGRAMCONFIGURATION_HPP