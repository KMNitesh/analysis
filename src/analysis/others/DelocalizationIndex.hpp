#ifndef TINKER_DELOCALIZATIONINDEX_HPP
#define TINKER_DELOCALIZATIONINDEX_HPP

#include <string_view>
#include <vector>

class DelocalizationIndex {
public:
    [[nodiscard]] static std::string_view title() { return "Delocalization Index"; }

    static void process();

    static void process_interactive();

private:
    static void process(const std::string &file, const std::vector<std::pair<int, int>> &bonds);
};

#endif  // TINKER_DELOCALIZATIONINDEX_HPP
