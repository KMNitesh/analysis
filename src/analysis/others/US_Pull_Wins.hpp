
#ifndef TINKER_OTHER_US_PULL_WINS
#define TINKER_OTHER_US_PULL_WINS

#include "utils/std.hpp"

class US_Pull_Wins {

public:
    [[nodiscard]] static std::string_view title() { return "Umbrella Sampling Windows"; }

    static void process();
};

#endif // TINKER_OTHER_US_PULL_WINS