//
// Created by xiamr on 7/24/19.
//

#ifndef TINKER_DYNFRAMEFIND_HPP
#define TINKER_DYNFRAMEFIND_HPP

#include "common.hpp"
#include "AbstractAnalysis.hpp"
#include "TinkerDynReader.hpp"

class DynFrameFind : public AbstractAnalysis {
public:
    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Find Frame Index base on Atomic Position of dyn file"; }

protected:

    std::shared_ptr<TinkerDynReader> reader;
    double eps;

    int nframe = 0;
    std::list<int> matched_frames;
};


#endif //TINKER_DYNFRAMEFIND_HPP
