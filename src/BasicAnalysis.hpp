//
// Created by xiamr on 6/14/19.
//

//  This Pure Virtrial Class is the interface for All trajectory analysis method

#ifndef TINKER_BASICANALYSIS_HPP
#define TINKER_BASICANALYSIS_HPP

#include <memory>
#include <string>
#include <ostream>

class Frame;

class BasicAnalysis {
public:
    virtual void processFirstFrame(std::shared_ptr<Frame> &/*frame*/) {};

    virtual void process(std::shared_ptr<Frame> &frame) = 0;

    virtual void print(std::ostream &os) = 0;

    virtual void readInfo() = 0;

    static const std::string title() { return "Base Class"; }


    virtual ~BasicAnalysis() = default;
};


#endif //TINKER_BASICANALYSIS_HPP
