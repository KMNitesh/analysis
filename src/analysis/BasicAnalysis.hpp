//
// Created by xiamr on 6/14/19.
//

//  This Pure Virtrial Class is the interface for All trajectory analysis method

#ifndef TINKER_BASICANALYSIS_HPP
#define TINKER_BASICANALYSIS_HPP

#include <utility>

#include "std.hpp"

class Frame;

class BasicAnalysis {
protected:
    std::string outfilename;
public:
    virtual void processFirstFrame(std::shared_ptr<Frame> &frame) {};

    virtual void process(std::shared_ptr<Frame> &frame) = 0;

    virtual void print(std::ostream &os) = 0;

    virtual void readInfo() = 0;

    virtual std::string getOutfileName() { return outfilename; }

    virtual std::string description() { return "Base Class"; }

    static const std::string title() { return "Base Class"; }

    void do_parallel_while(std::function<std::shared_ptr<Frame>()> func) {
        do_parallel_while_impl(std::move(func));
    }

    bool enable_parralle_while() {
        return enable_paralel_while_impl();
    }

    virtual ~BasicAnalysis() = default;

protected:
    virtual void do_parallel_while_impl(std::function<std::shared_ptr<Frame>()> func) {};

    virtual bool enable_paralel_while_impl() { return false; }
};


#endif //TINKER_BASICANALYSIS_HPP
