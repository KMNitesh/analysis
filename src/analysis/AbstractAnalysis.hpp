//
// Created by xiamr on 6/14/19.
//

//  This Pure Virtrial Class is the interface for All trajectory analysis method

#ifndef TINKER_ABSTRACTANALYSIS_HPP
#define TINKER_ABSTRACTANALYSIS_HPP

#include "std.hpp"
#include <boost/core/ignore_unused.hpp>

class Frame;

class AbstractAnalysis {
public:
    virtual void processFirstFrame(std::shared_ptr<Frame> &frame) { boost::ignore_unused(frame); };

    virtual void process(std::shared_ptr<Frame> &frame) = 0;

    virtual void print(std::ostream &os) = 0;

    virtual void readInfo() = 0;

    virtual std::string getOutfileName() { return outfilename; }

    virtual std::string description() { return "Abstract Class"; }

    static std::string title() { return "Abstract Class"; }

    void do_parallel_while(std::function<std::shared_ptr<Frame>()> func) {
        do_parallel_while_impl(std::move(func));
    }

    bool enable_parralle_while() {
        return enable_paralel_while_impl();
    }

    virtual ~AbstractAnalysis() = default;

protected:
    virtual void do_parallel_while_impl(std::function<std::shared_ptr<Frame>()> func) { boost::ignore_unused(func); };

    virtual bool enable_paralel_while_impl() { return false; }

    std::string outfilename;
};


#endif //TINKER_ABSTRACTANALYSIS_HPP
