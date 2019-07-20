//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_FIRSTCOORDEXCHANGESEARCH_HPP
#define TINKER_FIRSTCOORDEXCHANGESEARCH_HPP

#include <memory>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <list>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "frame.hpp"
#include "atom.hpp"

class Frame;

class FirstCoordExchangeSearch : public BasicAnalysis {

public:
    FirstCoordExchangeSearch() {
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "Water exchange analysis";
    }


private:

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    double dist_cutoff, tol_dist, time_cutoff;
    enum class Direction {
        IN, OUT
    };
    typedef struct {
        Direction direction;
        int seq;
        int exchange_frame;
    } ExchangeItem;

    std::list<ExchangeItem> exchange_list;

    int step = 0;

    typedef struct {
        bool inner;
    } State;
    std::map<int, State> state_machine;

    std::set<int> init_seq_in_shell;

};


#endif //TINKER_FIRSTCOORDEXCHANGESEARCH_HPP
