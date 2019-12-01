//
// Created by xiamr on 6/14/19.
//

#include "FirstCoordExchangeSearch.hpp"
#include "frame.hpp"
#include "common.hpp"

FirstCoordExchangeSearch::FirstCoordExchangeSearch() {
    enable_outfile = true;
}

void FirstCoordExchangeSearch::process(std::shared_ptr<Frame> &frame) {

    step++;

    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            if (step == 1) {
                State state;
                state.inner = atom_distance(atom1, atom2, frame) <= this->dist_cutoff;
                state_machine[atom2->seq] = state;
                if (state.inner) init_seq_in_shell.insert(atom2->seq);
            } else {
                auto &state = state_machine[atom2->seq];
                if (state.inner) {
                    if (atom_distance(atom1, atom2, frame) >= this->dist_cutoff + tol_dist) {
                        state.inner = false;
                        ExchangeItem item;
                        item.seq = atom2->seq;
                        item.direction = Direction::OUT;
                        item.exchange_frame = step;
                        exchange_list.push_back(item);
                    }
                } else {
                    if (atom_distance(atom1, atom2, frame) <= this->dist_cutoff - tol_dist) {
                        state.inner = true;
                        ExchangeItem item;
                        item.seq = atom2->seq;
                        item.direction = Direction::IN;
                        item.exchange_frame = step;
                        exchange_list.push_back(item);
                    }
                }
            }
        }

    }
}

void FirstCoordExchangeSearch::print(std::ostream &os) {
    os << "***************************" << '\n';
    os << "***** Exchange Search *****" << '\n';
    os << "type 1 :" << ids1 << '\n';
    os << "type 2 :" << ids2 << '\n';

    os << "cutoff :" << dist_cutoff << '\n';
    os << "tol dist :" << tol_dist << '\n';
    os << "time_cutoff :" << time_cutoff << '\n';
    os << "***************************" << '\n';
    for (int seq : init_seq_in_shell) {
        os << "  " << seq;
    }
    os << "\n***************************" << '\n';

    os << "** seq **  direction ***** exchange frame *****" << '\n';
    for (auto it = exchange_list.begin(); it != exchange_list.end(); it++) {
        os << boost::format("%10d%6s%15d   !   ")
              % it->seq
              % (it->direction == Direction::IN ? "IN" : "OUT")
              % it->exchange_frame;
        if (it->direction == Direction::IN) {
            init_seq_in_shell.insert(it->seq);
        } else {
            init_seq_in_shell.erase(it->seq);
        }
        for (int seq : init_seq_in_shell) {
            os << "  " << seq;
        }
        os << '\n';
    }
    os << "***************************\n";


}

void FirstCoordExchangeSearch::readInfo() {
    Atom::select2group(ids1, ids2);
    dist_cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
    tol_dist = choose(0.0, std::numeric_limits<double>::max(), "Please enter tol dist:");
    time_cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter timecutoff:");
}

void FirstCoordExchangeSearch::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](std::shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}