//
// Created by xiamr on 6/14/19.
//

#include "SearchInteractionResidue.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

using namespace std;

SearchInteractionResidue::SearchInteractionResidue() { enable_outfile = true; }

void SearchInteractionResidue::process(std::shared_ptr<Frame> &frame) {

    std::unordered_set<std::string> residue_set; // resname:no

    for (auto &atomA : group1) {
        for (auto &atomB : group2) {
            if (atom_distance(atomA, atomB, frame) <= cutoff) {
                residue_set.insert(
                        atomB->residue_name.get() + "-" + boost::lexical_cast<std::string>(atomB->residue_num.get()));
            }
        }
    }
    interaction_residues.push_back(residue_set);
    total_frames++;
}


void SearchInteractionResidue::print(std::ostream &os) {
    os << "************************************************\n";
    os << "*****" << SearchInteractionResidue::title() << " ****\n";

    os << "First Group  " << ids1 << " Second Group  " << ids2 << '\n';
    os << "cutoff :" << cutoff << " Ang\n";

    os << "************************************************\n";


    using ResItem = struct {
        std::string name;
        int count;
    };

    std::unordered_map<std::string, ResItem *> map;
    for (auto &set : interaction_residues) {
        for (auto &item : set) {
            auto it = map.find(item);
            if (it == map.end()) {
                map[item] = new ResItem{item, 1};
            } else {
                it->second->count++;
            }
        }
    }

    std::vector<ResItem *> itemVec;
    for (auto &item : map) {
        itemVec.push_back(item.second);
    }

    std::sort(itemVec.begin(), itemVec.end(), [](ResItem *i1, ResItem *i2) { return i1->count > i2->count; });

    std::size_t nframe = 1;
    os << boost::format("%10s") % "name";
    for (auto &item : itemVec) {
        os << boost::format("%10s") % item->name;
    }
    os << boost::format("\n%10s") % "Freq";;
    for (auto &item : itemVec) {
        os << boost::format("%9.1f%%") % (item->count * 100.0 / total_frames);
    }
    os << '\n';
    for (auto &set : interaction_residues) {
        os << boost::format("%10d") % nframe;
        int index = 1;
        for (auto &item : itemVec) {
            os << boost::format("%10d") % (style == OutputStyle::NUMBER ?
                                           (set.count(item->name) ? index : 0) : set.count(item->name));
            index++;
        }
        os << '\n';
        nframe++;
    }

    os << "************************************************" << endl;


}

void SearchInteractionResidue::readInfo() {
    std::cout << "The output residues is in the second group\n";
    Atom::select2group(ids1, ids2);
    cutoff = choose(0.0, std::numeric_limits<double>::max(), "cutoff [Ang]:");
    std::cout << "(0) bool style\n(1) number style\n";
    style = static_cast<OutputStyle>(choose(0, 1, "which style ? [ 0 ] ", Default(0)));
}

void SearchInteractionResidue::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
                      if (Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
                  });
}