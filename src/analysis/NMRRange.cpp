//
// Created by xiamr on 6/14/19.
//

#include "NMRRange.hpp"

#include "frame.hpp"

using namespace std;

void NMRRange::process(std::shared_ptr<Frame> &frame) {
    if (first_frame) {
        loadTop();
        recognize_amino_acid(frame);
        first_frame = false;
        vector<int> need_to_calc_list;
        for (auto &aminoacid : amino_acid_list) {
            for (auto &item : aminoacid->atom_no_map)
                need_to_calc_list.push_back(item.first);
        }
        for (size_t i = 0; i < need_to_calc_list.size() - 1; i++) {
            for (size_t j = i + 1; j < need_to_calc_list.size(); j++) {
                dist_range_map[make_pair(need_to_calc_list[i], need_to_calc_list[j])] = list<double>();
                dist_range_map[make_pair(need_to_calc_list[i], need_to_calc_list[j])].push_back(
                        atom_distance(frame->atom_map[need_to_calc_list[i]], frame->atom_map[need_to_calc_list[j]],
                                      frame));
            }
        }
    } else {
        for (auto &item : dist_range_map) {
            item.second.push_back(
                    atom_distance(frame->atom_map[item.first.first], frame->atom_map[item.first.second], frame));
        }
    }

}


void NMRRange::print(std::ostream &os) {

    for (auto &item : dist_range_map) {

        double sum = 0.0;
        int count = 0;
        for (auto v : item.second) {
            sum += pow(v, -6);
            count++;
        }
        double avg = sum / count;
        double e = -1.0 / 6;
        double value = pow(avg, e);
        os << name_map[item.first.first] << " <->\t"
           << name_map[item.first.second] << "\t" << value << endl;
    }
}

void NMRRange::readInfo() {
    enable_forcefield = true;

}


void NMRRange::recognize_amino_acid(std::shared_ptr<Frame> &frame) {
    cout << "First Frame, Analyze..." << endl;
    int amino_seq_no = 0;
    for (auto atom : frame->atom_list) {
        if (which(atom) == Symbol::Nitrogen) {
            for (auto &amino_top : amino_top_list) {
                list<shared_ptr<Atom>> child_atom_list;
                list<shared_ptr<AminoTop::AminoItem>> child_item_list;
                bool match = recognize_walk(atom, amino_top.topmap[1], amino_top, frame,
                                            child_atom_list, child_item_list);


                if (match) {
                    amino_seq_no++;
                    cout << aminotype_str_bimap.left.find(amino_top.type)->second << "   " << endl;

                    auto new_amino = make_shared<AminoAcid>();
                    new_amino->type = amino_top.type;
                    new_amino->sequence_no = amino_seq_no;

                    for (auto &item :amino_top.topmap) {
                        if (item.second->symbol == Symbol::Hydrogen) {
                            new_amino->atom_no_map[item.second->atom->seq] = item.second->H_;

                            name_map[item.second->atom->seq] = to_string(item.second->atom->seq) + ":"
                                                               + aminotype_str_bimap.left.find(amino_top.type)->second
                                                               + ":" + to_string(amino_seq_no) + " " + item.second->H_;
                        }
                        cout << item.first << " " << item.second->H_ <<
                             "," << item.second->atom->seq << " " << item.second->atom->atom_name << endl;
                        //if (item.second->symbol == Symbol::X)
                        item.second->atom->mark = false;
                    }
                    amino_acid_list.push_back(new_amino);
                    amino_top.atom_null();
                    break;
                }


            }
        }
    }
    cout << endl;

}

bool NMRRange::recognize_walk(shared_ptr<Atom> atom, shared_ptr<AminoTop::AminoItem> item,
                              AminoTop &top, shared_ptr<Frame> &frame,
                              list<shared_ptr<Atom>> &atom_list,
                              list<shared_ptr<AminoTop::AminoItem>> &item_list) {
    if (item->symbol == Symbol::X) {
        atom->mark = true;
        item->atom = atom;
        atom_list.push_back(atom);
        item_list.push_back(item);
        return true;
    }
    if (which(atom) != item->symbol) {
        return false;
    }
    if (atom->con_list.size() != item->linked_atom_nos.size()) {

        return false;
    }
    atom->mark = true;
    item->atom = atom;

    bool match = true;
    list<shared_ptr<Atom>> child_atom_list;
    list<shared_ptr<AminoTop::AminoItem>> child_item_list;
    for (int i : atom->con_list) {
        auto next_atom = frame->atom_map[i];
        if (next_atom->mark) continue;
        bool ok = false;

        for (int j : item->linked_atom_nos) {
            auto next_item = top.topmap[j];
            if (next_item->atom) continue;
            ok = recognize_walk(next_atom, next_item, top, frame, child_atom_list, child_item_list);
            if (ok) break;
        }

        match = ok && match;
        if (!match) break;
    }
    if (match) {
        atom_list.push_back(atom);
        item_list.push_back(item);
        for (auto &child_atom: child_atom_list) {
            atom_list.push_back(child_atom);
        }
        for (auto &child_item : child_item_list) {
            item_list.push_back(child_item);
        }

        return true;
    }

    atom->mark = false;
    item->atom.reset();

    for (auto &child_atom: child_atom_list) {
        child_atom->mark = false;
    }
    for (auto &child_item : child_item_list) {
        child_item->atom.reset();
    }

    return false;
}