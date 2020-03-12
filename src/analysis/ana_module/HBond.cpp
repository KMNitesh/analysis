//
// Created by xiamr on 6/14/19.
//

#include "HBond.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

using namespace std;

HBond::HBond() {
    enable_outfile = true;
    enable_forcefield = true;
}

Symbol which(const std::shared_ptr<Atom> &atom) {
    double mass = forcefield.find_mass(atom);
    if (mass < 2.0)
        return Symbol::Hydrogen;
    else if (mass > 11.5 and mass < 13.0)
        return Symbol::Carbon;
    else if (mass > 13.0 and mass < 15.0)
        return Symbol::Nitrogen;
    else if (mass >= 15.0 and mass < 17.0)
        return Symbol::Oxygen;
    else if (mass >= 30.0 and mass < 31.5)
        return Symbol::Phosphorus;
    else if (mass >= 31.5 and mass < 33.0)
        return Symbol::Sulfur;
    else if (mass >= 22.5 and mass < 23.5)
        return Symbol::Sodium;
    else
        return Symbol::Unknown;
}

void HBond::print(std::ostream &os) {
    os << "***************************" << endl;
    os << "****** Hydrogen Bond ******" << endl;
    os << "TYPE:";
    switch (mode) {
        case Selector::Both:
            os << "Both" << endl;
            os << "SET1:" << ids1 << endl;
            break;
        case Selector::Donor:
            os << "Donor" << endl;
            os << "SET1:" << ids1 << endl;
            os << "SET2:" << ids2 << endl;
            break;
        case Selector::Acceptor:
            os << "Acceptor" << endl;
            os << "SET1:" << ids1 << endl;
            os << "SET2:" << ids2 << endl;
            break;
    }
    switch (hbond_type) {
        case HBondType::VMDVerion:
            os << "HBond criteria : VMD version" << endl;
            break;
        case HBondType::GMXVersion:
            os << "HBond criteria : GMX version" << endl;
            break;
    }
    os << "distance cutoff : " << this->donor_acceptor_dist_cutoff << endl;
    os << "angle cutoff : " << this->angle_cutoff << endl;
    os << "***************************" << endl;
    for (auto cyc : range(1, steps + 1)) {
        os << cyc << "   " << hbonds[cyc] << endl;
    }
    os << "***************************" << endl;
}

void HBond::process(std::shared_ptr<Frame> &frame) {
    steps++;
    hbonds[steps] = 0;
    switch (mode) {
        case Selector::Both:
            Selector_Both(frame);
            break;
        case Selector::Donor:
            Selector_Donor(frame);
            break;
        case Selector::Acceptor:
            Selector_Acceptor(frame);
            break;
    }
}

void HBond::readInfo() {
    Atom::select1group(ids1, "Please enter group:");

    auto input_line = input("Which selector: [(1)Acceptor | (2)Donor | (3)Both]:");
    auto mode_no = stoi(input_line);
    switch (mode_no) {
        case 1:
            this->mode = Selector::Acceptor;
            break;
        case 2:
            this->mode = Selector::Donor;
            break;
        case 3:
            this->mode = Selector::Both;
            break;
        default:
            cerr << "Error type" << endl;
            exit(1);
    }
    if (this->mode not_eq Selector::Both) {
        Atom::select1group(ids2, "Please enter group2:");
    }
    cout << "HBond criteria:\n(1)VMD version\n(2)GMX version\n";
    switch (choose(1, 2, "choose:")) {
        case 1:
            hbond_type = HBondType::VMDVerion;
            break;
        case 2:
            hbond_type = HBondType::GMXVersion;
            break;
        default:
            cerr << "wrong type! \n";
            exit(2);
    }

    this->donor_acceptor_dist_cutoff =
        choose(0.0, static_cast<double>(std::numeric_limits<int>::max()), "Donor-Acceptor Distance:");
    this->angle_cutoff = choose(0.0, static_cast<double>(std::numeric_limits<int>::max()), "Angle cutoff:");
}

void HBond::Selector_Acceptor(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group2) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3 : group1) {
                    // Atom3
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            auto vx1 = atom2->x - atom1->x;
                            auto vy1 = atom2->y - atom1->y;
                            auto vz1 = atom2->z - atom1->z;
                            frame->image(vx1, vy1, vz1);
                            auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                            double vx2, vy2, vz2;
                            switch (hbond_type) {
                                case HBondType::VMDVerion:
                                    // atom2 -> atom3
                                    vx2 = atom3->x - atom2->x;
                                    vy2 = atom3->y - atom2->y;
                                    vz2 = atom3->z - atom2->z;
                                    break;
                                case HBondType::GMXVersion:
                                    // atom1 -> atom3
                                    vx2 = atom3->x - atom1->x;
                                    vy2 = atom3->y - atom1->y;
                                    vz2 = atom3->z - atom1->z;
                                    break;
                            }

                            frame->image(vx2, vy2, vz2);
                            auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
                            auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                            auto ang = radian * std::acos(cosine);
                            if (std::abs(ang) <= this->angle_cutoff) {
                                hbonds[steps]++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::Selector_Donor(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group1) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3 : group2) {
                    // Atom3
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            auto vx1 = atom2->x - atom1->x;
                            auto vy1 = atom2->y - atom1->y;
                            auto vz1 = atom2->z - atom1->z;
                            frame->image(vx1, vy1, vz1);
                            auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                            double vx2, vy2, vz2;
                            switch (hbond_type) {
                                case HBondType::VMDVerion:
                                    // atom2 -> atom3
                                    vx2 = atom3->x - atom2->x;
                                    vy2 = atom3->y - atom2->y;
                                    vz2 = atom3->z - atom2->z;
                                    break;
                                case HBondType::GMXVersion:
                                    // atom1 -> atom3
                                    vx2 = atom3->x - atom1->x;
                                    vy2 = atom3->y - atom1->y;
                                    vz2 = atom3->z - atom1->z;
                                    break;
                            }

                            frame->image(vx2, vy2, vz2);
                            auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);
                            auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                            auto ang = radian * std::acos(cosine);
                            if (std::abs(ang) <= this->angle_cutoff) {
                                hbonds[steps]++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::Selector_Both(std::shared_ptr<Frame> &frame) {
    for (auto &atom2 : group1) {
        // Atom2
        if (which(atom2) == Symbol::Hydrogen) {
            // Atom1
            auto &atom1 = frame->atom_map[atom2->con_list.front()];
            auto atom1_symbol = which(atom1);
            if (atom1_symbol == Symbol::Nitrogen or atom1_symbol == Symbol::Oxygen) {
                for (auto &atom3 : group1) {
                    // atom1 and atom3 can not the same atom
                    if (atom1 == atom3) continue;
                    auto atom3_symbol = which(atom3);
                    if (atom3_symbol == Symbol::Nitrogen or atom3_symbol == Symbol::Oxygen) {
                        auto distance = atom_distance(atom1, atom3, frame);
                        if (distance <= this->donor_acceptor_dist_cutoff) {
                            if (atom1 != atom3 && !atom1->adj(atom3)) {
                                auto vx1 = atom2->x - atom1->x;
                                auto vy1 = atom2->y - atom1->y;
                                auto vz1 = atom2->z - atom1->z;
                                // atom1 -> atom2
                                frame->image(vx1, vy1, vz1);
                                auto len1 = std::sqrt(vx1 * vx1 + vy1 * vy1 + vz1 * vz1);

                                double vx2, vy2, vz2;
                                switch (hbond_type) {
                                    case HBondType::VMDVerion:
                                        // atom2 -> atom3
                                        vx2 = atom3->x - atom2->x;
                                        vy2 = atom3->y - atom2->y;
                                        vz2 = atom3->z - atom2->z;
                                        break;
                                    case HBondType::GMXVersion:
                                        // atom1 -> atom3
                                        vx2 = atom3->x - atom1->x;
                                        vy2 = atom3->y - atom1->y;
                                        vz2 = atom3->z - atom1->z;
                                        break;
                                }

                                frame->image(vx2, vy2, vz2);
                                auto len2 = std::sqrt(vx2 * vx2 + vy2 * vy2 + vz2 * vz2);

                                auto cosine = (vx1 * vx2 + vy1 * vy2 + vz1 * vz2) / (len1 * len2);
                                auto ang = radian * std::acos(cosine);
                                if (std::abs(ang) <= this->angle_cutoff) {
                                    hbonds[steps]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void HBond::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(), [this](shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, this->ids1)) this->group1.insert(atom);
        if (this->mode not_eq Selector::Both && Atom::is_match(atom, this->ids2)) this->group2.insert(atom);
    });
}
