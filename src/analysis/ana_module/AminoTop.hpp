

#ifndef TINKER_AMINOTOP_HPP
#define TINKER_AMINOTOP_HPP

#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include "HBond.hpp"
#include "data_structure/atom.hpp"
#include "utils/std.hpp"

enum class AminoAcidType {
    H3N_Ala,
    H3N_Gly,
    H3N_Ile,
    H3N_Leu,
    H3N_Pro,
    H3N_Val,
    H3N_Phe,
    H3N_Trp,
    H3N_Tyr,
    H3N_Asp,
    H3N_Glu,
    H3N_Arg,
    H3N_His,
    H3N_Lys,
    H3N_Ser,
    H3N_Thr,
    H3N_Cys,
    H3N_Met,
    H3N_Asn,
    H3N_Gln,
    Ala,
    Gly,
    Ile,
    Leu,
    Pro,
    Val,
    Phe,
    Trp,
    Tyr,
    Asp,
    Glu,
    Arg,
    His,
    Lys,
    Ser,
    Thr,
    Cys,
    Met,
    Asn,
    Gln
};

class AminoAcid {
public:
    AminoAcidType type;
    std::map<int, std::string> atom_no_map;
    int sequence_no;
};

extern boost::bimap<boost::bimaps::set_of<AminoAcidType>, boost::bimaps::set_of<std::string>> aminotype_str_bimap;

class AminoTop {
public:
    class AminoItem {
    public:
        int no;
        Symbol symbol;
        std::list<int> linked_atom_nos;
        std::shared_ptr<Atom> atom;
        std::string H_;
    };

    std::map<int, std::shared_ptr<AminoItem>> topmap;

    void atom_null() {
        for (auto &item : topmap) {
            item.second->atom.reset();
        }
    }

    AminoAcidType type;

    void readTop(const std::string &filename);
};

#endif  // TINKER_AMINOTOP_HPP
