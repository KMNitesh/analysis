//
// Created by xiamr on 6/14/19.
//

#include "AminoTop.hpp"

#include <boost/assign.hpp>

#include "utils/common.hpp"

boost::bimap<boost::bimaps::set_of<AminoAcidType>, boost::bimaps::set_of<std::string>> aminotype_str_bimap =
    boost::assign::list_of<
        boost::bimap<boost::bimaps::set_of<AminoAcidType>, boost::bimaps::set_of<std::string>>::relation>(
        AminoAcidType::H3N_Ala, "H3N_Ala")(AminoAcidType::H3N_Gly, "H3N_Gly")(AminoAcidType::H3N_Pro, "H3N_Pro")(
        AminoAcidType::H3N_Trp, "H3N_Trp")(AminoAcidType::H3N_Asp, "H3N_Asp")(AminoAcidType::H3N_Glu, "H3N_Glu")(
        AminoAcidType::H3N_Arg, "H3N_Arg")

        (AminoAcidType::Ala, "Ala")(AminoAcidType::Gly, "Gly")(AminoAcidType::Pro, "Pro")(AminoAcidType::Trp, "Trp")(
            AminoAcidType::Asp, "Asp")(AminoAcidType::Glu, "Glu")(AminoAcidType::Arg, "Arg");

void AminoTop::readTop(const std::string &filename) {
    std::fstream f(filename);
    if (!f) {
        throw std::runtime_error("Can not open file : " + filename);
    }
    std::string line;
    while (true) {
        std::getline(f, line);
        auto field = split(line);
        if (field.empty())
            break;
        auto item = std::make_shared<AminoItem>();
        if (field.size() < 3) {
            throw std::runtime_error("file <" + filename + "> content sytanx error \n");
        }
        item->no = stoi(field[0]);
        auto sym = field[1];
        if (sym == "H")
            item->symbol = Symbol::Hydrogen;
        else if (sym == "C")
            item->symbol = Symbol::Carbon;
        else if (sym == "N")
            item->symbol = Symbol::Nitrogen;
        else if (sym == "O")
            item->symbol = Symbol::Oxygen;
        else if (sym == "X")
            item->symbol = Symbol::X;
        else {
            throw std::runtime_error("file <" + filename + "> has unrecognized element symbol " + sym + '\n');
        }
        item->H_ = field[2];
        for (unsigned int i = 3; i < field.size(); i++) {
            item->linked_atom_nos.push_back(stoi(field[i]));
        }
        topmap[item->no] = item;
    }
}