//
// Created by xiamr on 6/14/19.
//

#include "Trajconv.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"
#include <boost/range/algorithm.hpp>

Trajconv::Trajconv(std::shared_ptr<TrajectoryWriterFactoryInterface> factory)
    : factory(std::move(factory)), pbc_utils(std::make_shared<PBCUtils>()) {}

void Trajconv::process(std::shared_ptr<Frame> &frame) {
    pbc_utils->doPBC(pbc_type, mask, frame);

    for (auto &[name, w] : writers) {
        w->write(frame, atoms_for_writetraj);
    }
    step++;
}

void Trajconv::print(std::ostream &) { CleanUp(); }

void Trajconv::CleanUp() {
    for (auto &[name, w] : writers) {
        w->close();
    }
}

void Trajconv::readInfo() {
    inputOutputFiles();
    selectPBCMode();
}

void Trajconv::inputOutputFiles(std::istream &in, std::ostream &out) {
    for (;;) {
        std::string filename = choose_file("output file [empty for next]: ", in, out).isExist(false).can_empty(true);
        if (filename.empty()) {
            if (!writers.empty()) {
                break;
            }
        }
        try {
            writers.emplace_back(filename, factory->make_instance(getFileType(filename)));
        } catch (std::exception &) {
            out << "ERROR !!  wrong type of target trajectory file (" << filename << ")\n";
        }
    }
}

void Trajconv::selectPBCMode() {
    std::cout << "PBC transform option\n";
    std::cout << "(0) Do nothing\n";
    std::cout << "(1) Make atom i as center\n";
    std::cout << "(2) Make molecule i as center\n";
    std::cout << "(3) Make atom group as center\n";
    std::cout << "(4) all atoms into box\n";

    while (true) {
        int option = choose(0, 4, "Choose : ");
        switch (option) {
        case 0:
            pbc_type = PBCType::None;
            break;
        case 1:
            pbc_type = PBCType::OneAtom;
            Atom::select1group(mask, "Please enter atom mask : ");
            break;
        case 2:
            pbc_type = PBCType::OneMol;
            Atom::select1group(mask, "Please enter one atom mask that the molecule include: ");
            break;
        case 3:
            pbc_type = PBCType::AtomGroup;
            Atom::select1group(mask, "Please enter atom group: ");
            break;
        case 4:
            pbc_type = PBCType::AllIntoBox;
            break;
        default:
            std::cerr << "option not found !\n";
            continue;
        }
        break;
    }
    Atom::select1group(mask_for_writetraj, "Enter mask for output > ", true);
}

void Trajconv::fastConvertTo(std::string target) noexcept(false) {
    boost::trim(target);
    if (target.empty()) {
        throw std::runtime_error("ERROR !! empty target trajectory file");
    }
    pbc_type = PBCType::None;

    try {
        writers.emplace_back(target, factory->make_instance(getFileType(target)));
    } catch (std::exception &) {
        throw std::runtime_error("ERROR !!  wrong type of target trajectory file (" + target + ")");
    }
}

void Trajconv::processFirstFrame(std::shared_ptr<Frame> &frame) {
    for (auto &[name, w] : writers) {
        w->open(name);
    }
    if (!Atom::isBlank(mask_for_writetraj)) {
        boost::for_each(frame->atom_list, [this](const std::shared_ptr<Atom> &atom) {
            if (Atom::is_match(atom, mask_for_writetraj)) {
                atoms_for_writetraj.push_back(atom);
            }
        });
    } else {
        atoms_for_writetraj = frame->atom_list;
    }
}
