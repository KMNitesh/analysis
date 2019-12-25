//
// Created by xiamr on 6/14/19.
//

#include "utils/common.hpp"
#include "Trajconv.hpp"

#include "data_structure/molecule.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/atom.hpp"
#include "utils/PBCUtils.hpp"

Trajconv::Trajconv(std::shared_ptr<TrajectoryWriterFactoryInterface> factory) :
        factory(std::move(factory)),
        pbc_utils(std::make_shared<PBCUtils>()) {}


void Trajconv::process(std::shared_ptr<Frame> &frame) {
    pbc_utils->doPBC(pbc_type, num, frame);

    for (auto &[name, w] : writers) {
        w->write(frame);
    }
    step++;
}


void Trajconv::print(std::ostream &) {
    CleanUp();
}

void Trajconv::CleanUp() {
    for (auto &[name, w]: writers) {
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

void Trajconv::initPBC(PBCType pbc_mode, int num) {
    pbc_type = pbc_mode;
    this->num = num;
}

void Trajconv::selectPBCMode() {
    std::cout << "PBC transform option\n";
    std::cout << "(0) Do nothing\n";
    std::cout << "(1) Make atom i as center\n";
    std::cout << "(2) Make molecule i as center\n";

    while (true) {
        int option = choose(0, 2, "Choose : ");
        switch (option) {
            case 0:
                pbc_type = PBCType::None;
                break;
            case 1:
                pbc_type = PBCType::OneAtom;
                num = choose(1, std::numeric_limits<int>::max(), "Please enter the atom NO. : ");
                break;
            case 2:
                pbc_type = PBCType::OneMol;
                num = choose(1, std::numeric_limits<int>::max(),
                             "Please enter one atom NO. that the molecule include: ");
                break;
            default:
                std::cerr << "option not found !\n";
                continue;
        }
        break;
    }
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

void Trajconv::processFirstFrame(std::shared_ptr<Frame> &) {
    for (auto &[name, w] : writers) {
        w->open(name);
    }
}


