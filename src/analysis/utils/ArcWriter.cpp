#include <iomanip>
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "ArcWriter.hpp"

void ArcWriter::open(const std::string &filename) {
    ofs.open(filename);
}

void ArcWriter::close() {
    ofs.close();
}

void ArcWriter::write(const std::shared_ptr<Frame> &frame) {
    ofs << std::setw(6) << frame->atom_list.size() << "  " << frame->title << '\n';
    ofs << std::fixed << std::setprecision(6);
    if (frame->enable_bound) {
        const auto &[a_axis, b_axis, c_axis, alpha, beta, gamma] = frame->box.getBoxParameter();
        ofs << std::setw(13) << a_axis
            << std::setw(12) << b_axis
            << std::setw(12) << c_axis
            << std::setw(12) << alpha
            << std::setw(12) << beta
            << std::setw(12) << gamma << '\n';
    }

    for (auto &atom : frame->atom_list) {

        ofs << std::setw(6) << atom->seq
            << ' '
            << std::left << std::setw(3) << atom->atom_name << std::right
            << std::setw(12) << atom->x
            << std::setw(12) << atom->y
            << std::setw(12) << atom->z
            << std::setw(6) << atom->typ;
        for (auto i : atom->con_list) ofs << std::setw(6) << i;
        ofs << '\n';
    }
}
