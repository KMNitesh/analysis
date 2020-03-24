
#include <boost/format.hpp>
#include <fstream>
#include <map>
#include <memory>
#include <string>

namespace gmx {

#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/utility/smalloc.h"

}  // namespace gmx

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "gro_writer.hpp"

void GROWriter::open(const std::string &filename) { os->open(filename, std::ios::out); }

void GROWriter::close() { os->close(); }

void GROWriter::write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) {
    *os << frame->title << '\n';
    *os << frame->atom_list.size() << '\n';
    const boost::format fmt("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n");
    for (auto &atom : atoms) {
        const auto &mol = atom->molecule.lock();
        *os << boost::format(fmt) %
                   ((atom->residue_num.has_value() ? atom->residue_num.get() : mol->sequence) % 100000) %
                   (atom->residue_name.has_value() ? atom->residue_name.get() : std::to_string(mol->sequence)) %
                   atom->atom_name % (atom->seq % 100000) % (atom->x / 10.0) % (atom->y / 10.0) % (atom->z / 10.0);
    }
    if (frame->enable_bound) {
        gmx::matrix box;
        frame->box.getBoxParameter(box);
        if (frame->box.get_box_type() == PBCBox::Type::orthogonal)
            *os << boost::format(" %9.5f %9.5f %9.5f\n") % box[0][0] % box[1][1] % box[2][2];
        else
            *os << boost::format(" %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n") % box[0][0] % box[1][1] %
                       box[2][2] % box[0][1] % box[0][2] % box[1][0] % box[1][2] % box[2][0] % box[2][1];
    }
}
