
#include <exception>

#include <boost/range/adaptors.hpp>

#include "CoordinateNumPerFrame.hpp"

#include "data_structure/frame.hpp"

CoordinateNumPerFrame::CoordinateNumPerFrame() { enable_outfile = true; }

void CoordinateNumPerFrame::process(std::shared_ptr<Frame> &frame) {
    int cn_sum = 0;
    for (auto &atom1 : group1) {
        for (auto &atom2 : group2) {
            if (atom_distance(atom1, atom2, frame) <= this->dist_cutoff) {
                cn_sum++;
            }
        }
    }
    cn_list.push_back(cn_sum);
}

void CoordinateNumPerFrame::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# mask for group1 > " << ids1 << '\n';
    os << "# mask for group2 > " << ids2 << '\n';
    os << "# cutoff(Ang) > " << dist_cutoff << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("#%15s %15s\n") % "Frame" % "CN";
    for (const auto &element : cn_list | boost::adaptors::indexed(1)) {
        os << boost::format(" %15d %15d\n") % element.index() % element.value();
    }
    os << std::string(50, '#') << '\n';
}

void CoordinateNumPerFrame::readInfo() {
    select2group(ids1, ids2);
    dist_cutoff = choose(0.0, std::numeric_limits<double>::max(), "Please enter distance cutoff:");
}

void CoordinateNumPerFrame::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(), [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, this->ids1)) this->group1.insert(atom);
        if (is_match(atom, this->ids2)) this->group2.insert(atom);
    });
}

void CoordinateNumPerFrame::setParameters(const AmberMask &M, const AmberMask &L, double cutoff,
                                          const std::string &out) {
    ids1 = M;
    ids2 = L;
    dist_cutoff = cutoff;

    outfilename = out;
    boost::trim(outfilename);
    if (outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}
