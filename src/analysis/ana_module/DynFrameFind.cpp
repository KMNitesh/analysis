//
// Created by xiamr on 7/24/19.
//

#include "DynFrameFind.hpp"

#include <boost/phoenix/phoenix.hpp>
#include <boost/phoenix/stl/algorithm.hpp>
#include <boost/range/combine.hpp>

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"

static bool double_tuple_equal(const std::tuple<double, double, double> &v1,
                               const std::tuple<double, double, double> &v2, double eps) {
    return std::abs(std::get<0>(v1) - std::get<0>(v2)) < eps and std::abs(std::get<1>(v1) - std::get<1>(v2)) < eps and
           std::abs(std::get<2>(v1) - std::get<2>(v2)) < eps;
}

static bool double_equal(double v1, double v2, double eps) { return std::abs(v1 - v2) < eps; }

void DynFrameFind::process(std::shared_ptr<Frame> &frame) {
    nframe++;

    const auto &[a_axis, b_axis, c_axis, alpha, beta, gamma] = frame->box.getBoxParameter();
    if (!double_equal(a_axis, reader->getPBCBox().xbox, eps)) return;
    if (!double_equal(b_axis, reader->getPBCBox().ybox, eps)) return;
    if (!double_equal(c_axis, reader->getPBCBox().zbox, eps)) return;
    if (!double_equal(alpha, reader->getPBCBox().alpha, eps)) return;
    if (!double_equal(beta, reader->getPBCBox().beta, eps)) return;
    if (!double_equal(gamma, reader->getPBCBox().gamma, eps)) return;

    std::shared_ptr<Atom> atom;
    std::tuple<double, double, double> coord;
    for (const auto &zipped : boost::combine(frame->atom_list, reader->getAtomicPosition().coordinates)) {
        boost::tie(atom, coord) = zipped;
        if (!double_tuple_equal(atom->getCoordinate(), coord, eps)) return;
    }

    matched_frames.push_back(nframe);
}

void DynFrameFind::print([[maybe_unused]] std::ostream &os) {
    std::cout << title() << '\n';

    std::cout << "Matched frame :";

    for (const auto matched : matched_frames) {
        std::cout << ' ' << matched;
    }
    std::cout << "\n!!!!!!!! Done !!!!!!!\n";
}

void DynFrameFind::readInfo() {
    reader = std::make_shared<TinkerDynReader>(
        std::make_shared<std::fstream>(choose_file("Enter dyn file :").isExist(true).extension("dyn")));
    reader->readContent();
    eps = choose<double>(0, 0.1, "Enter eps [ 0.0001 ] : ", Default(0.0001));
}
