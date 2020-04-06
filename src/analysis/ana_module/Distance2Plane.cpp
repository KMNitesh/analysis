

#include "Distance2Plane.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include "utils/std.hpp"

Distance2Plane::Distance2Plane() { enable_outfile = true; }

void Distance2Plane::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        for (std::size_t i = 0; i < 3; ++i) {
            if (is_match(atom, plane_marks[i])) {
                plane_atoms[i] = atom;
                return;
            }
        }
        for (std::size_t i = 0; i < outplane_marks.size(); ++i) {
            if (is_match(atom, outplane_marks[i])) {
                outplane_atoms[i] = atom;
                return;
            }
        }
    });

    mol = PBCUtils::calculate_intermol(join(plane_atoms, outplane_atoms), frame);
}

void Distance2Plane::process([[maybe_unused]] std::shared_ptr<Frame> &frame) {
    PBCUtils::move(mol, frame);
    Point p1 = plane_atoms[0]->getCoordinate();
    Point p2 = plane_atoms[1]->getCoordinate();
    Point p3 = plane_atoms[2]->getCoordinate();

    auto parameter = get_panel(p1, p2, p3);

    std::tuple<double, double, double> middle{};
    for (auto &atom : outplane_atoms) {
        middle += atom->getCoordinate();
    }
    middle /= outplane_atoms.size();

    auto dist = dis_pt2panel(middle, parameter);
    acc(dist);

    distances.push_back(dist);
}

void Distance2Plane::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';

    os << "# plane atoms : \n";
    for (const auto &item : plane_marks | boost::adaptors::indexed(0)) {
        os << "#  atom [" << item.index() << "] = " << item.value() << '\n';
    }

    os << "\n# outplane atoms : \n";
    for (const auto &item : outplane_marks | boost::adaptors::indexed(0)) {
        os << "#  atom [" << item.index() << "] = " << item.value() << '\n';
    }
    os << '\n';
    os << "# mean : " << boost::accumulators::mean(acc) << '\n';
    os << "# standard deviation : " << std::sqrt(boost::accumulators::variance(acc)) << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("%15s %15s\n") % "Frame" % "Distance(Ang)";

    const auto fmt = boost::format("%15.3f %15.3f\n");
    for (const auto &element : distances | boost::adaptors::indexed(1)) {
        os << boost::format(fmt) % element.index() % element.value();
    }

    os << std::string(50, '#') << '\n';
}

void Distance2Plane::readInfo() {
    select1group(plane_marks[0], "Enter atom1 for plane > ");
    select1group(plane_marks[1], "Enter atom2 for plane > ");
    select1group(plane_marks[2], "Enter atom3 for plane > ");

    auto num = choose(1, "Enter number of outplanar points > ");
    outplane_atoms.resize(num);
    outplane_marks.resize(num);

    for (int i = 0; i < num; ++i) {
        select1group(outplane_marks[i], "Enter atom(" + std::to_string(i + 1) + ") for plane > ");
    }
}

std::tuple<double, double, double, double> Distance2Plane::get_panel(Point p1, Point p2, Point p3) {
    auto a = ((p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y));

    auto b = ((p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z));

    auto c = ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x));

    auto d = (0 - (a * p1.x + b * p1.y + c * p1.z));

    return {a, b, c, d};
}

double Distance2Plane::dis_pt2panel(Point pt, std::tuple<double, double, double, double> parameter) {
    auto &[a, b, c, d] = parameter;
    return std::abs(a * pt.x + b * pt.y + c * pt.z + d) / std::sqrt(a * a + b * b + c * c);
}
