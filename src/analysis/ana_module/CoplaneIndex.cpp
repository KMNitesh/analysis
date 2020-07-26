
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "CoplaneIndex.hpp"

#include "data_structure/frame.hpp"
#include "nlohmann/json.hpp"
#include "utils/common.hpp"
#include "utils/std.hpp"

CoplaneIndex::CoplaneIndex() { enable_outfile = true; }

void CoplaneIndex::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        for (std::size_t i = 0; i < mask_arrays.size(); ++i) {
            for (std::size_t j = 0; j < mask_arrays[i].size(); ++j) {
                if (is_match(atom, mask_arrays[i][j])) {
                    atom_arrays[i][j] = atom;
                    return;
                }
            }
        }
    });
}

void CoplaneIndex::process(std::shared_ptr<Frame> &frame) {
    static std::vector<std::tuple<double, double, double>> plane_normal_vectors;
    plane_normal_vectors.clear();

    for (const auto &plane : atom_arrays) {
        auto atom1_coord = plane[0]->getCoordinate();
        auto atom2_coord = plane[1]->getCoordinate();
        auto atom3_coord = plane[2]->getCoordinate();

        auto c2_c1 = atom2_coord - atom1_coord;
        auto c3_c2 = atom3_coord - atom2_coord;

        frame->image(c2_c1);
        frame->image(c3_c2);

        auto v = cross_multiplication(c2_c1, c3_c2);
        v /= vector_norm(v);
        plane_normal_vectors.push_back(v);
    }

    using namespace boost::accumulators;
    accumulator_set<double, features<tag::mean, tag::variance>> acc;

    for (auto it1 = begin(plane_normal_vectors); it1 != end(plane_normal_vectors) - 1; ++it1) {
        for (auto it2 = it1 + 1; it2 != end(plane_normal_vectors); ++it2) {
            acc(std::abs(dot_multiplication(*it1, *it2)));
        }
    }

    coplaneIndex.emplace_back(mean(acc), std::sqrt(variance(acc)));
}

void CoplaneIndex::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    nlohmann::json json;
    json["title"] = title();

    for (const auto &plane : mask_arrays | boost::adaptors::indexed(1)) {
        nlohmann::json plane_json;
        os << "# plane <" << plane.index() << ">\n";
        for (const auto &item : plane.value() | boost::adaptors::indexed(0)) {
            os << "#  atom [" << item.index() << "] = " << item.value() << '\n';
            plane_json.push_back(to_string(item.value()));
        }
        os << '\n';
        json["plane"].push_back(plane_json);
    }

    os << std::string(50, '#') << '\n';
    os << boost::format("%15s %15s %15s\n") % "Frame" % "CoplaneIndex" % "STD";

    const auto fmt = boost::format("%15.3f %15.3f %15.3f\n");
    for (const auto &element : coplaneIndex | boost::adaptors::indexed(1)) {
        os << boost::format(fmt) % element.index() % element.value().first % element.value().second;
        json["Frame"].push_back(element.index());
        json["CoplaneIndex"].push_back(element.value().first);
        json["STD"].push_back(element.value().second);
    }

    os << std::string(50, '#') << '\n';

    os << ">>>JSON<<<\n";
    os << json;
    os << "<<<JSON>>>\n";
}

void CoplaneIndex::readInfo() {
    std::size_t num = choose(2, "Enter the number of plane > ");
    mask_arrays.resize(num);
    atom_arrays.resize(num);

    for (std::size_t i = 0; i < num; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            select1group(mask_arrays[i][j], "Enter mark for atom (" + std::to_string(j + 1) + ") of plane (" +
                                                std::to_string(i + 1) + ") > ");
        }
    }
}
