//
// Created by xiamr on 6/14/19.
//
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include "Distance.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include "nlohmann/json.hpp"

Distance::Distance() {
    enable_outfile = true;
    enable_forcefield = true;
}

template<typename SinglePassRange>
std::tuple<double, double, double> Distance::calculate_mass_center(const SinglePassRange &atoms_group) {
    std::tuple<double, double, double> coord{};
    double weigh{};

    for (auto &atom : atoms_group) {
        double mass = atom->mass.value();
        coord += atom->getCoordinate() * mass;
        weigh += mass;
    }
    return coord / weigh;
}


void Distance::process(std::shared_ptr<Frame> &frame) {

    auto coord1 = calculate_mass_center(atoms_for_group1);
    auto coord2 = calculate_mass_center(atoms_for_group2);

    auto r = coord2 - coord1;
    frame->image(r);

    auto dist = vector_norm(r);
    acc(dist);

    distances.push_back(dist);
}


void Distance::print(std::ostream &os) {

    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# mask for group1 > " << mask_for_group1 << '\n';
    os << "# mask for group2 > " << mask_for_group2 << '\n';
    os << std::string(50, '#') << '\n';
    os << "# mean : " << boost::accumulators::mean(acc) << '\n';
    os << "# standard deviation : " << std::sqrt(boost::accumulators::variance(acc)) << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("#%15s %15s\n") % "Frame" % "Distance(Ang)";
    for (const auto &element : distances | boost::adaptors::indexed(1)) {
        os << boost::format(" %15d %15.8f\n") % element.index() % element.value();
    }
    os << std::string(50, '#') << '\n';

    os << ">>>JSON<<<\n";
    saveJson(os);
    os << "<<<JSON>>>\n";
}

void Distance::readInfo() {
    Atom::select2group(mask_for_group1, mask_for_group2);
}

void Distance::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](std::shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, mask_for_group1)) atoms_for_group1.insert(atom);
                        if (Atom::is_match(atom, mask_for_group2)) atoms_for_group2.insert(atom);
                    });
}

void Distance::saveJson(std::ostream &os) const {

    nlohmann::json json;

    json["title"] = title();
    json["group1"] = to_string(mask_for_group1);
    json["group2"] = to_string(mask_for_group2);
    json["mean"] = boost::accumulators::mean(acc);
    json["std"] = std::sqrt(boost::accumulators::variance(acc));
    json["unit"] = "Ang";
    json["distances"] = distances;

    os << json;
}

