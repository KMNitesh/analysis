//
// Created by xiamr on 6/14/19.
//

#include "DemixIndexOfTwoGroup.hpp"

#include <boost/range/adaptors.hpp>

#include "data_structure/frame.hpp"
#include "nlohmann/json.hpp"

DemixIndexOfTwoGroup::DemixIndexOfTwoGroup() { enable_outfile = true; }

auto DemixIndexOfTwoGroup::calculate_grid_index(const std::shared_ptr<Atom> &atom,
                                                const std::shared_ptr<Frame> &frame) {
    auto &[a_axis, b_axis, c_axis] = frame->box.getAxis();
    auto box_index_x = int(atom->x / (a_axis / grid_x)) % grid_x;
    while (box_index_x < 0) box_index_x += grid_x;
    auto box_index_y = int(atom->y / (b_axis / grid_y)) % grid_y;
    while (box_index_y < 0) box_index_y += grid_y;
    auto box_index_z = int(atom->z / (c_axis / grid_z)) % grid_z;
    while (box_index_z < 0) box_index_z += grid_z;

    if (!atom->mass) {
        std::cerr << "ERROR !!  Atom mass not available !\n";
        exit(EXIT_FAILURE);
    }

    assert(box_index_x >= 0 && box_index_x < grid_x);
    assert(box_index_y >= 0 && box_index_y < grid_y);
    assert(box_index_z >= 0 && box_index_z < grid_z);
    return std::make_tuple(box_index_x, box_index_y, box_index_z);
}

void DemixIndexOfTwoGroup::process(std::shared_ptr<Frame> &frame) {
    double group1_dens[grid_x][grid_y][grid_z];
    double group2_dens[grid_x][grid_y][grid_z];

    bzero(group1_dens, grid_x * grid_y * grid_z * sizeof(double));
    bzero(group2_dens, grid_x * grid_y * grid_z * sizeof(double));

    for (auto &atom : group1) {
        auto [box_index_x, box_index_y, box_index_z] = calculate_grid_index(atom, frame);
        group1_dens[box_index_x][box_index_y][box_index_z] += atom->mass.get();
    }
    for (auto &atom : group2) {
        auto [box_index_x, box_index_y, box_index_z] = calculate_grid_index(atom, frame);
        group2_dens[box_index_x][box_index_y][box_index_z] += atom->mass.get();
    }

    double d_sum = 0.0;
    double d1_sum = 0.0;
    double d2_sum = 0.0;
    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                auto d1 = group1_dens[i][j][k];
                d1_sum += d1;
                auto d2 = group2_dens[i][j][k];
                d2_sum += d2;
                if (d1 == 0.0 or d2 == 0.0) continue;
                d_sum += 1 / ((1 / d1) + (1 / d2));
            }
        }
    }

    d_sum /= frame->volume();

    auto d_ideal = 1 / (1 / d1_sum + 1 / d2_sum);
    d_ideal /= frame->volume();

    demix_index_list.emplace_back(d_sum, d_ideal);
}

void DemixIndexOfTwoGroup::readInfo() {
    Atom::select2group(mask1, mask2, "Please select group1 > ", "Please select group2 > ");
    grid_x = choose<int>(0, 1000, "Grid in X dememsion  :  ");
    grid_y = choose<int>(0, 1000, "Grid in Y dememsion  :  ");
    grid_z = choose<int>(0, 1000, "Grid in Z dememsion  :  ");
}

void DemixIndexOfTwoGroup::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "#    Demix Rate (normalization)\n";
    os << "#    Group1  " << mask1 << '\n';
    os << "#    Group2  " << mask2 << '\n';
    os << "#    Grid    X = " << grid_x << "  Y = " << grid_y << "  Z = " << grid_z << '\n';
    os << std::string(50, '#') << '\n';

    os << "@   title \"Demix Rate\"\n";
    os << "@    xaxis  label \"Frame Number\"\n";
    os << "@    yaxis  label \"Demix Rate\"\n";
    os << "@TYPE xy\n";
    os << "@ legend on\n";
    os << "@ legend length 1\n";
    os << "@ s0 legend \"Demix Rate\"\n";

    os << "# Frame      Demix Rate \n";
    for (const auto &element : demix_index_list | boost::adaptors::indexed(1)) {
        os << element.index() << "        " << std::get<0>(element.value()) / std::get<1>(element.value()) << '\n';
    }
    os << std::string(50, '#') << '\n';

    nlohmann::json json;

    json["title"] = "Demix Rate (normalization)";
    json["group1"] = to_string(mask1);
    json["group2"] = to_string(mask2);

    json["grid"] = {{"X", grid_x}, {"Y", grid_y}, {"Z", grid_z}};

    json["DemixRate"] = demix_index_list;

    os << ">>>JSON<<<\n";
    os << json;
    os << "<<<JSON>>>\n";
}

void DemixIndexOfTwoGroup::processFirstFrame(std::shared_ptr<Frame> &frame) {
    if (frame->box.get_box_type() != PBCBox::Type::orthogonal) {
        std::cerr << "Not supported \n";
        std::exit(1);
    }
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(), [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, this->mask1)) this->group1.insert(atom);
        if (Atom::is_match(atom, this->mask2)) this->group2.insert(atom);
    });
}

void DemixIndexOfTwoGroup::setParameters(const AmberMask &id1, const AmberMask &id2, const Grid &grid,
                                         const std::string &outfilename) {
    this->mask1 = id1;
    this->mask2 = id2;

    this->grid_x = grid.x;
    this->grid_y = grid.y;
    this->grid_z = grid.z;

    if (grid.x < 1 or grid.y < 1 or grid.z < 1) {
        throw std::runtime_error("grid component must large than one");
    }

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }
}
