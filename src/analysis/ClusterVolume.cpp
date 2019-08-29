//
// Created by xiamr on 8/26/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include "ClusterVolume.hpp"
#include "common.hpp"
#include "frame.hpp"

ClusterVolume::ClusterVolume() {
    enable_outfile = true;
}

void ClusterVolume::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](std::shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, this->atom_mask)) this->atom_group.insert(atom);
                        else this->other_atoms.insert(atom);
                    });
}

void ClusterVolume::fill_atom(double grid_x_step, double grid_y_step, double grid_z_step,
                              const std::shared_ptr<Atom> &atom, ATOM_Category category) {

    double x = atom->x;
    double y = atom->y;
    double z = atom->z;

    double radii = getVdwRadii(atom);
    auto radii2 = radii * radii;

    int radii_x = std::ceil(radii / grid_x_step);
    int radii_y = std::ceil(radii / grid_y_step);
    int radii_z = std::ceil(radii / grid_z_step);

    for (int i = -radii_x; i <= radii_x; ++i) {
        for (int j = -radii_y; j <= radii_y; ++j) {
            for (int k = -radii_z; k <= radii_z; ++k) {

                auto xr = i * grid_x_step;
                auto yr = j * grid_y_step;
                auto zr = k * grid_z_step;

                if ((xr * xr + yr * yr + zr * zr) < radii2) {

                    auto box_index_x = int((x + xr) / grid_x_step) % grid_x;
                    while (box_index_x < 0) box_index_x += grid_x;
                    auto box_index_y = int((y + yr) / grid_y_step) % grid_y;
                    while (box_index_y < 0) box_index_y += grid_y;
                    auto box_index_z = int((z + zr) / grid_z_step) % grid_z;
                    while (box_index_z < 0) box_index_z += grid_z;

                    assert(box_index_x >= 0 && box_index_x < grid_x);
                    assert(box_index_y >= 0 && box_index_y < grid_y);
                    assert(box_index_z >= 0 && box_index_z < grid_z);

                    grid[box_index_x][box_index_y][box_index_z] = category;
                }
            }
        }
    }
}

void ClusterVolume::process(std::shared_ptr<Frame> &frame) {

    std::fill_n(grid.data(), grid.num_elements(), ATOM_Category::EMPTY);

    double grid_x_step = frame->a_axis / grid_x;
    double grid_y_step = frame->b_axis / grid_y;
    double grid_z_step = frame->c_axis / grid_z;

    for (auto &atom : other_atoms) {
        fill_atom(grid_x_step, grid_y_step, grid_z_step, atom, ATOM_Category::OTHER);
    }

    for (auto &atom : atom_group) {
        fill_atom(grid_x_step, grid_y_step, grid_z_step, atom, ATOM_Category::DEST);
    }

    std::size_t num_grid_point_before_fill = countFilledGridPoints();

    while (fill_space(grid_x_step, grid_y_step, grid_z_step));

    std::size_t num_grid_point_after_fill = countFilledGridPoints();


    auto total_volume = frame->a_axis * frame->b_axis * frame->c_axis; // Ang^3

    double total_grid_points = grid_x * grid_y * grid_z;

    auto percentage_before_fill = num_grid_point_before_fill / total_grid_points;
    auto volume_before_fill = total_volume * percentage_before_fill;

    auto percentage_after_fill = num_grid_point_after_fill / total_grid_points;
    auto volume_after_fill = total_volume * percentage_after_fill;

    volumes.emplace_back(percentage_before_fill, volume_before_fill, percentage_after_fill, volume_after_fill);

}


bool ClusterVolume::fill_space(double grid_x_step, double grid_y_step, double grid_z_step) {

    bool updated = false;

    auto neighs = generate_neighbor_grids(grid_x_step, grid_y_step, grid_z_step);

    auto total_size = neighs.size() * 3;
    int neighs_vec[total_size];


    int *p = neighs_vec;
    for (auto[neigh_x, neigh_y, neigh_z] : neighs) {
        *p = neigh_x;
        ++p;
        *p = neigh_y;
        ++p;
        *p = neigh_z;
        ++p;
    }


    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                if (grid[i][j][k] == ATOM_Category::EMPTY) {
                    bool no_water_surround = true;
                    bool found_dest_atom = false;
                    p = neighs_vec;
                    while (p < (neighs_vec + total_size)) {

                        auto neigh_x = *p;
                        ++p;
                        auto neigh_y = *p;
                        ++p;
                        auto neigh_z = *p;
                        ++p;

                        neigh_x = (neigh_x + i) % grid_x;
                        neigh_x = neigh_x < 0 ? neigh_x + grid_x : neigh_x;

                        neigh_y = (neigh_y + j) % grid_y;
                        neigh_y = neigh_y < 0 ? neigh_y + grid_y : neigh_y;

                        neigh_z = (neigh_z + k) % grid_z;
                        neigh_z = neigh_z < 0 ? neigh_z + grid_z : neigh_z;

                        if (grid[neigh_x][neigh_y][neigh_z] == ATOM_Category::OTHER) {
                            no_water_surround = false;
                            break;
                        } else if (grid[neigh_x][neigh_y][neigh_z] == ATOM_Category::DEST) {
                            found_dest_atom = true;
                        }
                    }
                    if (no_water_surround and found_dest_atom) {
                        grid[i][j][k] = ATOM_Category::DEST;
                        updated = true;
                    }
                }
            }
        }
    }
    return updated;
}

size_t ClusterVolume::countFilledGridPoints() const {
    std::size_t num_grid_point = 0;
    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                if (grid[i][j][k] == ATOM_Category::DEST) {
                    num_grid_point++;
                }
            }
        }
    }
    return num_grid_point;
}

void ClusterVolume::print(std::ostream &os) {
    os << "#####################################\n";
    os << "#    " << title() << '\n';
    os << "#    Group  " << atom_mask << '\n';
    os << "#    Grid    X = " << grid_x << "  Y = " << grid_y << "  Z = " << grid_z << '\n';
    os << "#####################################\n";

    os << boost::format("%10s %20s %20s %20s %20s\n")
          % "Frame"
          % "vdW Volume Pencentage(%)"
          % "vdW Volume(Ang^3)"
          % "Fill Space Volume Pencentage(%)"
          % "Fill Space Volume(Ang^3)";
    for (const auto &element: volumes | boost::adaptors::indexed(1)) {
        os << boost::format("%10d %20.5f %20.5f %20.5f %20.5f\n")
              % element.index()
              % (100 * std::get<0>(element.value()))
              % std::get<1>(element.value())
              % (100 * std::get<2>(element.value()))
              % std::get<3>(element.value());
    }
}

void ClusterVolume::readInfo() {
    Atom::select1group(atom_mask, "Please enter AmberMask for selected group > ");
    grid_x = choose<int>(1, 10000, "Grid in X dememsion  :  ");
    grid_y = choose<int>(1, 10000, "Grid in Y dememsion  :  ");
    grid_z = choose<int>(1, 10000, "Grid in Z dememsion  :  ");

    grid.resize(boost::extents[grid_x][grid_y][grid_z]);
}

double ClusterVolume::getVdwRadii(const std::shared_ptr<Atom> &atom) {
    auto sym = which(atom);
    auto it = vdWRadiis.find(sym);
    if (it != vdWRadiis.end()) {
        return it->second;
    } else {
        std::cerr << "vdW radii not available for atom "
                     "( name = " << atom->atom_name << ", mass = " << atom->mass << ")\n";
        exit(1);
    }
}

std::vector<std::tuple<int, int, int>>
ClusterVolume::generate_neighbor_grids(double grid_x_step, double grid_y_step, double grid_z_step) const {

    std::vector<std::tuple<int, int, int>> neighs;
    const double radii = 3.0;
    auto radii2 = radii * radii;

    int radii_x = std::ceil(radii / grid_x_step);
    int radii_y = std::ceil(radii / grid_y_step);
    int radii_z = std::ceil(radii / grid_z_step);

    for (int i = -radii_x; i <= radii_x; ++i) {
        for (int j = -radii_y; j <= radii_y; ++j) {
            for (int k = -radii_z; k <= radii_z; ++k) {

                auto xr = i * grid_x_step;
                auto yr = j * grid_y_step;
                auto zr = k * grid_z_step;

                if ((xr * xr + yr * yr + zr * zr) < radii2) {

                    auto box_index_x = int(xr / grid_x_step) % grid_x;
                    auto box_index_y = int(yr / grid_y_step) % grid_y;
                    auto box_index_z = int(zr / grid_z_step) % grid_z;

                    neighs.emplace_back(box_index_x, box_index_y, box_index_z);
                }
            }
        }
    }

    return neighs;
}
