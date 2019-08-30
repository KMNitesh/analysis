//
// Created by xiamr on 8/26/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>
#include <tbb/parallel_while.h>
#include "ClusterVolume.hpp"
#include "common.hpp"
#include "frame.hpp"

ClusterVolume::ClusterVolume() {
    enable_outfile = true;
    enable_tbb = true;
}

void ClusterVolume::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list,
                    [this](std::shared_ptr<Atom> &atom) {
                        if (Atom::is_match(atom, this->atom_mask)) this->atom_group.insert(atom);
                        else this->other_atoms.insert(atom);
                    });

    for (auto &atom : atom_group) {
        radii_for_atom_group.push_back(getVdwRadii(atom));
    }
    for (auto &atom : other_atoms) {
        radii_for_other_atoms.push_back(getVdwRadii(atom));
    }
}

void
ClusterVolume::fill_atom(boost::multi_array<ATOM_Category, 3> *grid,
                         double grid_x_step, double grid_y_step, double grid_z_step,
                         ATOM_Category category,
                         double x, double y, double z, double radii) const {

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

                    (*grid)[box_index_x][box_index_y][box_index_z] = category;
                }
            }
        }
    }
}

void ClusterVolume::process(std::shared_ptr<Frame> &frame) {


    auto[grid, atom_group_array, other_atom_array,
    grid_x_step, grid_y_step, grid_z_step,
    total_volume, nframe] = preprocess(frame);

    work_body(grid, atom_group_array, other_atom_array,
              grid_x_step, grid_y_step, grid_z_step,
              total_volume, nframe);
}

void ClusterVolume::work_body(boost::multi_array<ATOM_Category, 3> *grid,
                              boost::multi_array<double, 2> *atom_group_array,
                              boost::multi_array<double, 2> *other_atom_array,
                              double grid_x_step,
                              double grid_y_step,
                              double grid_z_step,
                              double total_volume,
                              int nframe) {

    double total_grid_points = grid_x * grid_y * grid_z;
    auto[num_grid_point_before_fill, num_grid_point_after_fill] = do_grid(grid, atom_group_array, other_atom_array,
                                                                          grid_x_step, grid_y_step, grid_z_step);

    auto percentage_before_fill = num_grid_point_before_fill / total_grid_points;
    auto volume_before_fill = total_volume * percentage_before_fill;
    auto percentage_after_fill = num_grid_point_after_fill / total_grid_points;
    auto volume_after_fill = total_volume * percentage_after_fill;

    std::decay_t<decltype(volumes)>::accessor accessor;
    volumes.insert(accessor, nframe);
    accessor->second = {percentage_before_fill, volume_before_fill, percentage_after_fill, volume_after_fill};
    accessor.release();
}

std::pair<std::size_t, std::size_t> ClusterVolume::do_grid(boost::multi_array<ATOM_Category, 3> *grid,
                                                           boost::multi_array<double, 2> *atom_group_array,
                                                           boost::multi_array<double, 2> *other_atom_array,
                                                           double grid_x_step,
                                                           double grid_y_step,
                                                           double grid_z_step) const {


    for (auto index : boost::irange(radii_for_other_atoms.size())) {
        fill_atom(grid, grid_x_step, grid_y_step, grid_z_step, ATOM_Category::OTHER,
                  (*other_atom_array)[index][0], (*other_atom_array)[index][1], (*other_atom_array)[index][2],
                  radii_for_other_atoms[index]);
    }

    for (auto index : boost::irange(radii_for_atom_group.size())) {
        fill_atom(grid, grid_x_step, grid_y_step, grid_z_step, ATOM_Category::DEST,
                  (*atom_group_array)[index][0], (*atom_group_array)[index][1], (*atom_group_array)[index][2],
                  radii_for_atom_group[index]);
    }
    auto num_grid_point_before_fill = countFilledGridPoints(grid);
    while (fill_space(grid, grid_x_step, grid_y_step, grid_z_step));
    auto num_grid_point_after_fill = countFilledGridPoints(grid);
    return {num_grid_point_before_fill, num_grid_point_after_fill};
}


bool
ClusterVolume::fill_space(boost::multi_array<ATOM_Category, 3> *grid,
                          double grid_x_step, double grid_y_step, double grid_z_step) const {

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
                if ((*grid)[i][j][k] == ATOM_Category::EMPTY) {
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

                        if ((*grid)[neigh_x][neigh_y][neigh_z] == ATOM_Category::OTHER) {
                            no_water_surround = false;
                            break;
                        } else if ((*grid)[neigh_x][neigh_y][neigh_z] == ATOM_Category::DEST) {
                            found_dest_atom = true;
                        }
                    }
                    if (no_water_surround and found_dest_atom) {
                        (*grid)[i][j][k] = ATOM_Category::DEST;
                        updated = true;
                    }
                }
            }
        }
    }
    return updated;
}

size_t ClusterVolume::countFilledGridPoints(boost::multi_array<ATOM_Category, 3> *grid) const {
    std::size_t num_grid_point = 0;
    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                if ((*grid)[i][j][k] == ATOM_Category::DEST) {
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
    for (int index : boost::irange(1, current_frame_num + 1)) {
        std::decay_t<decltype(volumes)>::const_accessor accessor;
        volumes.find(accessor, index);
        os << boost::format("%10d %20.5f %20.5f %20.5f %20.5f\n")
              % accessor->first
              % (100 * std::get<0>(accessor->second))
              % std::get<1>(accessor->second)
              % (100 * std::get<2>(accessor->second))
              % std::get<3>(accessor->second);
        accessor.release();
    }
}

void ClusterVolume::readInfo() {
    Atom::select1group(atom_mask, "Please enter AmberMask for selected group > ");
    grid_x = choose<int>(1, 10000, "Grid in X dememsion  :  ");
    grid_y = choose<int>(1, 10000, "Grid in Y dememsion  :  ");
    grid_z = choose<int>(1, 10000, "Grid in Z dememsion  :  ");


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

bool ClusterVolume::enable_paralel_while_impl() {
    return true;
}


ClusterVolume::argument_type ClusterVolume::preprocess(std::shared_ptr<Frame> &frame) {
    current_frame_num++;

    auto grid = new boost::multi_array<ATOM_Category, 3>(boost::extents[grid_x][grid_y][grid_z]);
    std::fill_n(grid->data(), grid->num_elements(), ATOM_Category::EMPTY);

    auto atom_group_array = new boost::multi_array<double, 2>(boost::extents[atom_group.size()][3]);
    auto other_atom_array = new boost::multi_array<double, 2>(boost::extents[other_atoms.size()][3]);

    int i = 0;
    for (auto &atom : atom_group) {
        (*atom_group_array)[i][0] = atom->x;
        (*atom_group_array)[i][1] = atom->y;
        (*atom_group_array)[i][2] = atom->z;
        ++i;
    }
    i = 0;
    for (auto &atom : other_atoms) {
        (*other_atom_array)[i][0] = atom->x;
        (*other_atom_array)[i][1] = atom->y;
        (*other_atom_array)[i][2] = atom->z;
        ++i;
    }

    double grid_x_step = frame->a_axis / grid_x;
    double grid_y_step = frame->b_axis / grid_y;
    double grid_z_step = frame->c_axis / grid_z;

    auto total_volume = frame->a_axis * frame->b_axis * frame->c_axis; // Ang^3
    return {grid, atom_group_array, other_atom_array, grid_x_step, grid_y_step, grid_z_step, total_volume,
            current_frame_num};
}

void ClusterVolume::do_parallel_while_impl(std::function<std::shared_ptr<Frame>()> func) {

    FrameStream stream(func, this);

    ApplyBody body(this);

    tbb::parallel_while<ApplyBody> w;
    w.run(stream, body);
}

bool ClusterVolume::FrameStream::pop_if_present(ClusterVolume::argument_type &item) {
    auto frame = func();
    if (frame) {
        if (parent->current_frame_num == 0) {
            parent->processFirstFrame(frame);
        }
        item = parent->preprocess(frame);
        return true;
    } else {
        return false;
    }
}

void ClusterVolume::ApplyBody::operator()(ClusterVolume::argument_type &item) const {
    auto&[grid, atom_group_array, other_atom_array,
    grid_x_step, grid_y_step, grid_z_step,
    total_volume, nframe] = item;

    parent->work_body(grid, atom_group_array, other_atom_array,
                      grid_x_step, grid_y_step, grid_z_step,
                      total_volume, nframe);
    delete std::get<0>(item);
    delete std::get<1>(item);
    delete std::get<2>(item);
}

