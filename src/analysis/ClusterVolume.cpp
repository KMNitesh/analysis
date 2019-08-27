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
                    });
}


void ClusterVolume::process(std::shared_ptr<Frame> &frame) {
    bool grid[grid_x][grid_y][grid_z];
    bzero(grid, grid_x * grid_y * grid_z * sizeof(bool));

    double grid_x_step = frame->a_axis / grid_x;
    double grid_y_step = frame->b_axis / grid_y;
    double grid_z_step = frame->c_axis / grid_z;

    for (auto &atom : atom_group) {

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

                        grid[box_index_x][box_index_y][box_index_z] = true;
                    }
                }
            }
        }
    }

    std::size_t num_grid_point = 0;
    for (int i = 0; i < grid_x; i++) {
        for (int j = 0; j < grid_y; j++) {
            for (int k = 0; k < grid_z; k++) {
                if (grid[i][j][k]) {
                    num_grid_point++;
                }
            }
        }
    }


    auto total_volume = frame->a_axis * frame->b_axis * frame->c_axis; // Ang^3

    double total_grid_points = grid_x * grid_y * grid_z;

    auto percentage = num_grid_point / total_grid_points;

    auto volume = total_volume * percentage;

    volumes.emplace_back(percentage, volume);

}

void ClusterVolume::print(std::ostream &os) {
    os << "#####################################\n";
    os << "#    " << title() << '\n';
    os << "#    Group  " << atom_mask << '\n';
    os << "#    Grid    X = " << grid_x << "  Y = " << grid_y << "  Z = " << grid_z << '\n';
    os << "#####################################\n";

    os << boost::format("%15s %15s %15s\n") % "Frame" % "Pencentage(%)" % "Volume(Ang^3)";
    for (const auto &element: volumes | boost::adaptors::indexed(1)) {
        os << boost::format("%15.3f %15.5f %15.5f\n")
              % element.index() % (100 * element.value().first) % element.value().second;
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
