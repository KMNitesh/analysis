//
// Created by xiamr on 6/14/19.
//
#include <boost/range/adaptors.hpp>
#include "RMSFCal.hpp"
#include "frame.hpp"

void RMSFCal::process(std::shared_ptr<Frame> &frame) {
    steps++;
    int nfit, n;
    nfit = n = static_cast<int>(this->group.size());
    BOOST_ASSERT_MSG(group.size() < ATOM_MAX, "need to increase ATOM_MAX");

    if (!first_frame) {
        first_frame = false;
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x1[index] = atom->x;
                first_y = y1[index] = atom->y;
                first_z = z1[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x1[index] = first_x + xr;
                y1[index] = first_y + yr;
                z1[index] = first_z + zr;
            }
            index++;
        }

        double mid[3];
        center(n, x1, y1, z1, mid, nfit);
        std::map<int, double> f1x, f1y, f1z;
        for (unsigned index = 0; index < this->group.size(); index++) {
            x_avg[index] = x1[index];
            y_avg[index] = y1[index];
            z_avg[index] = z1[index];
            f1x[index] = x1[index];
            f1y[index] = y1[index];
            f1z[index] = z1[index];
        }
        x[steps] = f1x;
        y[steps] = f1y;
        z[steps] = f1z;
    } else {
        int index = 0;
        bool first_atom = true;
        double first_x, first_y, first_z;
        for (auto &atom : this->group) {
            if (first_atom) {
                first_atom = false;
                first_x = x2[index] = atom->x;
                first_y = y2[index] = atom->y;
                first_z = z2[index] = atom->z;
            } else {
                double xr = atom->x - first_x;
                double yr = atom->y - first_y;
                double zr = atom->z - first_z;
                frame->image(xr, yr, zr);
                x2[index] = first_x + xr;
                y2[index] = first_y + yr;
                z2[index] = first_z + zr;
            }
            index++;
        }
        double mid[3];
        center(n, x2, y2, z2, mid, nfit);
        quatfit(n, x1, y1, z1, n, x2, y2, z2, nfit);
        std::map<int, double> f2x, f2y, f2z;
        for (unsigned index = 0; index < this->group.size(); index++) {
            x_avg[index] += x2[index];
            y_avg[index] += y2[index];
            z_avg[index] += z2[index];
            f2x[index] = x2[index];
            f2y[index] = y2[index];
            f2z[index] = z2[index];
        }
        x[steps] = f2x;
        y[steps] = f2y;
        z[steps] = f2z;
    }

}

void RMSFCal::print(std::ostream &os) {

    for (unsigned int index = 0; index < this->group.size(); index++) {
        x_avg[index] /= this->steps;
        y_avg[index] /= this->steps;
        z_avg[index] /= this->steps;
    }

    os << "***************************\n";
    os << "****** RMSF Calculator ****\n";
    os << "SET:" << ids << '\n';
    os << "***************************\n";
    for (const auto &element : group | boost::adaptors::indexed(1)) {
        os << element.value()->seq << "     " << rmsvalue(element.index()) << '\n';
    }
    os << "***************************\n";
}

void RMSFCal::readInfo() {
    Atom::select1group(ids, "Please enter group:");
}

double RMSFCal::rmsvalue(int index) {

    double dx2_y2_z2 = 0.0;
    double dx, dy, dz;
    for (int frame = 1; frame <= steps; frame++) {
        dx = x[frame][index] - x_avg[index];
        dy = y[frame][index] - y_avg[index];
        dz = z[frame][index] - z_avg[index];
        dx2_y2_z2 += dx * dx + dy * dy + dz * dz;
    }

    return sqrt(dx2_y2_z2 / steps);

}

void RMSFCal::find_matched_atoms(std::shared_ptr<Frame> &frame) {
    if (first_frame) {
        std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                      [this](std::shared_ptr<Atom> &atom) {
                          if (Atom::is_match(atom, ids)) group.insert(atom);
                      });
    }
}
