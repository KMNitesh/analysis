//
// Created by xiamr on 9/4/19.
//

#include <boost/range/algorithm.hpp>
#include "OrientationResolvedRadialDistributionFunction.hpp"
#include "common.hpp"
#include "frame.hpp"
#include "VectorSelectorFactory.hpp"

OrientationResolvedRadialDistributionFunction::OrientationResolvedRadialDistributionFunction() {
    enable_outfile = true;
}

void OrientationResolvedRadialDistributionFunction::processFirstFrame(std::shared_ptr<Frame> &frame) {
    water_OW_atoms.reserve(frame->atom_list.size() / 3);
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, water_Ow_atom_mask)) water_OW_atoms.push_back(atom);
        else if (Atom::is_match(atom, reference_atom_mask)) {
            if (reference_atom) {
                std::cerr << "ERROR !! Only one reference atom allowed for module <" << title() << ">\n";
                exit(1);
            } else {
                reference_atom = atom;
            }
        }
    });
    water_OW_atoms.shrink_to_fit();
}

void OrientationResolvedRadialDistributionFunction::process(std::shared_ptr<Frame> &frame) {
    nframes++;
    auto volume = frame->volume();
    for (auto &ow : water_OW_atoms) {
        auto v1 = vectorSelector->calculateVector(ow->molecule.lock(), frame);
        auto v2 = ow->getCoordinate() - reference_atom->getCoordinate();
        frame->image(v2);

        auto d_ref = vector_norm(v2);

        auto theta = acos(dot_multiplication(v1, v2) / (d_ref * vector_norm(v1)));

        // determine bin indices
        const int bin_r = int(d_ref / distance_width);
        const int bin_th = int(theta / angle_width);
        if ((0 <= bin_r) && (bin_r < distance_bins) && (bin_th >= 0.0) && (bin_th < angle_bins)) {
            // add increment
            hist[bin_th][bin_r] += volume;
        }
    }
}

void OrientationResolvedRadialDistributionFunction::print(std::ostream &os) {
    normalize();
    os << std::string(50, '#') << '\n';
    os << title() << '\n';
    os << "# reference atom > " << reference_atom_mask << '\n';
    os << "# water Ow atoms > " << water_Ow_atom_mask << '\n';
    os << "# distance_width(Ang) > " << distance_width << '\n';
    os << "# angle_width(deg) > " << angle_width * 180 / M_PI << '\n';
    os << "# max_distance(Ang) > " << max_distance << '\n';
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Distance(Ang)", "Angle(degree)", "g(r,theta)");
    write(os);
}

void OrientationResolvedRadialDistributionFunction::write(std::ostream &os) const {
    // write r bin centers
    for (int i = 0; i < this->distance_bins; ++i) {
        for (int j = 0; j < this->angle_bins; ++j) {
            os << boost::format("%15.6f %15.6f %15.6f\n")
                  % ((i + 0.5) * this->distance_width)
                  % ((j + 0.5) * 180.0 * this->angle_width / M_PI)
                  % hist[j][i];
        }
    }
}

void OrientationResolvedRadialDistributionFunction::normalize() {
    auto n_factor = static_cast<double>(nframes) * water_OW_atoms.size();
    std::vector<double> normr, normth;

    normr.resize(distance_bins, 0);
    normth.resize(angle_bins, 0);

    // prepare volume element factors in r

    auto r = 0.5 * distance_width;
    for (int j = 0; j < distance_bins; ++j) {
        normr[j] = (4.0 / 3.0) * M_PI * distance_width * (3 * r * r + 0.25 * distance_width * distance_width);
        r += distance_width;
    }

    // prepare volume element factors in theta
    auto th = 0.5 * angle_width;
    for (int i = 0; i < angle_bins; ++i) {
        normth[i] = angle_width * 0.5 * sin(th);
        th += angle_width;
    }

    // apply normalization
    for (int i = 0; i < angle_bins; ++i) {
        for (int j = 0; j < distance_bins; ++j) {
            hist[i][j] /= normth[i] * normr[j] * n_factor;
        }
    }
}

void OrientationResolvedRadialDistributionFunction::readInfo() {

    Atom::select1group(reference_atom_mask, "Eneter mask for reference atom > ");
    Atom::select1group(water_Ow_atom_mask, "Enter mask for OW atoms > ");

    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();

    distance_width = choose<double>(0, 1000, "Binwidth in r (Ang) [ 0.01 ] > ", true, 0.01);
    angle_width = choose<double>(0, 360, "Binwidth in theta (deg) [ 0.5 ] > ", true, 0.5);
    angle_width = angle_width / 180 * M_PI;

    max_distance = choose<double>(0, 1000, "Maximum value of r (Ang) [ 10.0 ] >", true, 10.0);

    distance_bins = static_cast<int>(max_distance / distance_width);
    angle_bins = M_PI / angle_width;

    hist.resize(boost::extents[angle_bins][distance_bins]);
    std::fill_n(hist.data(), hist.num_elements(), 0);
}
