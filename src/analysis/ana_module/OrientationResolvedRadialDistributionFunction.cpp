//
// Created by xiamr on 9/4/19.
//

#include "OrientationResolvedRadialDistributionFunction.hpp"

#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "utils/VectorSelectorFactory.hpp"
#include "utils/common.hpp"

OrientationResolvedRadialDistributionFunction::OrientationResolvedRadialDistributionFunction() {
    enable_outfile = true;
}

void OrientationResolvedRadialDistributionFunction::processFirstFrame(std::shared_ptr<Frame> &frame) {
    water_OW_atoms.reserve(frame->atom_list.size() / 3);
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, water_Ow_atom_mask))
            water_OW_atoms.push_back(atom);
        else if (is_match(atom, reference_atom_mask)) {
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
    os << vectorSelector;
    os << "# distance_width(Ang) > " << distance_width << '\n';
    os << "# angle_width(deg) > " << angle_width * 180 / M_PI << '\n';
    os << "# max_distance(Ang) > " << max_distance << '\n';
    os << "# Temperature(K) > " << temperature << '\n';
    os << std::string(50, '#') << '\n';

    const double factor = -kb * temperature * avogadro_constant / 4184.0;
    double max_value = *std::max_element(hist.origin(), hist.origin() + hist.num_elements());
    os << format("#%15s %15s %15s %15s\n", "Distance(Ang)", "Angle(degree)", "g(r,theta)", "Energy(kcal/mol)");
    for (int i = 0; i < distance_bins; ++i) {
        for (int j = 0; j < angle_bins; ++j) {
            double pop = hist[j][i] / max_value;
            os << boost::format("%15.6f %15.6f %15.6f %15.6f\n") % ((i + 0.5) * distance_width) %
                      ((j + 0.5) * 180.0 * angle_width / M_PI) % hist[j][i] % (pop == 0.0 ? 100.0 : factor * log(pop));
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
    select1group(reference_atom_mask, "Eneter mask for reference atom > ");
    select1group(water_Ow_atom_mask, "Enter mask for OW atoms > ");

    vectorSelector = VectorSelectorFactory::getVectorSelector();
    vectorSelector->readInfo();

    distance_width = choose<double>(0, 1000, "Binwidth in r (Ang) [ 0.01 ] > ", Default(0.01));
    angle_width = choose<double>(0, 360, "Binwidth in theta (deg) [ 0.5 ] > ", Default(0.5));
    angle_width = angle_width / 180 * M_PI;

    max_distance = choose<double>(0, 1000, "Maximum value of r (Ang) [ 10.0 ] >", Default(10.0));

    distance_bins = static_cast<int>(max_distance / distance_width);
    angle_bins = M_PI / angle_width;

    hist.resize(boost::extents[angle_bins][distance_bins]);
    std::fill_n(hist.data(), hist.num_elements(), 0);

    temperature = choose(0.0, 10000.0, "Temperature [298] (K):", Default(298.0));
}
