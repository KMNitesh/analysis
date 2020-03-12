
#include "AngleDistributionBetweenTwoVectorWithCutoff.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "data_structure/frame.hpp"
#include "nlohmann/json.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/VectorSelectorFactory.hpp"
#include "utils/common.hpp"

AngleDistributionBetweenTwoVectorWithCutoff::AngleDistributionBetweenTwoVectorWithCutoff() { enable_outfile = true; }

void AngleDistributionBetweenTwoVectorWithCutoff::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, this->metal_mask))
            this->group1.insert(atom);
        if (Atom::is_match(atom, this->ligand_mask))
            this->group2.insert(atom);
    });
    throw_assert(this->group1.size() == 1, "Group1 must has only one atom");
    throw_assert(!this->group1.empty(), "Group2 must has at least one atom");
}

void AngleDistributionBetweenTwoVectorWithCutoff::process(std::shared_ptr<Frame> &frame) {
    nframe++;
    for (auto &ref : group1) {
        for (auto &atom : group2) {
            double distance = atom_distance(ref, atom, frame);
            if (cutoff1 <= distance and distance < cutoff2) {
                auto v1 = vector1->calculateVector(ref->molecule.lock(), frame);
                auto v2 = vector2->calculateVector(atom->molecule.lock(), frame);

                double cos_ = dot_multiplication(v1, v2) / sqrt(vector_norm2(v1) * vector_norm2(v2));

                cos_hist.update(cos_);

                double angle = radian * acos(cos_);

                angle_evolution.emplace(nframe, std::make_pair<int>(atom->seq, angle));

                angle_hist.update(angle);
            }
        }
    }
}

void AngleDistributionBetweenTwoVectorWithCutoff::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# Group1 > " << metal_mask << '\n';
    os << "# Group2 > " << ligand_mask << '\n';

    os << "# vector1 > ";
    os << vector1;
    os << "# vector2 > ";
    os << vector2;

    os << "# angle max   > " << angle_hist.dimension_range.second << '\n';
    os << "# angle width > " << angle_hist.getWidth() << '\n';

    os << "# cutoff1 > " << cutoff1 << '\n';
    os << "# cutoff2 > " << cutoff2 << '\n';

    std::set<int> appeared_atoms;

    for (auto &it : angle_evolution) {
        appeared_atoms.insert(it.second.first);
    }

    os << std::string(25, '#') << "Angle(degree)" << std::string(25, '#') << '\n';
    os << format("#%15s", "Frame");
    for (auto i : appeared_atoms) {
        os << format(" [ %10d ]", i);
    }
    os << '\n';

    for (auto i = 1; i <= nframe; ++i) {
        os << format("%15.3d", i);
        std::map<int, double> mapping;
        for (auto it = angle_evolution.lower_bound(i); it != angle_evolution.upper_bound(i); ++it) {
            mapping.insert(it->second);
        }
        for (auto seq : appeared_atoms) {
            auto it = mapping.find(seq);
            if (it != mapping.end()) {
                os << format(" %15.3f", it->second);
            } else {
                os << std::string(16, ' ');
            }
        }
        os << '\n';
    }

    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Angle(degree)", "Probability Density(% degree-1)");
    for (auto [grid, value] : angle_hist.getDistribution()) {
        os << format("%15.3f %15.3f\n", grid, 100 * value);
    }
    os << std::string(50, '#') << '\n';

    os << format("#%15s %15s\n", "cos(theta)", "Probability Density(%)");
    for (auto [grid, value] : cos_hist.getDistribution()) {
        os << format("%15.3f %15.3f\n", grid, 100 * value);
    }
    os << std::string(50, '#') << '\n';

    os << ">>>JSON<<<\n";
    saveJson(os);
    os << "<<<JSON>>>\n";
}

void AngleDistributionBetweenTwoVectorWithCutoff::readInfo() {
    Atom::select2group(metal_mask, ligand_mask);

    std::cout << "For first Vector\n";
    vector1 = VectorSelectorFactory::getVectorSelector();
    vector1->readInfo();

    std::cout << "For second Vector\n";
    vector2 = VectorSelectorFactory::getVectorSelector();
    vector2->readInfo();

    double angle_max = choose(0.0, 180.0, "Enter Maximum Angle to Accumulate[180.0 degree]:", Default(180.0));
    double angle_width = choose(0.0, 180.0, "Enter Width of Angle Bins [0.5 degree]:", Default(0.5));

    cutoff1 = choose(0.0, 100.0, "Cutoff1 [Angstrom]:");
    cutoff2 = choose(0.0, 100.0, "Cutoff2 [Angstrom]:");

    throw_assert(cutoff1 < cutoff2, "Cutoff1 must less than Cutoff2");

    angle_hist.initialize(angle_max, angle_width);

    init_cos_hist(angle_max, angle_width);
}

std::string AngleDistributionBetweenTwoVectorWithCutoff::description() {
    std::stringstream ss;
    std::string title_line = "------ " + std::string(title()) + " ------";
    ss << title_line << "\n";
    ss << " M                 = [ " << metal_mask << " ]\n";
    ss << " L                 = [ " << ligand_mask << " ]\n";
    ss << " vector1           = " << vector1->description() << "\n";
    ss << " vector2           = " << vector2->description() << "\n";
    ss << " angle_max         = " << angle_hist.dimension_range.second << " (degree)\n";
    ss << " angle_width       = " << angle_hist.getWidth() << " (degree)\n";
    ss << " cutoff1           = " << cutoff1 << " (Ang)\n";
    ss << " cutoff2           = " << cutoff2 << " (Ang)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << std::string(title_line.size(), '-') << '\n';
    return ss.str();
}

void AngleDistributionBetweenTwoVectorWithCutoff::setParameters(const AmberMask &M, const AmberMask &L,
                                                                std::shared_ptr<VectorSelector> vector1,
                                                                std::shared_ptr<VectorSelector> vector2,
                                                                double angle_max, double angle_width, double cutoff1,
                                                                double cutoff2, const std::string &outfilename) {
    this->metal_mask = M;
    this->ligand_mask = L;

    if (!vector1) {
        throw std::runtime_error("vector1 not vaild");
    } else {
        this->vector1 = vector1;
    }

    if (!vector2) {
        throw std::runtime_error("vector2 not vaild");
    } else {
        this->vector2 = vector2;
    }

    if (angle_max < 0) {
        throw std::runtime_error("`angle_max` must not be negative");
    }
    if (angle_width < 0) {
        throw std::runtime_error("`cutoff1` must not be negative");
    }

    if (cutoff1 < 0) {
        throw std::runtime_error("`cutoff1` must not be negative");
    } else {
        this->cutoff1 = cutoff1;
    }

    if (cutoff2 <= 0) {
        throw std::runtime_error("`cutoff2` must be postive");
    } else {
        this->cutoff2 = cutoff2;
    }

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw std::runtime_error("outfilename cannot empty");
    }

    angle_hist.initialize(angle_max, angle_width);
    init_cos_hist(angle_max, angle_width);
}

void AngleDistributionBetweenTwoVectorWithCutoff::init_cos_hist(double angle_max, double angle_width) {
    cos_hist.initialize({cos(angle_max / radian), 1}, abs(1 - cos(angle_max / radian)) / (angle_max / angle_width));
}

void AngleDistributionBetweenTwoVectorWithCutoff::saveJson(std::ostream &os) const {
    nlohmann::json json;

    json["title"] = title();
    json["M"] = to_string(metal_mask);
    json["L"] = to_string(ligand_mask);
    json["vector1"] = vector1->description();
    json["vector2"] = vector2->description();
    json["angle_max"] = {{"value", angle_hist.dimension_range.second}, {"unit", "degree"}};

    json["angle_width"] = {{"value", angle_hist.getWidth()}, {"unit", "degree"}};

    json["cutoff1"] = {{"value", cutoff1}, {"unit", "Ang"}};

    json["cutoff2"] = {{"value", cutoff2}, {"unit", "Ang"}};

    json["Probability Density"] = {{"meta", {{"X", "Angle(degree)"}, {"Y", "Probability Density(degree-1)"}}},
                                   {"values", angle_hist.getDistribution()}};

    json["Probability Density by cosine"] = {{"meta", {{"X", "cos(theta)"}, {"Y", "Probability Density"}}},
                                             {"values", cos_hist.getDistribution()}};

    os << json;
}