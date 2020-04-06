#include "EnergyCalculator.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"

#include <boost/range/algorithm.hpp>

EnergyCalculator::EnergyTerm EnergyCalculator::calculate_energy(const std::shared_ptr<Frame> &frame) {

    const auto &[a_axis, b_axis, c_axis] = frame->box.getAxis();

    constexpr auto coulomb = 332.063714;
    double e_ele = 0;
    double e_lj = 0;
    for (auto &atom1 : group1) {
        auto atom1_charge = atom1->charge.get();
        auto [sigma1, epsilon1] = atom1->get_lj_p();
        sigma1 *= 10;
        epsilon1 /= 4.184;
        for (auto &atom2 : group2) {
            auto atom2_charge = atom2->charge.get();
            auto [sigma2, epsilon2] = atom2->get_lj_p();
            sigma2 *= 10;
            epsilon2 /= 4.184;

            const auto sigma = 0.5 * (sigma1 + sigma2);
            const auto epsilon = std::sqrt(epsilon1 * epsilon2);

            auto r = atom1->getCoordinate() - atom2->getCoordinate();
            frame->image(r);

            auto r2 = vector_norm2(r);
            auto r6 = r2 * r2 * r2;
            auto r12 = r6 * r6;

            const auto sigma6 = std::pow(sigma, 6);
            const auto sigma12 = sigma6 * sigma6;

            e_lj += 4 * epsilon * (sigma12 / r12 - sigma6 / r6);

            double factor = 1 / vector_norm(r);

            constexpr auto n = 5;
            for (int i = -n; i <= n; ++i) {
                for (int j = -n; j <= n; ++j) {
                    for (int k = -n; k <= n; ++k) {
                        auto r_ = r;
                        std::get<0>(r_) += i * a_axis;
                        std::get<1>(r_) += j * b_axis;
                        std::get<2>(r_) += k * c_axis;

                        auto d = vector_norm(r_);
                        factor += 1 / d;
                    }
                }
            }
            e_ele += coulomb * atom1_charge * atom2_charge * factor;
        }
    }
    return {e_ele, e_lj};
}

void EnergyCalculator::setMask(AmberMask &mask1, AmberMask &mask2, const std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [&](const std::shared_ptr<Atom> &atom) {
        if (is_match(atom, mask1)) {
            group1.push_back(atom);
        } else if (is_match(atom, mask2)) {
            group2.push_back(atom);
        }
    });
}
