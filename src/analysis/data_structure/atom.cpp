//
// Created by xiamr on 3/17/19.
//

#include "atom.hpp"

std::array<double, 2> Atom::get_lj_p() const {
    const auto c6 = (*lj_param).c6;
    const auto c12 = (*lj_param).c12;
    const auto sigma = c6 != 0.0 ? std::pow(c12 / c6, 1.0 / 6) : 0.0;
    const auto epsilon = c12 != 0.0 ? c6 * c6 / (4 * c12) : 0.0;
    return {sigma, epsilon};
}

bool Atom::adj(const std::shared_ptr<Atom> &atom) {
    for (auto i : con_list) {
        if (atom->seq == i)
            return true;
    }
    return false;
}
