//
// Created by xiamr on 10/12/19.
//

#include "CoordinationStructureMatch.hpp"

#include <Eigen/Eigen>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>

#include "data_structure/frame.hpp"
#include "utils/common.hpp"

CoordinationStructureMatch::CoordinationStructureMatch() { enable_outfile = true; }

void CoordinationStructureMatch::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, metal_mask))
            metal = atom;
        else if (Atom::is_match(atom, Ow_atom_mask))
            Ow_atoms.push_back(atom);
    });
    assert(metal);
    assert(Ow_atoms.size() > 2);
}

void CoordinationStructureMatch::process(std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> coord;
    coord.reserve(10);
    for (auto &atom : Ow_atoms) {
        auto r = atom->getCoordinate() - metal->getCoordinate();
        frame->image(r);
        if (vector_norm2(r) < cutoff2) {
            coord.push_back(r);
        }
    }
    if (coord.size() == 9) {
        r_list.emplace_back(testTCTP(coord), testCASP(coord));
    } else {
        r_list.emplace_back(NAN, NAN);
    }
}

void CoordinationStructureMatch::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# metal atom mask > " << metal_mask << '\n';
    os << "# coodination atom mask > " << Ow_atom_mask << '\n';
    os << "# first hydration shell cutoff(Ang) = " << std::sqrt(cutoff2) << '\n';
    os << std::string(50, '#') << '\n';
    std::unordered_map<std::string, std::size_t> statatics;

    os << format("#%15s %15s %15s %15s\n", "frame", "TCTP", "CASP", "Match");
    for (const auto &element : r_list | boost::adaptors::indexed(1)) {
        if (std::isnan(std::get<0>(element.value()))) {
            os << boost::format(" %15d\n") % element.index();
        } else {
            std::string coordStruct = std::get<0>(element.value()) < std::get<1>(element.value()) ? "TCTP" : "CASP";
            os << boost::format(" %15d %15.8f %15.8f %15s\n") % element.index() % std::get<0>(element.value()) %
                      std::get<1>(element.value()) % coordStruct;
            ++statatics[coordStruct];
        }
    }
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s\n", "Coord TYPE", "Count");
    for (const auto &[key, count] : statatics) {
        os << format("%15s %15s\n", key, count);
    }
}

void CoordinationStructureMatch::readInfo() {
    Atom::select2group(metal_mask, Ow_atom_mask, "Enter mask for center metal > ",
                       "Enter mask for coordination atom > ");
    auto cutoff = choose(0.0, "Enter cutoff for first hydration shell (Ang) [ 3.0 ] > ", Default(3.0));
    cutoff2 = cutoff * cutoff;
}

namespace {
struct Constraint {
    std::size_t lhs, rhs;
    std::size_t variable;

    Constraint(std::size_t lhs, std::size_t rhs, std::size_t variable) : lhs(lhs), rhs(rhs), variable(variable) {}

    double length(const std::vector<std::tuple<double, double, double>> &coord) {
        return vector_norm(coord[lhs] - coord[rhs]);
    }
};
}  // namespace

double CoordinationStructureMatch::testCASP(std::vector<std::tuple<double, double, double>> &coord) {
    Eigen::Matrix<double, 5, 1> x0 = Eigen::Matrix<double, 5, 1>::Ones();
    //    Eigen::Matrix<double, 21, 1> R;

    // residuals
    Eigen::Matrix<double, 20, 1> r;

    // Jacobi Matrix
    Eigen::Matrix<double, 20, 5> Dr = Eigen::Matrix<double, 20, 5>::Zero();

    std::vector<Constraint> constraints{
        {0, 1, 0},  // a
        {2, 3, 0},

        {1, 2, 0},  // b
        {0, 3, 0},

        {4, 5, 1},  // c
        {6, 7, 1},

        {5, 6, 1},  // d
        {4, 7, 1},

        {0, 4, 2},  // e
        {1, 5, 2}, {2, 6, 2}, {3, 7, 2},

        {0, 5, 3},  // f
        {1, 6, 3}, {2, 7, 3}, {3, 8, 3},

        {5, 8, 4},  // g
        {6, 8, 4},

        {4, 8, 4},  // h
        {6, 8, 4},
    };

    for (std::size_t i = 0; i < constraints.size(); ++i) {
        Dr(i, constraints[i].variable) = 1.0;
    }

    int iteration = 0;
    double dv;

    boost::sort(coord);

    std::vector<double> residuals;
    residuals.reserve(9 * 8 * 7 * 6 * 5 * 5 * 4 * 3 * 2 * 1);
    do {
        for (;;) {
            ++iteration;
            for (std::size_t i = 0; i < constraints.size(); ++i) {
                r[i] = x0[constraints[i].variable] - constraints[i].length(coord);
            }
            auto DrT = Dr.transpose();
            auto LeftA = DrT * Dr;
            auto b = -DrT * r;
            auto v = LeftA.colPivHouseholderQr().solve(b);
            x0 += v;

            dv = v.norm() / coord.size();
            //            std::cout << format("it = %5d |r| = %15.6f |v| = %f\n",
            //                                iteration, r.norm() / coord.size(), dv);
            if (dv < 1e-8) {
                for (std::size_t i = 0; i < constraints.size(); ++i) {
                    r[i] = x0[constraints[i].variable] - constraints[i].length(coord);
                }
                auto dr = r.norm() / coord.size();
                residuals.push_back(dr);
                break;
            }
        };
    } while (boost::next_permutation(coord));

    //    for (int i = 0; i < 5; ++i) {
    //        std::cout << "x[" << std::to_string(i) << "] = " << x0[i] << '\n';
    //    }
    //    std::cout << format("Total it = %5d |r| = %15.6f |v| = %f\n",
    //                        iteration, r.norm() / coord.size(), dv);

    auto min = std::min_element(std::begin(residuals), std::end(residuals));
    //    std::cout << "min residual = " << *min << '\n';
    //    std::cout << "max residual = " << *max << '\n';
    return *min;
}

double CoordinationStructureMatch::testTCTP(std::vector<std::tuple<double, double, double>> &coord) {
    Eigen::Matrix<double, 5, 1> x0 = Eigen::Matrix<double, 5, 1>::Ones();

    // residuals
    Eigen::Matrix<double, 21, 1> r;

    // Jacobi Matrix
    Eigen::Matrix<double, 21, 5> Dr = Eigen::Matrix<double, 21, 5>::Zero();

    std::vector<Constraint> constraints{{0, 1, 0},  // a
                                        {0, 2, 0}, {3, 4, 0}, {3, 5, 0},

                                        {1, 2, 1},  // b
                                        {4, 5, 1},

                                        {0, 3, 2},  // c
                                        {1, 4, 2}, {2, 5, 2},

                                        {6, 0, 3},  // d
                                        {6, 1, 3}, {6, 3, 3}, {6, 4, 3}, {7, 0, 3}, {7, 2, 3}, {7, 3, 3}, {7, 5, 3},

                                        {8, 1, 4},  // e
                                        {8, 2, 4}, {8, 4, 4}, {8, 5, 4}};

    for (std::size_t i = 0; i < constraints.size(); ++i) {
        Dr(i, constraints[i].variable) = 1.0;
    }

    int iteration = 0;
    double dv;

    boost::sort(coord);

    std::vector<double> residuals;
    residuals.reserve(9 * 8 * 7 * 6 * 5 * 5 * 4 * 3 * 2 * 1);
    do {
        for (;;) {
            ++iteration;
            for (std::size_t i = 0; i < constraints.size(); ++i) {
                r[i] = x0[constraints[i].variable] - constraints[i].length(coord);
            }
            auto DrT = Dr.transpose();
            auto LeftA = DrT * Dr;
            auto b = -DrT * r;
            auto v = LeftA.colPivHouseholderQr().solve(b);
            x0 += v;

            dv = v.norm() / coord.size();
            //            std::cout << format("it = %5d |r| = %15.6f |v| = %f\n",
            //                                iteration, r.norm() / coord.size(), dv);
            if (dv < 1e-8) {
                for (std::size_t i = 0; i < constraints.size(); ++i) {
                    r[i] = x0[constraints[i].variable] - constraints[i].length(coord);
                }
                auto dr = r.norm() / coord.size();
                residuals.push_back(dr);
                break;
            }
        };
    } while (boost::next_permutation(coord));

    //    for (int i = 0; i < 5; ++i) {
    //        std::cout << "x[" << std::to_string(i) << "] = " << x0[i] << '\n';
    //    }
    //    std::cout << format("Total it = %5d |r| = %15.6f |v| = %f\n",
    //                        iteration, r.norm() / coord.size(), dv);

    auto min = std::min_element(std::begin(residuals), std::end(residuals));
    //    std::cout << "min residual = " << *min << '\n';
    //    std::cout << "max residual = " << *max << '\n';
    return *min;
}
