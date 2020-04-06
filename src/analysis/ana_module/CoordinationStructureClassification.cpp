//
// Created by xiamr on 9/20/19.
//

#include "CoordinationStructureClassification.hpp"

#include <tbb/tbb.h>
#include <tbb/tick_count.h>

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>

#include "RMSDCal.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

CoordinationStructureClassification::CoordinationStructureClassification() {
    enable_outfile = true;
    enable_tbb = true;
}

void CoordinationStructureClassification::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (is_match(atom, metal_mask))
            metal = atom;
        else if (is_match(atom, Ow_atom_mask))
            Ow_atoms.push_back(atom);
    });
    assert(metal);
    assert(Ow_atoms.size() > 2);
}

void CoordinationStructureClassification::process(std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> coord;
    for (auto &atom : Ow_atoms) {
        auto r = atom->getCoordinate() - metal->getCoordinate();
        frame->image(r);
        if (vector_norm2(r) < cutoff2) {
            coord.push_back(r);
        }
    }
    frame_cn_mapping[nframe] = coord.size();
    systems[coord.size()].emplace_back(nframe, coord);
    nframe++;
}

void CoordinationStructureClassification::print(std::ostream &os) {
    auto rmsd_list_map = do_calculate_rmsd_list_parallel();

    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# metal atom mask > " << metal_mask << '\n';
    os << "# coodination atom mask > " << Ow_atom_mask << '\n';
    os << "# first hydration shell cutoff(Ang) = " << std::sqrt(cutoff2) << '\n';
    os << "# rmsd cutoff(Ang) = " << rmsd_cutoff << '\n';
    os << std::string(50, '#') << '\n';

    std::list<Cluster::rmsd_matrix> rmsd_list;

    boost::for_each(rmsd_list_map, [&rmsd_list](auto &element) {
        rmsd_list.merge(element.second, [](const Cluster::rmsd_matrix &m1, const Cluster::rmsd_matrix &m2) {
            return (m1.rms < m2.rms);
        });
    });

    if (output_rms_matrix) {
        for (auto &element : rmsd_list) {
            os << format("%10d %10d %10.5f\n", element.i, element.j, element.rms);
        }
    }

    std::vector<Cluster::conf_clust> c = Cluster::do_cluster(rmsd_list, nframe, rmsd_cutoff);

    int cid = Cluster::do_sort_and_renumber_parallel(c);
    os << "# Total cluster number : " << cid << '\n';
    os << std::string(50, '#') << '\n';
    os << format("#%15s %15s %15s\n", "Frame", "C.N.", "Clust No.");
    for (auto k : boost::irange<int>(0, c.size())) {
        os << format(" %15d %15d %15d\n", k + 1, frame_cn_mapping[k], c[k].clust);
    }
    os << std::string(50, '#') << '\n';

    // do statistics
    std::unordered_map<int, std::vector<int>> mm = do_find_frames_in_same_clust(c);
    os << format("#%-10s %-5s %-13s %-13s %-13s %-13s", "Clust No.", "C.N.", "Frame Count", "MediumFrame", "AvgRMSD",
                 "Frames enumeration");
    std::unordered_map<int, std::pair<int, double>> mm2 = do_find_medium_in_clust(c, rmsd_list);
    for (auto i_clust : range(1, cid + 1)) {
        auto &s = mm[i_clust];
        std::size_t string_length = 0;
        os << format("\n %-10d %-5d %-13d %-13d %-13g ", i_clust, frame_cn_mapping[s.front()], s.size(),
                     mm2[i_clust].first + 1, mm2[i_clust].second);
        auto seq = s | boost::adaptors::transformed([](auto i) { return i + 1; });
        for (const auto &frame : combine_seq(seq)) {
            if (string_length > 80) {
                os << '\n' << std::string(60, ' ');
                string_length = 0;
            }
            os << ' ' << frame << ' ';
            string_length += frame.length() + 2;
        }
    }
    os << '\n';
    os << std::string(50, '#') << '\n';
}

std::map<int, std::list<Cluster::rmsd_matrix>> CoordinationStructureClassification::do_calculate_rmsd_list_parallel() {
    class Body {
    public:
        std::list<Cluster::rmsd_matrix> local_rms_list;
        std::deque<std::pair<int, std::vector<std::tuple<double, double, double>>>> &tables;
        CoordinationStructureClassification *parent;
        std::atomic<std::size_t> &current_completed_compute_amount;
        std::size_t &total_compute_amount;

        Body(std::deque<std::pair<int, std::vector<std::tuple<double, double, double>>>> &tables,
             CoordinationStructureClassification *parent, std::atomic<std::size_t> &current_completed_compute_amount,
             std::size_t &total_compute_amount)
            : tables(tables), parent(parent), current_completed_compute_amount(current_completed_compute_amount),
              total_compute_amount(total_compute_amount){};

        Body(Body &c, tbb::split)
            : tables(c.tables), parent(c.parent), current_completed_compute_amount(c.current_completed_compute_amount),
              total_compute_amount(c.total_compute_amount) {}

        void join(Body &c) {
            local_rms_list.merge(c.local_rms_list, [](const Cluster::rmsd_matrix &m1, const Cluster::rmsd_matrix &m2) {
                return (m1.rms < m2.rms);
            });
        }

        void operator()(const tbb::blocked_range<std::size_t> &range) {
            auto report_amount = total_compute_amount / 10;
            std::size_t current_com_amount = 0;
            for (std::size_t index1 = range.begin(); index1 != range.end(); ++index1) {
                for (std::size_t index2 = index1 + 1; index2 < tables.size(); ++index2) {
                    auto &item1 = tables[index1];
                    auto &item2 = tables[index2];
                    local_rms_list.emplace_back(
                        item1.first, item2.first,
                        CoordinationStructureClassification::calculateRmsdOfTwoStructs(item1.second, item2.second));
                }
                current_com_amount += tables.size() - index1 - 1;
                if (current_com_amount >= report_amount or (index1 + 1) == range.end()) {
                    current_completed_compute_amount += current_com_amount;
                    std::cout << "\rTBB parallel block Complete " << std::setw(3)
                              << current_completed_compute_amount / total_compute_amount << " %    " << std::flush;
                    current_com_amount = 0;
                }
            }
            local_rms_list.sort(
                [](const Cluster::rmsd_matrix &m1, const Cluster::rmsd_matrix &m2) { return (m1.rms < m2.rms); });
        }
    };

    std::map<int, Body> bodys;

    std::size_t total_compute_amount = 0;
    std::atomic<std::size_t> current_completed_compute_amount = 0;

    tbb::task_group taskGroup;
    for (auto &element : systems) {
        bodys.insert(
            {element.first, Body(element.second, this, current_completed_compute_amount, total_compute_amount)});
        auto b = &bodys.at(element.first);
        taskGroup.run(
            [b, &element] { tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, element.second.size()), *b); });
        total_compute_amount += element.second.size() * (element.second.size() - 1) / 2;
    }
    total_compute_amount /= 100;
    auto start_time = std::chrono::steady_clock::now();
    std::cout << "\rTBB parallel block Complete   0 %    " << std::flush;
    taskGroup.wait();
    std::cout << "\rTBB parallel block Complete 100 %    elapsed time "
              << chrono_cast(std::chrono::steady_clock::now() - start_time) << std::endl;

    std::map<int, std::list<Cluster::rmsd_matrix>> rmsd_list_map;

    for (auto &b : bodys) {
        rmsd_list_map.insert({b.first, std::move(b.second.local_rms_list)});
    }

    return rmsd_list_map;
}

namespace {
int rms_max_index(double x1[], double y1[], double z1[], double x2[], double y2[], double z2[], int n_rms_calc) {
    double xr, yr, zr, dist2;
    double rms2 = 0.0;
    int index = 0;
    for (int i = 0; i < n_rms_calc; ++i) {
        xr = x1[i] - x2[i];
        yr = y1[i] - y2[i];
        zr = z1[i] - z2[i];
        dist2 = xr * xr + yr * yr + zr * zr;
        if (dist2 > rms2) {
            index = i;
            rms2 = dist2;
        }
    }
    return index;
}
} // namespace

double
CoordinationStructureClassification::calculateRmsdOfTwoStructs(std::vector<std::tuple<double, double, double>> &c1,
                                                               std::vector<std::tuple<double, double, double>> &c2) {
    assert(c1.size() == c2.size());
    assert(c1.size() > 2);

    double x1[c1.size() + 1], y1[c1.size() + 1], z1[c1.size() + 1];
    double x2[c1.size() + 1], y2[c1.size() + 1], z2[c1.size() + 1];

    auto [c1_index1, c1_index2] = find_min_distance_pair(c1);

    double rmsd_value2 = std::numeric_limits<double>::max();

    for (auto c2_index1 : boost::irange<int>(0, c2.size())) {
        for (auto c2_index2 : boost::irange<int>(0, c2.size())) {
            if (c2_index1 == c2_index2)
                continue;

            x1[0] = y1[0] = z1[0] = x2[0] = y2[0] = z2[0] = 0.0;

            fill_coord(x1, y1, z1, c1_index1, c1_index2, c1);
            fill_coord(x2, y2, z2, c2_index1, c2_index2, c2);

            // superpose with 3 atoms
            RMSDCal::quatfit(c1.size() + 1, x1, y1, z1, c1.size() + 1, x2, y2, z2, 3);

            //排列index >=3 之后的原子，使相对应的距离最小
            permutation(x1, y1, z1, x2, y2, z2, 3, c1.size() + 1);

            RMSDCal::quatfit(c1.size() + 1, x1, y1, z1, c1.size() + 1, x2, y2, z2, c1.size() + 1);

            auto rmax_index = rms_max_index(x1, y1, z1, x2, y2, z2, c1.size() + 1);

            std::swap(x1[rmax_index], x1[c1.size()]);
            std::swap(y1[rmax_index], y1[c1.size()]);
            std::swap(z1[rmax_index], z1[c1.size()]);

            std::swap(x2[rmax_index], x2[c1.size()]);
            std::swap(y2[rmax_index], y2[c1.size()]);
            std::swap(z2[rmax_index], z2[c1.size()]);

            RMSDCal::quatfit(c1.size() + 1, x1, y1, z1, c1.size() + 1, x2, y2, z2, c1.size());
            auto xr = x1[c1.size()] - x2[c1.size()];
            auto yr = y1[c1.size()] - y2[c1.size()];
            auto zr = z1[c1.size()] - z2[c1.size()];

            rmsd_value2 = std::min(rmsd_value2, xr * xr + yr * yr + zr * zr);
        }
    }
    return std::sqrt(rmsd_value2);
}

std::pair<int, int>
CoordinationStructureClassification::find_min_distance_pair(std::vector<std::tuple<double, double, double>> &c) {
    double min_distance2 = std::numeric_limits<double>::max();
    int index1, index2;
    index1 = index2 = std::numeric_limits<int>::max();

    for (auto i : boost::irange<int>(0, c.size() - 1)) {
        for (auto j : boost::irange<int>(i + 1, c.size())) {
            auto distance2 = vector_norm2(c[i] - c[j]);
            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                index1 = i;
                index2 = j;
            }
        }
    }
    return {index1, index2};
}

void CoordinationStructureClassification::fill_coord(double x[], double y[], double z[], int index1, int index2,
                                                     std::vector<std::tuple<double, double, double>> &c) {
    std::tie(x[1], y[1], z[1]) = c[index1];
    std::tie(x[2], y[2], z[2]) = c[index2];

    int index = 3;
    for (auto i : boost::irange<int>(0, c.size())) {
        if (i != index1 and i != index2) {
            std::tie(x[index], y[index], z[index]) = c[i];
            index++;
        }
    }
}

void CoordinationStructureClassification::permutation(double x1[], double y1[], double z1[], double x2[], double y2[],
                                                      double z2[], int start, int end /*not included*/) {
    for (auto i : boost::irange(start, end)) {
        double min_distance2 = std::numeric_limits<double>::max();
        int min_index;
        for (auto j : boost::irange(i, end)) {
            auto xr = x1[i] - x2[j];
            auto yr = y1[i] - y2[j];
            auto zr = z1[i] - z2[j];

            auto distance2 = xr * xr + yr * yr + zr * zr;

            if (distance2 < min_distance2) {
                min_distance2 = distance2;
                min_index = j;
            }
        }
        if (i != min_index) {
            std::swap(x2[i], x2[min_index]);
            std::swap(y2[i], y2[min_index]);
            std::swap(z2[i], z2[min_index]);
        }
    }
}

void CoordinationStructureClassification::readInfo() {
    select2group(metal_mask, Ow_atom_mask, "Enter mask for center metal > ", "Enter mask for coordination atom > ");
    auto cutoff = choose(0.0, "Enter cutoff for first hydration shell (Ang) [ 3.0 ] > ", Default(3.0));
    cutoff2 = cutoff * cutoff;

    rmsd_cutoff = choose(0.0, "Enter cutoff for RMSD (Ang) > ");

    output_rms_matrix = choose_bool("Output RMS Matrix [N] > ", Default(false));
}
