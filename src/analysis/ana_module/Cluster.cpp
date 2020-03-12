//
// Created by xiamr on 6/14/19.
//
#include "Cluster.hpp"

#include <tbb/tbb.h>

#include <boost/container_hash/hash.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <type_traits>

#include "data_structure/frame.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

Cluster::Cluster() {
    enable_tbb = true;
    enable_outfile = true;
}

void Cluster::process(std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> coord;
    coord.reserve(group.size());

    for (const auto &element : group | boost::adaptors::indexed()) {
        if (element.index() == 0) {
            coord.emplace_back(element.value()->x, element.value()->y, element.value()->z);
        } else {
            auto shift = element.value()->getCoordinate() - coord[0];
            frame->image(shift);
            coord.push_back(coord[0] + shift);
        }
    }
    auto center = boost::accumulate(coord, std::tuple<double, double, double>{},
                                    [](const auto &init, const auto &rhs) { return init + rhs; }) /
                  coord.size();
    for (auto &element : coord) {
        element -= center;
    }
    coords.push_back(std::move(coord));
    steps++;
}

double Cluster::rmsvalue(int index1, int index2) {
    int nfit, n;
    nfit = n = group.size();

    double x1[n], y1[n], z1[n];
    double x2[n], y2[n], z2[n];

    for (int i = 0; i < n; i++) {
        std::tie(x1[i], y1[i], z1[i]) = coords[index1][i];
        std::tie(x2[i], y2[i], z2[i]) = coords[index2][i];
    }

    RMSDCal::quatfit(n, x1, y1, z1, n, x2, y2, z2, nfit);
    return RMSDCal::rmsfit(x1, y1, z1, x2, y2, z2, nfit);
}

std::list<Cluster::rmsd_matrix> Cluster::do_calculate_rmsd_list_parallel() {
    class CalCore {
    public:
        std::list<rmsd_matrix> local_rms_list;
        int steps;
        double cutoff;

        Cluster *parent;

        CalCore(int steps, double cutoff, Cluster *parent) : steps(steps), cutoff(cutoff), parent(parent){};

        CalCore(CalCore &c, tbb::split) {
            steps = c.steps;
            cutoff = c.cutoff;
            parent = c.parent;
        }

        void join(CalCore &c) {
            local_rms_list.merge(c.local_rms_list,
                                 [](const rmsd_matrix &m1, const rmsd_matrix &m2) { return (m1.rms < m2.rms); });
        }

        void operator()(const tbb::blocked_range<int> &range) {
            double rms;
            for (int index1 = range.begin(); index1 != range.end(); ++index1) {
                for (int index2 = index1 + 1; index2 < steps; ++index2) {
                    rms = parent->rmsvalue(index1, index2);
                    local_rms_list.emplace_back(index1, index2, rms);
                }
            }
            local_rms_list.sort([](const rmsd_matrix &m1, const rmsd_matrix &m2) { return (m1.rms < m2.rms); });
        }
    } core(steps, cutoff, this);

    tbb::parallel_reduce(tbb::blocked_range<int>(0, steps - 1), core, tbb::auto_partitioner());

    return core.local_rms_list;
}

void Cluster::print(std::ostream &os) {
    auto rmsd_list = do_calculate_rmsd_list_parallel();

    int n = 0;
    for (auto &v : rmsd_list) {
        std::cout << v.i << "  " << v.j << "  " << v.rms << '\n';
        n++;
        if (n > 99) break;
    }

    auto c = do_cluster(rmsd_list, steps, this->cutoff);

    std::cout << "Sorting and renumbering clusters...\n";
    int cid = do_sort_and_renumber_parallel(c);

    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "SET:" << ids << '\n';
    os << "cutoff : " << this->cutoff << '\n';
    os << "Total cluster number : " << cid << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("#%15s %15s\n") % "Frame" % "ClusterNo";
    for (const auto &element : c | boost::adaptors::indexed(1)) {
        os << boost::format(" %15d %15d\n") % element.index() % element.value().clust;
    }

    auto mm = do_find_frames_in_same_clust(c);
    os << "# Clust No.   Count      Frames";
    for (auto i_clust : range(1, cid + 1)) {
        auto &s = mm[i_clust];
        int index = 0;
        auto seq = s | boost::adaptors::transformed([](auto i) { return i + 1; });
        for (const auto &frame : combine_seq(seq)) {
            if (index % 10 == 0) {
                if (index == 0) {
                    os << format("\n%-10d      %-10d ", i_clust, s.size());
                } else {
                    os << '\n' << std::string(27, ' ');
                }
            }
            os << ' ' << frame << ' ';
            index++;
        }
    }
    os << std::string(50, '#') << '\n';

    auto mm2 = do_find_medium_in_clust(c, rmsd_list);
    os << boost::format("#%15s %15s %15s %15s\n") % "Clust No." % "Fame Count" % "Medium_Frame" % "AvgRMSD";
    for (int i_clust : range(1, cid + 1)) {
        os << format(" %15d %15d %15d %15g\n", i_clust, mm[i_clust].size(), mm2[i_clust].first + 1,
                     mm2[i_clust].second);
    }
    os << std::string(50, '#') << '\n';
}

/*
 * Sort and Renumber the cluster vector, make cluster ID start with 1 and return the total cluster amount
 *
 */

int Cluster::do_sort_and_renumber_parallel(std::vector<conf_clust> &conf_clust_vector) {
    do_sort_clust_parallel(conf_clust_vector);
    int total_clust = do_renumber_clust(conf_clust_vector);
    do_sort_conf_parallel(conf_clust_vector);
    return total_clust;
}

int Cluster::do_renumber_clust(std::vector<Cluster::conf_clust> &conf_clust_vector) {
    int cid = 1;
    unsigned int k;
    for (k = 1; k < conf_clust_vector.size(); k++) {
        if (conf_clust_vector[k].clust != conf_clust_vector[k - 1].clust) {
            conf_clust_vector[k - 1].clust = cid;
            cid++;
        } else {
            conf_clust_vector[k - 1].clust = cid;
        }
    }
    conf_clust_vector[k - 1].clust = cid;

    return cid;
}

void Cluster::do_sort_conf_parallel(std::vector<Cluster::conf_clust> &conf_clust_vector) {
    tbb::parallel_sort(conf_clust_vector.begin(), conf_clust_vector.end(),
                       [](const conf_clust &i, const conf_clust &j) { return (i.conf < j.conf); });
}

void Cluster::do_sort_clust_parallel(std::vector<Cluster::conf_clust> &conf_clust_vector) {
    tbb::parallel_sort(conf_clust_vector.begin(), conf_clust_vector.end(),
                       [](const conf_clust &i, const conf_clust &j) { return (i.clust < j.clust); });
}

/*
 *  give the rmsd matrix between frames, return cluster number assgined to frames
 *  the cluster numbers are not have order or continuous
 *
 */

template <typename ForwardRange>
std::vector<Cluster::conf_clust> Cluster::do_cluster(const ForwardRange &rmsd_list, int conf_size, double cutoff) {
    std::vector<Cluster::conf_clust> c = initialize_conf_clust_vector(conf_size);

    // The algorithm of below code block comes from gromacs
    bool bChange;
    do {
        bChange = false;
        for (auto &k : rmsd_list) {
            if (k.rms >= cutoff) break;
            int diff = c[k.j].clust - c[k.i].clust;
            if (diff) {
                bChange = true;
                if (diff > 0) {
                    c[k.j].clust = c[k.i].clust;
                } else {
                    c[k.i].clust = c[k.j].clust;
                }
            }
        }
    } while (bChange);
    /////////////////////////////////////////
    return c;
}

template std::vector<Cluster::conf_clust> Cluster::do_cluster(const std::list<rmsd_matrix> &, int, double);

std::vector<Cluster::conf_clust> Cluster::initialize_conf_clust_vector(int conf_size) {
    std::vector<conf_clust> c;
    for (int i = 0; i < conf_size; i++) {
        c.emplace_back(i, i);
    }
    return c;
}

void Cluster::readInfo() {
    Atom::AmberMask atomIndenter;
    Atom::select1group(atomIndenter, "Please enter group > ");
    setSetting(atomIndenter, choose(0.0, static_cast<double>(std::numeric_limits<int>::max()), "Cutoff > "));
}

void Cluster::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        if (Atom::is_match(atom, this->ids)) this->group.push_back(atom);
    });
}

void Cluster::setSetting(const Atom::AmberMask &atomIndenter, double cutoff) {
    this->ids = atomIndenter;
    throw_assert(cutoff > 0, "cutoff must postive (cutoff = " << cutoff << ")");
    this->cutoff = cutoff;
}

/*
 *
 *
 * frame start from 0
 * clust start from 1
 */

std::unordered_map<int, std::vector<int>> do_find_frames_in_same_clust(const std::vector<Cluster::conf_clust> &clusts) {
    std::unordered_map<int, std::vector<int>> ret;
    for (const auto &ct : clusts) {
        ret[ct.clust].push_back(ct.conf);
    }
    return ret;
}

std::unordered_map<int, std::pair<int, double>> do_find_medium_in_clust(
    const std::vector<Cluster::conf_clust> &clusts, const std::list<Cluster::rmsd_matrix> &rmsd_list) {
    std::unordered_map<int, std::pair<int, double>> ret;
    assert(!clusts.empty());
    std::vector<double> rmsd_sum(clusts.size(), 0.0);
    std::unordered_map<int, std::unordered_set<int>> s;

    std::vector<int> clust_map(clusts.size());

    for (const auto &i : clusts) {
        s[i.clust].insert(i.conf);
        assert(i.conf >= 0 && i.conf < clusts.size());
        clust_map[i.conf] = i.clust;
    }
    for (auto &i : rmsd_list) {
        assert(i.i >= 0 && i.i < clusts.size());
        assert(i.j >= 0 && i.j < clusts.size());
        if (clust_map[i.i] == clust_map[i.j]) {
            rmsd_sum[i.i] += i.rms;
            rmsd_sum[i.j] += i.rms;
        }
    }
    for (auto &j : s) {
        auto min =
            *boost::min_element(j.second, [&rmsd_sum](auto &lhs, auto &rhs) { return rmsd_sum[lhs] < rmsd_sum[rhs]; });
        ret[j.first] = {min, j.second.size() == 1 ? NAN : (rmsd_sum[min] / (j.second.size() - 1))};
    }
    return ret;
}
