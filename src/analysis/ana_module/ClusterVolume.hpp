//
// Created by xiamr on 8/26/19.
//

#ifndef TINKER_CLUSTERVOLUME_HPP
#define TINKER_CLUSTERVOLUME_HPP

#include <tbb/tbb.h>

#include <boost/multi_array.hpp>
#include <utility>

#include "AbstractAnalysis.hpp"
#include "HBond.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "utils/std.hpp"

class ClusterVolume : public AbstractAnalysis {
public:
    ClusterVolume();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "vdW Volume"; }

    enum class ATOM_Category : uint8_t { EMPTY, DEST, OTHER };

    using argument_type = std::tuple<boost::multi_array<ATOM_Category, 3> *, boost::multi_array<double, 2> *,
                                     boost::multi_array<double, 2> *, double, double, double, double, int>;

    argument_type preprocess(std::shared_ptr<Frame> &frame);

    void work_body(boost::multi_array<ATOM_Category, 3> *grid, boost::multi_array<double, 2> *atom_group_array,
                   boost::multi_array<double, 2> *other_atom_array, double grid_x_step, double grid_y_step,
                   double grid_z_step, double total_volume, int nframe);

protected:
    bool enable_paralel_while_impl() override;

    void do_parallel_while_impl(std::function<std::shared_ptr<Frame>()> func) override;

    static double getVdwRadii(const std::shared_ptr<Atom> &atom);

    void fill_atom(boost::multi_array<ATOM_Category, 3> *grid, double grid_x_step, double grid_y_step,
                   double grid_z_step, ATOM_Category category, double x, double y, double z, double radii) const;

    std::vector<std::tuple<int, int, int>> generate_neighbor_grids(double grid_x_step, double grid_y_step,
                                                                   double grid_z_step) const;

    AmberMask atom_mask;

    std::unordered_set<std::shared_ptr<Atom>> atom_group, other_atoms;

    int grid_x;
    int grid_y;
    int grid_z;

    tbb::concurrent_hash_map<int, std::tuple<double, double, double, double>> volumes;

    inline static std::unordered_map<Symbol, double> vdWRadiis{
        {Symbol::Hydrogen, 1.10},   {Symbol::Carbon, 1.70}, {Symbol::Nitrogen, 1.55}, {Symbol::Oxygen, 1.52},
        {Symbol::Phosphorus, 1.80}, {Symbol::Sulfur, 1.80}, {Symbol::Sodium, 2.27}};

    size_t countFilledGridPoints(boost::multi_array<ATOM_Category, 3> *grid) const;

    bool fill_space(boost::multi_array<ATOM_Category, 3> *grid, double grid_x_step, double grid_y_step,
                    double grid_z_step) const;

    std::vector<double> radii_for_atom_group, radii_for_other_atoms;

    std::pair<std::size_t, std::size_t> do_grid(boost::multi_array<ATOM_Category, 3> *grid,
                                                boost::multi_array<double, 2> *atom_group_array,
                                                boost::multi_array<double, 2> *other_atom_array, double grid_x_step,
                                                double grid_y_step, double grid_z_step) const;

    int current_frame_num = 0;

    class FrameStream {
        std::function<std::shared_ptr<Frame>()> func;
        ClusterVolume *parent;

    public:
        FrameStream(std::function<std::shared_ptr<Frame>()> func, ClusterVolume *parent)
            : func(std::move(func)), parent(parent) {}

        bool pop_if_present(ClusterVolume::argument_type &item);
    };

    friend class FrameStream;

    class ApplyBody {
        ClusterVolume *parent;

    public:
        explicit ApplyBody(ClusterVolume *parent) : parent(parent) {}

        using argument_type = ClusterVolume::argument_type;

        void operator()(argument_type &item) const;
    };

    friend class ApplyBody;
};

#endif  // TINKER_CLUSTERVOLUME_HPP
