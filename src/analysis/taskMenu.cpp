//
// Created by xiamr on 6/14/19.
//

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 100

#if (BOOST_MPL_LIMIT_VECTOR_SIZE <= 50)
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#include <boost/mpl/vector.hpp>
#else

#include <boost/preprocessor/iterate.hpp>
#include <boost/mpl/vector/vector50.hpp>

namespace boost {
    namespace mpl {
#define BOOST_PP_ITERATION_PARAMS_1 \
        (3,(51, BOOST_MPL_LIMIT_VECTOR_SIZE, <boost/mpl/vector/aux_/numbered.hpp>))

#include BOOST_PP_ITERATE()

    } // namespace mpl
} // namespace boost

#define BOOST_MPL_PREPROCESSING_MODE

#include <boost/mpl/vector.hpp>

#undef BOOST_MPL_PREPROCESSING_MODE
#endif

#include <boost/phoenix.hpp>

// Boost metaprogramming library
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/string.hpp>

namespace mpl = boost::mpl;

#include "taskMenu.hpp"

#include "Trajconv.hpp"
#include "Distance.hpp"
#include "CoordinateNumPerFrame.hpp"
#include "FirstCoordExchangeSearch.hpp"
#include "RadicalDistribtuionFunction.hpp"
#include "RMSDCal.hpp"
#include "RMSFCal.hpp"
#include "HBond.hpp"
#include "ResidenceTime.hpp"
#include "GreenKubo.hpp"
#include "Cluster.hpp"
#include "NMRRange.hpp"
#include "RotAcfCutoff.hpp"
#include "DiffuseCutoff.hpp"
#include "Diffuse.hpp"
#include "DipoleAngle.hpp"
#include "DipoleAngleSingleDistanceNormal.hpp"
#include "DipoleAngle2Gibbs.hpp"
#include "DipoleAngleVolumeNormal.hpp"
#include "ShellDensity.hpp"
#include "SearchInteractionResidue.hpp"
#include "FindMinBetweenTwoGroups.hpp"
#include "DemixIndexOfTwoGroup.hpp"
#include "DipoleAngleWithDistanceRange.hpp"
#include "DipoleAxisDistribution.hpp"
#include "EquatorialAngle.hpp"
#include "PlaneAngle.hpp"
#include "AngleWat.hpp"
#include "RotAcf.hpp"
#include "DistanceAngle.hpp"
#include "DipoleAngleAxis3D.hpp"
#include "DynFrameFind.hpp"
#include "SpatialOrientationDistribution.hpp"
#include "AngleDistributionBetweenTwoVectorWithCutoff.hpp"
#include "HBondLifeTime.hpp"
#include "HBondLifeTimeCutoff.hpp"
#include "HBondSpread.hpp"
#include "IRSpectrum.hpp"
#include "VelocityAutocorrelationFunction.hpp"
#include "IRSpectrumElectricalFlux.hpp"
#include "DensityOfStates.hpp"
#include "ConvertVelocityToVelocityCharge.hpp"
#include "IRSpectrumDeltaDipole.hpp"
#include "ClusterVolume.hpp"

using namespace std;

template<typename T>
struct boost_type {
    typedef T type;

    boost_type(boost::type<T>) {};
};

std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> getTasks() {
    auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();

    using components = mpl::vector<
            Trajconv,
            Distance,
            CoordinateNumPerFrame,
            RadicalDistribtuionFunction,
            ResidenceTime,
            GreenKubo,
            HBond,
            RMSDCal,
            RMSFCal,
            Cluster,
            NMRRange,
            Diffuse,
            DiffuseCutoff,
            FirstCoordExchangeSearch,
            RotAcfCutoff,
            DipoleAngle,
            DipoleAngle2Gibbs,
            DipoleAngleSingleDistanceNormal,
            DipoleAngleVolumeNormal,
            ShellDensity,
            SearchInteractionResidue,
            FindMinBetweenTwoGroups,
            DemixIndexOfTwoGroup,
            DipoleAngleWithDistanceRange,
            DipoleAxisDistribution,
            EquatorialAngle,
            PlaneAngle,
            AngleWat,
            RotAcf,
            DistanceAngle,
            DipoleAngleAxis3D,
            DynFrameFind,
            SpatialOrientationDistribution,
            AngleDistributionBetweenTwoVectorWithCutoff,
            HBondLifeTime,
            HBondLifeTimeCutoff,
            HBondSpread,
            IRSpectrum,
            VelocityAutocorrelationFunction,
            IRSpectrumElectricalFlux,
            DensityOfStates,
            ConvertVelocityToVelocityCharge,
            IRSpectrumDeltaDipole,
            ClusterVolume
    >;

    BOOST_MPL_ASSERT((mpl::equal<mpl::unique<components, is_same<mpl::_1, mpl::_2> >::type, components>));

    std::vector<std::function<shared_ptr<BasicAnalysis>()>> task_vec;
    std::vector<string> item_menu;

    mpl::for_each<components, boost::type<mpl::_>>([&task_vec, &item_menu](auto t) {
        using T = typename decltype(boost_type(t))::type;
        task_vec.emplace_back(bind(make_shared<T>));
        item_menu.emplace_back((boost::format("(%d) %s") % task_vec.size() % T::title()).str());
    });

    auto menu1 = [&item_menu]() {
        std::cout << "Please select the desired operation (Trajectory Analysis)" << std::endl;
        std::cout << "(0) Start\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << boost::phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<components>::value, "select : ");
    };
    while (true) {
        int num = menu1();
        if (num == 0) return task_list;
        shared_ptr<BasicAnalysis> task = task_vec[num - 1]();

        string line(item_menu[num - 1].size() + 6, '-');

        std::cout << line << "\n";
        std::cout << "<- " << item_menu[num - 1] << " ->\n";
        std::cout << line << "\n";
        task->readInfo();
        task_list->push_back(task);
    }
}
