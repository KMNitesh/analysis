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
#include "OrientationResolvedRadialDistributionFunction.hpp"
#include "ConditionalTimeCorrelationFunction.hpp"
#include "RadiusOfGyration.hpp"
#include "HBondLifeTimeContinuous.hpp"
#include "HBondLifeTimeCutoffContinuous.hpp"
#include "SspResidenceTime.hpp"
#include "CoordinationStructureClassification.hpp"
#include "LocalStructureIndex.hpp"
#include "LocalStructureIndexForLiquid.hpp"

using namespace std;
namespace {
    template<typename T>
    struct boost_type {
        typedef T type;

        boost_type(boost::type<T>) {};
    };

    template<typename components>
    std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> subMenu(const std::string &title) {

        auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();

        BOOST_MPL_ASSERT((mpl::equal<typename mpl::unique<components, is_same<mpl::_1, mpl::_2> >::type, components>));

        std::vector<std::function<shared_ptr<BasicAnalysis>()>> task_vec;
        std::vector<string> item_menu;

        mpl::for_each<components, boost::type<mpl::_>>([&task_vec, &item_menu](auto t) {
            using T = typename decltype(boost_type(t))::type;
            task_vec.emplace_back(bind(make_shared<T>));
            item_menu.emplace_back((boost::format("(%d) %s") % task_vec.size() % T::title()).str());
        });

        auto menu1 = [&item_menu, &title]() {
            std::cout << title << std::endl;
            std::cout << "(0) Return to Main Menu\n";
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
}


std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> getTasks() {
    auto task_list = make_shared<list<shared_ptr<BasicAnalysis>>>();

    using angleDistributionMenu = mpl::vector<
            DipoleAngle,
            DipoleAngle2Gibbs,
            DipoleAngleSingleDistanceNormal,
            DipoleAngleVolumeNormal,
            DipoleAngleWithDistanceRange,
            DipoleAxisDistribution,
            EquatorialAngle,
            PlaneAngle,
            AngleWat,
            DistanceAngle,
            DipoleAngleAxis3D,
            SpatialOrientationDistribution,
            AngleDistributionBetweenTwoVectorWithCutoff
    >;

    using structurePropertyMenu = mpl::vector<
            Distance,
            CoordinateNumPerFrame,
            RadicalDistribtuionFunction,
            RMSDCal,
            RMSFCal,
            Cluster,
            ShellDensity,
            OrientationResolvedRadialDistributionFunction,
            RadiusOfGyration,
            CoordinationStructureClassification,
            LocalStructureIndex,
            LocalStructureIndexForLiquid
    >;

    using hydrogenBondMenu  = mpl::vector<
            HBond,
            HBondLifeTime,
            HBondLifeTimeContinuous,
            HBondLifeTimeCutoff,
            HBondLifeTimeCutoffContinuous,
            HBondSpread
    >;

    using dynamicsPropertyMenu = mpl::vector<
            ResidenceTime,
            SspResidenceTime,
            Diffuse,
            GreenKubo,
            DiffuseCutoff,
            VelocityAutocorrelationFunction,
            RotAcf,
            RotAcfCutoff,
            ConditionalTimeCorrelationFunction
    >;

    using biphaseSystemMenu = mpl::vector<
            DemixIndexOfTwoGroup,
            ClusterVolume
    >;


    using spectraMenu = mpl::vector<
            IRSpectrum,
            IRSpectrumElectricalFlux,
            IRSpectrumDeltaDipole,
            DensityOfStates
    >;

    using trajectoryTransformationMenu = mpl::vector<
            Trajconv,
            ConvertVelocityToVelocityCharge
    >;

    using noeMenu = mpl::vector<NMRRange>;

    using otherUtilsMenu = mpl::vector<
            DynFrameFind,
            SearchInteractionResidue,
            FindMinBetweenTwoGroups,
            FirstCoordExchangeSearch
    >;

    using mainMenu = mpl::vector<
            trajectoryTransformationMenu,
            structurePropertyMenu,
            dynamicsPropertyMenu,
            angleDistributionMenu,
            hydrogenBondMenu,
            biphaseSystemMenu,
            spectraMenu,
            noeMenu,
            otherUtilsMenu
    >;

    std::vector<std::string> menuString{
            "Trajectory Transformation",
            "Structure Properties",
            "Dynamic Properties",
            "Various Angle Distribution",
            "Hydrogen Bond",
            "Biphase System",
            "Spetra",
            "NMR",
            "Other Utils"
    };

    std::vector<std::function<std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>>()>> task_vec;
    std::vector<string> item_menu;

    mpl::for_each<mainMenu, boost::type<mpl::_>>([&task_vec, &item_menu, &menuString](auto t) {
        using T = typename decltype(boost_type(t))::type;
        auto title = menuString[task_vec.size()];
        task_vec.emplace_back([title] { return subMenu<T>(title); });
        item_menu.emplace_back((boost::format("(%d) %s") % task_vec.size() % title).str());
    });

    auto menu1 = [&item_menu]() {
        std::cout << "--- Trajectory Analysis ---" << std::endl;
        std::cout << "(0) Start...\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << boost::phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<mainMenu>::value, "select : ");
    };

    while (true) {
        int num = menu1();
        if (num == 0) return task_list;
        auto tasks = task_vec[num - 1]();
        task_list->splice(task_list->end(), *tasks);
    }
}
