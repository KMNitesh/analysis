
#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3

#include <boost/phoenix.hpp>

// Boost metaprogramming library
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/vector.hpp>

namespace mpl = boost::mpl;

#include "ana_module/Angle.hpp"
#include "ana_module/AngleDistributionBetweenTwoVectorWithCutoff.hpp"
#include "ana_module/AngleWat.hpp"
#include "ana_module/Cluster.hpp"
#include "ana_module/ClusterVolume.hpp"
#include "ana_module/ConditionalTimeCorrelationFunction.hpp"
#include "ana_module/ConvertVelocityToVelocityCharge.hpp"
#include "ana_module/CoordinateNumPerFrame.hpp"
#include "ana_module/CoordinationStructureClassification.hpp"
#include "ana_module/CoordinationStructureMatch.hpp"
#include "ana_module/CoplaneIndex.hpp"
#include "ana_module/DemixIndexOfTwoGroup.hpp"
#include "ana_module/DensityOfStates.hpp"
#include "ana_module/Diffuse.hpp"
#include "ana_module/DiffuseCutoff.hpp"
#include "ana_module/DipoleAngle.hpp"
#include "ana_module/DipoleAngle2Gibbs.hpp"
#include "ana_module/DipoleAngleAxis3D.hpp"
#include "ana_module/DipoleAngleSingleDistanceNormal.hpp"
#include "ana_module/DipoleAngleVolumeNormal.hpp"
#include "ana_module/DipoleAngleWithDistanceRange.hpp"
#include "ana_module/DipoleAxisDistribution.hpp"
#include "ana_module/Distance.hpp"
#include "ana_module/Distance2Plane.hpp"
#include "ana_module/DistanceAngle.hpp"
#include "ana_module/DynFrameFind.hpp"
#include "ana_module/EquatorialAngle.hpp"
#include "ana_module/FindMinBetweenTwoGroups.hpp"
#include "ana_module/FirstCoordExchangeSearch.hpp"
#include "ana_module/GreenKubo.hpp"
#include "ana_module/HBond.hpp"
#include "ana_module/HBondLifeTime.hpp"
#include "ana_module/HBondLifeTimeContinuous.hpp"
#include "ana_module/HBondLifeTimeCutoff.hpp"
#include "ana_module/HBondLifeTimeCutoffContinuous.hpp"
#include "ana_module/HBondSpread.hpp"
#include "ana_module/IRSpectrum.hpp"
#include "ana_module/IRSpectrumDeltaDipole.hpp"
#include "ana_module/IRSpectrumElectricalFlux.hpp"
#include "ana_module/LocalStructureIndex.hpp"
#include "ana_module/LocalStructureIndexForLiquid.hpp"
#include "ana_module/NMRRange.hpp"
#include "ana_module/OrientationResolvedRadialDistributionFunction.hpp"
#include "ana_module/PlaneAngle.hpp"
#include "ana_module/RMSDCal.hpp"
#include "ana_module/RMSFCal.hpp"
#include "ana_module/RadicalDistribtuionFunction.hpp"
#include "ana_module/RadiusOfGyration.hpp"
#include "ana_module/ResidenceTime.hpp"
#include "ana_module/RotAcf.hpp"
#include "ana_module/RotAcfCutoff.hpp"
#include "ana_module/SearchInteractionResidue.hpp"
#include "ana_module/ShellDensity.hpp"
#include "ana_module/SpatialOrientationDistribution.hpp"
#include "ana_module/SspResidenceTime.hpp"
#include "ana_module/Trajconv.hpp"
#include "ana_module/VelocityAutocorrelationFunction.hpp"
#include "taskMenu.hpp"

using namespace std;
namespace {

template <typename components>
std::shared_ptr<std::list<std::shared_ptr<AbstractAnalysis>>> subMenu(const std::string &title) {
    auto task_list = make_shared<list<shared_ptr<AbstractAnalysis>>>();

    BOOST_MPL_ASSERT((mpl::equal<typename mpl::unique<components, is_same<mpl::_1, mpl::_2>>::type, components>));

    std::vector<std::function<shared_ptr<AbstractAnalysis>()>> task_vec;
    std::vector<string> item_menu;

    mpl::for_each<components, boost::type<mpl::_>>([&task_vec, &item_menu ]<typename T>(boost::type<T>) {
        task_vec.emplace_back(bind(make_shared<T>));
        item_menu.emplace_back((boost::format("(%d) %s") % task_vec.size() % T::title()).str());
    });

    auto menu1 = [&item_menu, &title]() {
        std::cout << title << std::endl;
        std::cout << "(0) Return to Main Menu\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << boost::phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<components>::value, "select : ");
    };
    int num = menu1();
    if (num == 0) return task_list;
    shared_ptr<AbstractAnalysis> task = task_vec[num - 1]();

    string line(item_menu[num - 1].size() + 6, '-');

    std::cout << line << "\n";
    std::cout << "<- " << item_menu[num - 1] << " ->\n";
    std::cout << line << "\n";
    task->readInfo();
    task_list->push_back(task);
    return task_list;
}
}  // namespace

std::shared_ptr<std::list<std::shared_ptr<AbstractAnalysis>>> getTasks() {
    auto task_list = make_shared<list<shared_ptr<AbstractAnalysis>>>();

    using angleDistributionMenu =
        mpl::vector<DipoleAngle, DipoleAngle2Gibbs, DipoleAngleSingleDistanceNormal, DipoleAngleVolumeNormal,
                    DipoleAngleWithDistanceRange, DipoleAxisDistribution, EquatorialAngle, PlaneAngle, AngleWat,
                    DistanceAngle, DipoleAngleAxis3D, SpatialOrientationDistribution,
                    AngleDistributionBetweenTwoVectorWithCutoff, CoplaneIndex>;

    using structurePropertyMenu =
        mpl::vector<Distance, Angle, Distance2Plane, CoordinateNumPerFrame, RadicalDistribtuionFunction, RMSDCal,
                    RMSFCal, Cluster, ShellDensity, OrientationResolvedRadialDistributionFunction, RadiusOfGyration,
                    CoordinationStructureClassification, CoordinationStructureMatch, LocalStructureIndex,
                    LocalStructureIndexForLiquid>;

    using hydrogenBondMenu = mpl::vector<HBond, HBondLifeTime, HBondLifeTimeContinuous, HBondLifeTimeCutoff,
                                         HBondLifeTimeCutoffContinuous, HBondSpread>;

    using dynamicsPropertyMenu =
        mpl::vector<ResidenceTime, SspResidenceTime, Diffuse, GreenKubo, DiffuseCutoff, VelocityAutocorrelationFunction,
                    RotAcf, RotAcfCutoff, ConditionalTimeCorrelationFunction>;

    using biphaseSystemMenu = mpl::vector<DemixIndexOfTwoGroup, ClusterVolume>;

    using spectraMenu = mpl::vector<IRSpectrum, IRSpectrumElectricalFlux, IRSpectrumDeltaDipole, DensityOfStates>;

    using trajectoryTransformationMenu = mpl::vector<Trajconv, ConvertVelocityToVelocityCharge>;

    using nmrMenu = mpl::vector<NMRRange>;

    using otherUtilsMenu =
        mpl::vector<DynFrameFind, SearchInteractionResidue, FindMinBetweenTwoGroups, FirstCoordExchangeSearch>;

    using mainMenu =
        mpl::vector<trajectoryTransformationMenu, structurePropertyMenu, dynamicsPropertyMenu, angleDistributionMenu,
                    hydrogenBondMenu, biphaseSystemMenu, spectraMenu, nmrMenu, otherUtilsMenu>;

    std::vector<std::string> menuString{"Trajectory Transformation",
                                        "Structure Properties",
                                        "Dynamic Properties",
                                        "Various Angle Distribution",
                                        "Hydrogen Bond",
                                        "Biphase System",
                                        "Spetra",
                                        "NMR",
                                        "Other Utils"};

    std::vector<std::function<std::shared_ptr<std::list<std::shared_ptr<AbstractAnalysis>>>()>> menu_functions;
    std::vector<string> item_menu;

    mpl::for_each<mainMenu, boost::type<mpl::_>>(
        [&menu_functions, &item_menu, &menuString ]<typename T>(boost::type<T>) {
            auto title = menuString[menu_functions.size()];
            menu_functions.emplace_back([title] { return subMenu<T>(title); });
            item_menu.emplace_back((boost::format("(%d) %s") % menu_functions.size() % title).str());
        });

    auto menu1 = [&item_menu]() {
        std::cout << "--- Trajectory Analysis ---" << std::endl;
        std::cout << "(0) Start...\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << boost::phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<mainMenu>::value, "select : ");
    };

    for (;;) {
        int num = menu1();
        if (num == 0) {
            if (task_list->empty()) {
                std::cerr << "You do not have any task yet !!\n";
                continue;
            }
            return task_list;
        }
        auto tasks = menu_functions[num - 1]();
        task_list->splice(task_list->end(), *tasks);
    }
}
