//
// Created by xiamr on 6/14/19.
//

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30


//#include <boost/spirit/include/qi.hpp>
//#include <boost/spirit/include/phoenix.hpp>
//#include <boost/bind.hpp>
//#include <boost/lambda/lambda.hpp>
#include <boost/phoenix.hpp>
//#include <boost/variant.hpp>
//#include <boost/optional.hpp>
//#include <boost/fusion/sequence/intrinsic/at_c.hpp>
//#include <boost/fusion/include/at_c.hpp>
//#include <boost/phoenix/function/adapt_function.hpp>

// Boost metaprogramming library
#include <boost/mpl/vector.hpp>
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


using namespace std;

template<typename T1, typename T2>
struct add_item {
    add_item(T1 &v1, T2 &v2) : v1(v1), v2(v2) {};

    template<typename T>
    void operator()(boost::type<T>) {
        v1.emplace_back(bind(make_shared<T>));
        v2.emplace_back((boost::format("(%d) %s") % v1.size() % T::title()).str());
    }

    T1 &v1;
    T2 &v2;
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
            EquatorialAngle
    >;

    BOOST_MPL_ASSERT((mpl::equal<mpl::unique<components, is_same<mpl::_1, mpl::_2> >::type, components>));

    std::vector<std::function<shared_ptr<BasicAnalysis>()>> task_vec;
    std::vector<string> item_menu;

    mpl::for_each<components, boost::type<mpl::_>>(add_item<
            std::vector<std::function<shared_ptr<BasicAnalysis>()>>,
            std::vector<string>
    >(task_vec, item_menu));

    auto menu1 = [&item_menu]() {
        std::cout << "Please select the desired operation (Trajectory Analysis)" << std::endl;
        std::cout << "(0) Start\n";
        for_each(item_menu.cbegin(), item_menu.cend(), std::cout << boost::phoenix::placeholders::_1 << '\n');
        return choose<int>(0, mpl::size<components>::value, "select :");
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
