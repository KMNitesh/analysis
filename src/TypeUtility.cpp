//
// Created by xiamr on 7/18/19.
//

// Boost metaprogramming library
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/string.hpp>

namespace mpl = boost::mpl;

#include "TypeUtility.hpp"

namespace TypeUtiltyInternal {
    template<typename T1>
    class AddItem {
    public:
        explicit AddItem(T1 &mapping) : mapping(mapping) {}

        template<typename T2>
        void operator()(boost::type<T2>) {
            mapping[typeid(T2)] = pretty_name<T2>();
        }

    private:
        T1 &mapping;
    };
}


std::string getPrettyName(const boost::any &v) {
    using components = mpl::vector<
            int,
            double,
            bool,
            Grid,
            Atom::Node,
            std::shared_ptr<VectorSelector>,
            std::shared_ptr<BasicAnalysis>
    >;

    BOOST_MPL_ASSERT((mpl::equal<mpl::unique<components, std::is_same<mpl::_1, mpl::_2> >::type, components>));

    std::unordered_map<std::type_index, std::string> mapping;

    mpl::for_each<components, boost::type<mpl::_>>(TypeUtiltyInternal::AddItem(mapping));

    auto it = mapping.find(v.type());
    return it != mapping.end() ? it->second : v.type().name();
}
