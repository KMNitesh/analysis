//
// Created by xiamr on 7/5/19.
//

#include <iostream>
#include <functional>
#include <unordered_map>
#include "common.hpp"
#include "VectorSelectorFactory.hpp"

#include "NormalVectorSelector.hpp"
#include "TwoAtomVectorSelector.hpp"
#include "DipoleVectorSelector.hpp"

using namespace std;


std::shared_ptr<VectorSelector> VectorSelectorFactory::getVectorSelector() {

    std::cout << "Vector Selector Menu\n";
    std::cout << " 1. " << NormalVectorSelector::title() << '\n';
    std::cout << " 2. " << TwoAtomVectorSelector::title() << '\n';
    std::cout << " 3. " << DipoleVectorSelector::title() << '\n';

    const std::unordered_map<int, std::function<std::shared_ptr<VectorSelector>()>> mapping = {
            {1, std::bind(std::make_shared<NormalVectorSelector>)},
            {2, std::bind(std::make_shared<TwoAtomVectorSelector>)},
            {3, std::bind(std::make_shared<DipoleVectorSelector>)}
    };

    return mapping.at(choose(1, static_cast<int>(mapping.size()), " select > "))();

}

std::shared_ptr<VectorSelector> VectorSelectorFactory::getVectorSelectorByAST(vectorSelectorNode &ast) {
    struct Selector : boost::static_visitor<std::shared_ptr<VectorSelector>> {
        std::shared_ptr<VectorSelector> operator()(const NormalVectorSelectorNode &ast) const {
            auto selector = make_shared<NormalVectorSelector>();
            selector->readAST(ast);
            return selector;
        }

        std::shared_ptr<VectorSelector> operator()(const TwoAtomVectorSelectorNode &ast) const {
            auto selector = make_shared<TwoAtomVectorSelector>();
            selector->readAST(ast);
            return selector;
        }

        std::shared_ptr<VectorSelector> operator()(const DipoleVectorSelectorNode &ast) const {
            auto selector = make_shared<DipoleVectorSelector>();
            selector->readAST(ast);
            return selector;
        }
    };

    return boost::apply_visitor(Selector(), ast);
}
