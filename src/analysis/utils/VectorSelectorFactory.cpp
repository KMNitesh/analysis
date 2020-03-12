//
// Created by xiamr on 7/5/19.
//

#include "VectorSelectorFactory.hpp"

#include <functional>
#include <iostream>
#include <unordered_map>

#include "NormalVectorSelector.hpp"
#include "TwoAtomVectorSelector.hpp"
#include "ana_module/DipoleVectorSelector.hpp"
#include "common.hpp"

using namespace std;

std::shared_ptr<VectorSelector> VectorSelectorFactory::getVectorSelector() {
    std::cout << "Vector Selector Menu\n";
    std::cout << " 1. " << NormalVectorSelector::title() << '\n';
    std::cout << " 2. " << TwoAtomVectorSelector::title() << '\n';
    std::cout << " 3. " << DipoleVectorSelector::title() << '\n';

    const std::unordered_map<int, std::function<std::shared_ptr<VectorSelector>()>> mapping = {
        {1, std::bind(std::make_shared<NormalVectorSelector>)},
        {2, std::bind(std::make_shared<TwoAtomVectorSelector>)},
        {3, std::bind(std::make_shared<DipoleVectorSelector>)}};

    return mapping.at(choose(1, static_cast<int>(mapping.size()), " select > "))();
}
