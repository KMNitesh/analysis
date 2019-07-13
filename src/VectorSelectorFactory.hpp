//
// Created by xiamr on 7/5/19.
//

#ifndef TINKER_VECTORSELECTORFACTORY_HPP
#define TINKER_VECTORSELECTORFACTORY_HPP


#include "VectorSelector.hpp"
#include "VectorSelectorFactoryGrammar.hpp"

class VectorSelectorFactory {
public:
    static std::shared_ptr<VectorSelector> getVectorSelector();

    static std::shared_ptr<VectorSelector> getVectorSelectorByAST(vectorSelectorNode &ast);

private:
    VectorSelectorFactory() = default;
};


#endif //TINKER_VECTORSELECTORFACTORY_HPP
