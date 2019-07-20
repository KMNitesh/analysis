//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_TASKMENU_HPP
#define TINKER_TASKMENU_HPP

#include <list>
#include <memory>

#include "BasicAnalysis.hpp"

std::shared_ptr<std::list<std::shared_ptr<BasicAnalysis>>> getTasks();

#endif //TINKER_TASKMENU_HPP
