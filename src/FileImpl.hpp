//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_FILEIMPL_HPP
#define TINKER_FILEIMPL_HPP

#include <fstream>
#include "FileInterface.hpp"

class FileImpl : public FileInterface, public std::ofstream {

public:
    void open(const std::string &__s, std::ios_base::openmode __mode) override {
        std::ofstream::open(__s, __mode);
    }

    void close() override {
        std::ofstream::close();
    }

};


#endif //TINKER_FILEIMPL_HPP
