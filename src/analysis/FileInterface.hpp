//
// Created by xiamr on 6/28/19.
//

#ifndef TINKER_FILEINTERFACE_HPP
#define TINKER_FILEINTERFACE_HPP

#include <string>
#include <ostream>

class FileInterface : public std::ostream {
public:
    virtual void open(const std::string &__s, std::ios_base::openmode __mode) = 0;

    virtual void close() = 0;
};


#endif //TINKER_FILEINTERFACE_HPP
