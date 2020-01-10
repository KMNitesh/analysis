//
// Created by xiamr on 3/17/19.
//

#ifndef TINKER_TRAJECTORYREADER_HPP
#define TINKER_TRAJECTORYREADER_HPP

#include "config.h"
#include <fstream>
#include <string>
#include <memory>
#include <list>
#include <vector>


namespace gmx {

#include "gromacs/fileio/xtcio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/utility/smalloc.h"

}


#include "utils/AmberNetcdf.h"
#include "utils/amber_netcdf.h"
#include "utils/gmxtrr.h"

#include "utils/common.hpp"

class Frame;


class TrajectoryReader {
    std::fstream position_file;
    std::fstream velocity_file; // The velocity file

    bool isbin = false;
    bool isnetcdf = false;
    bool istrr = false;
    bool isxtc = false;

    bool enable_binaray_file = false;

    std::string topology_filename;

    enum class TOPOLOGY_TYPE {
        ARC, MOL2, TPR, PRMTOP
    } topology_type;

    struct AmberNetcdf NC;
    std::shared_ptr<Frame> frame;

    gmx::t_fileio *fio = nullptr;

    bool first_time = true;
    bool openvel = false;

    int natoms, step;
    gmx::rvec *x = nullptr;
    gmx::real prec, time;

    std::list<std::string> arc_filename_list; // the continuous trajectory files

    void close();

    void readOneFrameVel();

    int readOneFrameNetCDF();

    int readOneFrameTrr();

    int readOneFrameXtc();

    void open(const std::string &filename);

protected:

    std::shared_ptr<Frame> readOneFrameTraj();

    std::shared_ptr<Frame> readOneFrameTpr();

    std::shared_ptr<Frame> readOneFramePrmtop();

    std::shared_ptr<Frame> readOneFrameMol2();

    std::shared_ptr<Frame> readOneFrameArc();

    void assignAtom2Molecule();

public:

    void add_filename(const std::string &filename);

    void add_topology(const std::string &filename);

    std::string getTopologyFilename() { return topology_filename; }


    std::shared_ptr<Frame> readOneFrame();

    std::shared_ptr<Frame> readTopology();
};


#endif //TINKER_TRAJECTORYREADER_HPP
