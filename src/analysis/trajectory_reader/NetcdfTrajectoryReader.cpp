#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "NetcdfTrajectoryReader.hpp"

bool NetcdfTrajectoryReader::open(const std::string &file) {
    NC = std::make_unique<struct AmberNetcdf>();
    if (netcdfLoad(NC.get(), file.c_str())) {
        std::cerr << "error open NETCDF file: " << file << std::endl;
        return false;
    }
    return true;
}

bool NetcdfTrajectoryReader::readOneFrame(std::shared_ptr<Frame> &frame) {
    auto coord = std::make_unique<double[]>(NC->ncatom3);
    double box[6];
    if (!netcdfGetNextFrame(NC.get(), coord.get(), box)) return false;
    int i = 0;
    for (auto &atom : frame->atom_list) {
        atom->x = coord[3 * i];
        atom->y = coord[3 * i + 1];
        atom->z = coord[3 * i + 2];
        i++;
    }
    frame->a_axis = box[0];
    frame->b_axis = box[1];
    frame->c_axis = box[2];
    frame->alpha = box[3];
    frame->beta = box[4];
    frame->gamma = box[5];

    if (frame->enable_bound) {
        frame->a_axis_half = frame->a_axis / 2;
        frame->b_axis_half = frame->b_axis / 2;
        frame->c_axis_half = frame->c_axis / 2;
    }
    return true;
}

void NetcdfTrajectoryReader::close() {
    if (NC) {
        netcdfClose(NC.get());
        NC.reset();
    }
}

NetcdfTrajectoryReader::~NetcdfTrajectoryReader() {
    NetcdfTrajectoryReader::close();
}
