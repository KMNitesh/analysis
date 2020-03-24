#include "NetcdfTrajectoryReader.hpp"

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"

bool NetcdfTrajectoryReader::open(const std::string &file) {
    NC = std::make_unique<struct AmberNetcdf>();
    if (netcdfLoad(NC.get(), file.c_str())) {
        std::cerr << "error open NETCDF file: " << file << std::endl;
        return false;
    }
    return true;
}

bool NetcdfTrajectoryReader::readOneFrameImpl(std::shared_ptr<Frame> &frame,
                                              const std::vector<std::shared_ptr<Atom>> &atoms) {
    auto coord = std::make_unique<double[]>(NC->ncatom3);
    double box[6];
    if (!netcdfGetNextFrame(NC.get(), coord.get(), box))
        return false;

    static bool has_Warning_d = false;
    if (!has_Warning_d and NC->ncatom != static_cast<int>(atoms.size())) {
        std::cerr << boost::format("WARNING: topology has %d atoms, whereas trajectory has %d\n") % atoms.size() %
                         NC->ncatom;
        has_Warning_d = true;
    }

    for (auto i : boost::irange(NC->ncatom)) {
        auto &atom = atoms[i];
        atom->x = coord[3 * i];
        atom->y = coord[3 * i + 1];
        atom->z = coord[3 * i + 2];
    }
    frame->box = PBCBox(box[0], box[1], box[2], box[3], box[4], box[5]);
    return true;
}

void NetcdfTrajectoryReader::close() {
    if (NC) {
        netcdfClose(NC.get());
        NC.reset();
    }
}

NetcdfTrajectoryReader::~NetcdfTrajectoryReader() { NetcdfTrajectoryReader::close(); }
