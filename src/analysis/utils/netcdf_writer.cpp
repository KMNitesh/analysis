
#include "netcdf_writer.hpp"

#include <iostream>

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"

void NetCDFWriter::open(const std::string &filename) { this->filename = filename; }

void NetCDFWriter::close() {
    if (_is_open) {
        getNetcdfImpl()->netcdfClose(&NC);
        _is_open = false;
    }
}

void NetCDFWriter::write(const std::shared_ptr<Frame> &frame, const std::vector<std::shared_ptr<Atom>> &atoms) {
    if (!_is_open) {
        if (getNetcdfImpl()->netcdfCreate(&NC, filename.c_str(), atoms.size(), 1)) {
            throw std::runtime_error("Error open " + filename);
        }
        _is_open = true;
    }

    double x[NC.ncatom3];
    auto box = frame->box.getBoxParameter();
    double *p = x;
    for (auto &atom : atoms) {
        *p = atom->x;
        p++;
        *p = atom->y;
        p++;
        *p = atom->z;
        p++;
    }
    if (getNetcdfImpl()->netcdfWriteNextFrame(&NC, x, box.data())) {
        std::cerr << "Error write  mdcrd frame " << std::endl;
    }
    step++;
}
