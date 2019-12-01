//
// Created by xiamr on 6/28/19.
//

#include <unordered_map>
#include <functional>
#include "TrajectoryWriterFactoryImpl.hpp"
#include "utils/netcdf_writer.hpp"
#include "utils/trr_writer.hpp"
#include "utils/xtc_writer.hpp"
#include "utils/gro_writer.hpp"

using namespace std;

shared_ptr<TrajectoryFormatWriter> TrajectoryWriterFactoryImpl::make_instance(FileType t) {
    unordered_map<FileType, function<shared_ptr<TrajectoryFormatWriter>()>> mapping = {
            {FileType::NC,  bind(make_shared<NetCDFWriter>)},
            {FileType::XTC, bind(make_shared<XTCWriter>)},
            {FileType::TRR, bind(make_shared<TRRWriter>)},
            {FileType::GRO, bind(make_shared<GROWriter>)},
    };
    auto it = mapping.find(t);
    if (it != mapping.end()) { return it->second(); }
    throw std::runtime_error("File Type Error");
}
