
#include "utils/common.hpp"
#include "ArcTrajectoryReader.hpp"
#include "NetcdfTrajectoryReader.hpp"
#include "TrrTrajectoryReader.hpp"
#include "XtcTrajectoryReader.hpp"

#include "Mol2Reader.hpp"
#include "PrmtopReader.hpp"
#include "TprReader.hpp"

#include "ReaderFactory.hpp"

std::shared_ptr<TopologyInterface> ReaderFactory::getTopology(const std::string &filename) {
    const std::map<FileType, std::function<std::shared_ptr<TopologyInterface>()>> mapping{
            {FileType::MOL2,   [] { return std::make_shared<Mol2Reader>(); }},
            {FileType::PRMTOP, [] { return std::make_shared<PrmtopReader>(); }},
            {FileType::TPR,    [] { return std::make_shared<TprReader>(); }},
            {FileType::ARC,    [] { return std::make_shared<ArcTrajectoryReader>(); }}
    };
    auto it = mapping.find(getFileType(filename));
    return it != std::end(mapping) ? it->second.operator()() : std::shared_ptr<TopologyInterface>{};
}

std::shared_ptr<TrajectoryInterface> ReaderFactory::getTrajectory(const std::string &filename) {
    const std::map<FileType, std::function<std::shared_ptr<TrajectoryInterface>()>> mapping{
            {FileType::NC,  [] { return std::make_shared<NetcdfTrajectoryReader>(); }},
            {FileType::TRR, [] { return std::make_shared<TrrTrajectoryReader>(); }},
            {FileType::XTC, [] { return std::make_shared<XtcTrajectoryReader>(); }},
            {FileType::ARC, [] { return std::make_shared<ArcTrajectoryReader>(); }}
    };
    auto it = mapping.find(getFileType(filename));
    return it != std::end(mapping) ? it->second.operator()() : std::shared_ptr<TrajectoryInterface>{};
}
