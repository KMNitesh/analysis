
#include "trajectoryreader.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <boost/spirit/include/qi.hpp>
#include <fstream>
#include <iostream>
#include <list>

#include "ReaderFactory.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "dsl/AmberMask.hpp"
#include "nlohmann/json.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"

bool TrajectoryReader::TrajectoryFile::is_in_range(uint pos) const {

    if (range.empty())
        return true;
    for (const auto &r : range) {
        const auto &first = boost::fusion::at_c<0>(r);
        const auto &second = boost::fusion::at_c<1>(r);
        if (second.has_value()) {
            if (first <= pos) {
                LessEqual less_equal_visitor(pos);
                if (boost::apply_visitor(less_equal_visitor, second.get()))
                    return true;
            }
        } else {
            if (first == pos)
                return true;
        }
    }
    return false;
}

bool TrajectoryReader::TrajectoryFile::is_end(uint pos) const {
    if (range.empty())
        return false;

    const auto last = range.back();
    const auto &first = boost::fusion::at_c<0>(last);
    const auto &second = boost::fusion::at_c<1>(last);
    if (second.has_value()) {
        LessEqual less_equal_visitor(pos);
        return !boost::apply_visitor(less_equal_visitor, second.get());
    } else {
        return pos > first;
    }
}

TrajectoryReader::TrajectoryFile::Ranges
TrajectoryReader::TrajectoryFile::parse_range(const std::string &range_string) {
    using namespace boost::spirit;
    Ranges ranges;
    if (auto it = std::begin(range_string);
        qi::phrase_parse(it, std::end(range_string),
                         ((qi::uint_ >> -('-' >> (qi::uint_ | qi::char_('$')))) % ',') >> qi::eoi, ascii::space,
                         ranges) and
        it == std::end(range_string)) {
        return ranges;
    } else {
        throw std::runtime_error("range syntax error <" + range_string + ">");
    }
}

void TrajectoryReader::add_trajectoy_file(const std::string &filename) {
    if (getFileType(filename) == FileType::JSON) {
        std::ifstream i(filename);
        nlohmann::json j;
        i >> j;
        for (auto &item : j) {
            std::string name = item.at("name");
            std::string mask = item.value("mask", "");
            traj_filenames.emplace(std::move(name),
                                   mask.empty() ? boost::blank{} : AmberMaskAST::parse_atoms(mask, true));
            try {
                std::string range_string = item.at("range");
                traj_filenames.back().range = TrajectoryFile::parse_range(range_string);
            } catch (...) {
                uint start = item.value("start", 1);
                uint end = item.value("end", 0);
                boost::variant<uint, char> end_tag;
                if (end == 0)
                    end_tag = '$';
                else
                    end_tag = end;
                traj_filenames.back().range.emplace_back(start, end_tag);
            }
        }
    } else
        traj_filenames.push(filename);
}

void TrajectoryReader::set_topology(const std::string &filename) { topology_filename = filename; }

void TrajectoryReader::set_mask(std::string mask_string) {
    if (!mask_string.empty()) {
        mask = AmberMaskAST::parse_atoms(mask_string, true);
    }
}
std::shared_ptr<Frame> TrajectoryReader::readOneFrame() {
    if (!frame) {
        readTopology();
        atoms_for_readtraj = isBlank(mask) ? frame->atom_list : PBCUtils::find_atoms(mask, frame);
    }
    for (;;) {
        if (!traj_reader) {
            if (traj_filenames.empty())
                return {};
            current_trajectory_file = std::move(traj_filenames.front());
            traj_filenames.pop();
            current_frame_pos = 0;
            traj_reader = ReaderFactory::getTrajectory(current_trajectory_file);
            traj_reader->open(current_trajectory_file);
            if (isBlank(mask)) {
                atoms_for_readtraj = isBlank(current_trajectory_file.mask)
                                         ? frame->atom_list
                                         : PBCUtils::find_atoms(current_trajectory_file.mask, frame);
            }
        }
        while (traj_reader->readOneFrame(frame, atoms_for_readtraj)) {
            ++current_frame_pos;
            if (current_trajectory_file.is_in_range(current_frame_pos)) {
                return frame;
            } else if (current_trajectory_file.is_end(current_frame_pos)) {
                break;
            }
        }
        traj_reader->close();
        traj_reader.reset();
    }
}

std::shared_ptr<Frame> TrajectoryReader::readTopology() {
    std::string file;
    if (topology_filename)
        file = topology_filename.get();
    else {
        if (!traj_filenames.empty() and getFileType(traj_filenames.front()) == FileType::ARC) {
            file = traj_filenames.front();
        } else {
            std::cerr << "Topology file not set !\n";
            std::exit(EXIT_FAILURE);
        }
    }

    auto reader = ReaderFactory::getTopology(file);
    frame = reader->read(file);
    frame->enable_bound = true; // TODO: detect PBC condtion
    for (const auto &element : frame->molecule_list | boost::adaptors::indexed(1)) {
        element.value()->sequence = element.index();
    }
    if (frame->atom_list.front()->residue_num.has_value()) {
        auto resideu_name = frame->atom_list.front()->residue_name.get();
        auto residue_num = frame->atom_list.front()->residue_num.get();
        auto molecule = frame->atom_list.front()->molecule.lock();

        uint current_real_residue_number = 1;
        for (auto &atom : frame->atom_list) {
            if (!(resideu_name == atom->residue_name.get() and residue_num == atom->residue_num.get() and
                  molecule == atom->molecule.lock())) {
                ++current_real_residue_number;
                resideu_name = atom->residue_name.get();
                residue_num = atom->residue_num.get();
                molecule = atom->molecule.lock();
            }
            atom->real_residue_number = current_real_residue_number;
        }
    }
    return frame;
}
