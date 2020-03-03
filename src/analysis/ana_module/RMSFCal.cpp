//
// Created by xiamr on 6/14/19.
//
#include "ana_module/RMSFCal.hpp"

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>

#include "data_structure/frame.hpp"
#include "nlohmann/json.hpp"
#include "utils/common.hpp"

RMSFCal::RMSFCal() {
    enable_outfile = true;
    enable_forcefield = true;
}

void RMSFCal::process(std::shared_ptr<Frame> &frame) {
    int nfit1 = static_cast<int>(atoms_for_superpose.size());
    int nfit2 = static_cast<int>(atoms_for_superpose_and_rmsfcalc.size());
    int nfit = nfit1 + nfit2;
    int n = static_cast<int>(atoms_for_rmsfcalc.size());

    if (steps == 0) {
        save_coord(x1.data(), y1.data(), z1.data(), frame);
        double mid[3];
        center(nfit + n, x1.data(), y1.data(), z1.data(), mid, nfit);
        update_accumulators(x1.data(), y1.data(), z1.data(), nfit1, nfit2, n);
        if (pdb_ostream) {
            for (auto &atom : atoms_for_first_frame_output) {
                const auto &coord = atom->getCoordinate();
                *pdb_ostream << boost::format("ATOM  %5s %-4s %3s  %4s    %8.3f%8.3f%8.3f  1.00  0.00            \n") %
                                    atom->seq % atom->atom_name % atom->residue_name.value() %
                                    atom->residue_num.value() % (std::get<0>(coord) - mid[0]) %
                                    (std::get<1>(coord) - mid[1]) % (std::get<2>(coord) - mid[2]);
            }
            *pdb_ostream << "TER\n";
        }
    } else {
        save_coord(x2.data(), y2.data(), z2.data(), frame);
        double mid[3];
        center(nfit + n, x2.data(), y2.data(), z2.data(), mid, nfit);
        quatfit(nfit + n, x1.data(), y1.data(), z1.data(), nfit + n, x2.data(), y2.data(), z2.data(), nfit);
        append_pdb();
        update_accumulators(x2.data(), y2.data(), z2.data(), nfit1, nfit2, n);
    }
    steps++;
}

void RMSFCal::update_accumulators(double x[], double y[], double z[], int nfit1, int nfit2, int n) {
    for (int index = 0; index < nfit2 + n; index++) {
        acc[index][0](x[index + nfit1]);
        acc[index][1](y[index + nfit1]);
        acc[index][2](z[index + nfit1]);
    }
}

void RMSFCal::save_coord(double *x, double *y, double *z, const std::shared_ptr<Frame> &frame) {
    PBCUtils::move(mols, frame);
    for (const auto &element :
         join(atoms_for_superpose, atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) | boost::adaptors::indexed()) {
        std::tie(x[element.index()], y[element.index()], z[element.index()]) = element.value()->getCoordinate();
    }
}

void RMSFCal::append_pdb() {
    if (!pdb_ostream) return;
    for (const auto &element :
         join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) | boost::adaptors::indexed()) {
        *pdb_ostream << boost::format("ATOM  %5s %-4s %3s  %4s    %8.3f%8.3f%8.3f  1.00  0.00            \n") %
                            element.value()->seq % element.value()->atom_name % element.value()->residue_name.value() %
                            element.value()->residue_num.value() % x2[element.index()] % y2[element.index()] %
                            z2[element.index()];
    }
    *pdb_ostream << "TER\n";
}

void RMSFCal::print(std::ostream &os) {
    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# mask_for_superpose  > " << mask_for_superpose << '\n';
    os << "# mask_for_rmsfcalc   > " << mask_for_rmsfcalc << '\n';
    if (pdb_ostream) os << "# mask_for_first_frame_output   > " << mask_for_first_frame_output << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("#%15s %15s\n") % "AtomSeq" % "RMSF(Ang)";
    std::map<std::shared_ptr<Atom>, double> rms_values;
    for (const auto &element :
         join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) | boost::adaptors::indexed()) {
        auto rms = rmsvalue(element.index());
        if (output_residue_average) rms_values[element.value()] = rms;
        os << boost::format(" %15d %15.8f\n") % element.value()->seq % rms;
    }
    os << std::string(50, '#') << '\n';

    if (output_residue_average) {
        os << "# Average for residues \n";
        os << boost::format("#%15s %15s %15s %15s\n") % "Residue Sequence" % "Resid" % "Resname" % "RMSF(Ang)";
        for (const auto &element : residues | boost::adaptors::indexed(1)) {
            const auto &residue = element.value();
            double rms = boost::accumulate(residue.atoms, 0.0,
                                           [&rms_values](auto sum, auto &atom) { return sum + rms_values[atom]; }) /
                         residue.atoms.size();
            os << boost::format(" %15d %15d %15s %15.8f\n") % element.index() % residue.num % residue.name % rms;
        }

        os << std::string(50, '#') << '\n';
    }
    os << ">>>JSON<<<\n";
    saveJson(os);
    os << "<<<JSON>>>\n";
}

void RMSFCal::saveJson(std::ostream &os) const {
    nlohmann::json json;
    json["title"] = title();
    json["mask_fro_superpose"] = to_string(mask_for_superpose);
    json["mask_for_rmsfcalc"] = to_string(mask_for_rmsfcalc);
    if (pdb_ostream) json["mask_for_first_frame_output"] = to_string(mask_for_first_frame_output);
    for (const auto &element :
         join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) | boost::adaptors::indexed()) {
        json["RMSF"] = {
            {"X", "AtomSeq"},
            {"Y", {{"name", "RMSF"}, {"unit", "Ang"}}},
        };
        json["RMSF"]["values"].push_back({element.value()->seq, rmsvalue(element.index())});
    }
    os << json;
}

void RMSFCal::readInfo() {
    Atom::select1group(mask_for_superpose, "Please enter atoms for superpose > ");
    Atom::select1group(mask_for_rmsfcalc, "Please enter atoms for rmsf calc > ");

    output_residue_average = choose_bool("average for each residue [N] > ", Default(false));

    if (choose_bool("Output superposed strcuture [N] > ", Default(false))) {
        Atom::select1group(mask_for_first_frame_output, "Please enter atoms for first frame output > ");
        std::string filename = choose_file("Enter pdb filename for output > ").extension("pdb");
        pdb_ostream = std::make_unique<std::ofstream>(filename);
        if (!(*pdb_ostream)) {
            std::cerr << "ERROR !! file \'" << filename << "\' can not open for writting\n";
            exit(EXIT_FAILURE);
        }
    }
}

double RMSFCal::rmsvalue(int index) const {
    using boost::accumulators::variance;
    return std::sqrt(variance(acc[index][0]) + variance(acc[index][1]) + variance(acc[index][2]));
}

void RMSFCal::center(int n_for_center, double *x, double *y, double *z, double mid[], int nfit) {
    mid[0] = mid[1] = mid[2] = 0.0;
    double norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        mid[0] += x[i];
        mid[1] += y[i];
        mid[2] += z[i];
        norm += 1.0;
    }
    mid[0] /= norm;
    mid[1] /= norm;
    mid[2] /= norm;

    for (int i = 0; i < n_for_center; ++i) {
        x[i] -= mid[0];
        y[i] -= mid[1];
        z[i] -= mid[2];
    }
}

void RMSFCal::add_atom_to_residue_vector(std::shared_ptr<Atom> &atom) {
    auto num = atom->residue_num.value();
    auto name = atom->residue_name.value();
    if (residues.empty()) {
        residues.emplace_back(num, std::move(name), atom->molecule.lock());
        residues.back().atoms.insert(atom);
    } else {
        auto &last_residue = residues.back();
        if (last_residue.num == num and last_residue.name == name and last_residue.mol == atom->molecule.lock()) {
            last_residue.atoms.insert(atom);
        } else {
            residues.emplace_back(num, std::move(name), atom->molecule.lock());
            residues.back().atoms.insert(atom);
        }
    }
}

void RMSFCal::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(frame->atom_list, [this](std::shared_ptr<Atom> &atom) {
        auto b_for_superpose = Atom::is_match(atom, mask_for_superpose);
        auto b_for_rmsf_calc = Atom::is_match(atom, mask_for_rmsfcalc);

        if (b_for_superpose and b_for_rmsf_calc)
            atoms_for_superpose_and_rmsfcalc.push_back(atom);
        else if (b_for_rmsf_calc)
            atoms_for_rmsfcalc.push_back(atom);
        else if (b_for_superpose)
            atoms_for_superpose.push_back(atom);
        if (output_residue_average and b_for_rmsf_calc) {
            add_atom_to_residue_vector(atom);
        }

        if (pdb_ostream and Atom::is_match(atom, mask_for_first_frame_output))
            atoms_for_first_frame_output.push_back(atom);
    });

    allocate_array_memory();

    mols = PBCUtils::calculate_intermol(join(atoms_for_superpose, atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc),
                                        frame);
}

void RMSFCal::allocate_array_memory() {
    auto nfit1 = atoms_for_superpose.size();
    auto nfit2 = atoms_for_superpose_and_rmsfcalc.size();
    auto n = atoms_for_rmsfcalc.size();

    auto total_size = nfit1 + nfit2 + n;

    x1.resize(total_size);
    y1.resize(total_size);
    z1.resize(total_size);

    x2.resize(total_size);
    y2.resize(total_size);
    z2.resize(total_size);

    acc.resize(total_size);
}

