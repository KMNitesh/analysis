//
// Created by xiamr on 6/14/19.
//
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/join.hpp>
#include "RMSFCal.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"
#include "nlohmann/json.hpp"

RMSFCal::RMSFCal() {
    enable_outfile = true;
    enable_forcefield = true;

    coords.reserve(getDefaultVectorReserve());
}

void RMSFCal::process(std::shared_ptr<Frame> &frame) {

    int nfit1 = static_cast<int>(atoms_for_superpose.size());
    int nfit2 = static_cast<int>(atoms_for_superpose_and_rmsfcalc.size());
    int nfit = nfit1 + nfit2;
    int n = static_cast<int>(atoms_for_rmsfcalc.size());

    if (steps == 0) {
        auto first = save_coord(x1, y1, z1, frame);
        double mid[3];
        center(nfit + n, x1, y1, z1, mid, nfit);
        coords.emplace_back(append_coord(x1, y1, z1, nfit1, nfit2, n));
        if (pdb_ostream) {
            for (auto &atom : atoms_for_first_frame_output) {
                auto shift = atom->getCoordinate();
                shift -= first;
                frame->image(shift);

                *pdb_ostream << boost::format("ATOM  %5s %-4s %3s  %4s    %8.3f%8.3f%8.3f  1.00  0.00            \n")
                                % atom->seq
                                % atom->atom_name
                                % atom->residue_name.value()
                                % atom->residue_num.value()
                                % (std::get<0>(first) + std::get<0>(shift) - mid[0])
                                % (std::get<1>(first) + std::get<1>(shift) - mid[1])
                                % (std::get<2>(first) + std::get<2>(shift) - mid[2]);
            }
            *pdb_ostream << "TER\n";
        }
    } else {
        save_coord(x2, y2, z2, frame);
        double mid[3];
        center(nfit + n, x2, y2, z2, mid, nfit);
        quatfit(nfit + n, x1, y1, z1, nfit + n, x2, y2, z2, nfit);
        auto f_coord = append_coord(x2, y2, z2, nfit1, nfit2, n);
        append_pdb(f_coord);
        coords.emplace_back(std::move(f_coord));
    }
    steps++;
}

std::vector<std::tuple<double, double, double>>
RMSFCal::append_coord(double x[], double y[], double z[], int nfit1, int nfit2, int n) {
    std::vector<std::tuple<double, double, double>> f_coord;
    f_coord.reserve(nfit2 + n);
    for (int index = 0; index < nfit2 + n; index++) {
        x_avg[index] += x[index + nfit1];
        y_avg[index] += y[index + nfit1];
        z_avg[index] += z[index + nfit1];
        f_coord.emplace_back(x[index + nfit1], y[index + nfit1], z[index + nfit1]);
    }
    return f_coord;
}

std::tuple<double, double, double>
RMSFCal::save_coord(double *x, double *y, double *z, const std::shared_ptr<Frame> &frame) {
    for (const auto &element : join(atoms_for_superpose, atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) |
                               boost::adaptors::indexed()) {
        if (element.index() == 0) {
            std::tie(x[0], y[0], z[0]) = element.value()->getCoordinate();
        } else {
            auto shift = element.value()->getCoordinate() - std::make_tuple(x[0], y[0], z[0]);
            frame->image(shift);
            std::tie(x[element.index()], y[element.index()], z[element.index()]) =
                    std::make_tuple(x[0], y[0], z[0]) + shift;
        }
    }
    return {x[0], y[0], z[0]};
}

void RMSFCal::append_pdb(const std::vector<std::tuple<double, double, double>> &f_coord) {
    if (!pdb_ostream) return;
    for (const auto &element : join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) |
                               boost::adaptors::indexed()) {
        *pdb_ostream << boost::format("ATOM  %5s %-4s %3s  %4s    %8.3f%8.3f%8.3f  1.00  0.00            \n")
                        % element.value()->seq
                        % element.value()->atom_name
                        % element.value()->residue_name.value()
                        % element.value()->residue_num.value()
                        % std::get<0>(f_coord[element.index()])
                        % std::get<1>(f_coord[element.index()])
                        % std::get<2>(f_coord[element.index()]);
    }
    *pdb_ostream << "TER\n";
}

void RMSFCal::print(std::ostream &os) {

    calculate_average_structure();

    os << std::string(50, '#') << '\n';
    os << "# " << title() << " # \n";
    os << "# mask_for_superpose  > " << mask_for_superpose << '\n';
    os << "# mask_for_rmsfcalc   > " << mask_for_rmsfcalc << '\n';
    if (pdb_ostream) os << "# mask_for_first_frame_output   > " << mask_for_first_frame_output << '\n';
    os << std::string(50, '#') << '\n';
    os << boost::format("#%15s %15s\n") % "AtomSeq" % "RMSF(Ang)";
    for (const auto &element : join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) |
                               boost::adaptors::indexed()) {
        os << boost::format(" %15d %15.8f\n") % element.value()->seq % rmsvalue(element.index());
    }
    os << std::string(50, '#') << '\n';
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
    for (const auto &element : join(atoms_for_superpose_and_rmsfcalc, atoms_for_rmsfcalc) |
                               boost::adaptors::indexed()) {
        json["RMSF"] = {
                {"X", "AtomSeq"},
                {"Y", {{"name", "RMSF"}, {"unit", "Ang"}}},
        };
        json["RMSF"]["values"].push_back({element.value()->seq, rmsvalue(element.index())});
    }
    os << json;
}

void RMSFCal::calculate_average_structure() {
    for (std::size_t index = 0; index < atoms_for_superpose_and_rmsfcalc.size() + atoms_for_rmsfcalc.size(); index++) {
        x_avg[index] /= steps;
        y_avg[index] /= steps;
        z_avg[index] /= steps;
    }
}

void RMSFCal::readInfo() {
    Atom::select1group(mask_for_superpose, "Please enter atoms for superpose > ");
    Atom::select1group(mask_for_rmsfcalc, "Please enter atoms for rmsf calc > ");
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

    double dx2_y2_z2 = 0.0;
    for (int frame = 0; frame < steps; frame++) {
        auto dx = std::get<0>(coords[frame][index]) - x_avg[index];
        auto dy = std::get<1>(coords[frame][index]) - y_avg[index];
        auto dz = std::get<2>(coords[frame][index]) - z_avg[index];
        dx2_y2_z2 += dx * dx + dy * dy + dz * dz;
    }
    return sqrt(dx2_y2_z2 / steps);
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

void RMSFCal::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(
            frame->atom_list,
            [this](std::shared_ptr<Atom> &atom) {
                auto b_for_superpose = Atom::is_match(atom, mask_for_superpose);
                auto b_for_rmsf_calc = Atom::is_match(atom, mask_for_rmsfcalc);

                if (b_for_superpose and b_for_rmsf_calc) atoms_for_superpose_and_rmsfcalc.push_back(atom);
                else if (b_for_rmsf_calc) atoms_for_rmsfcalc.push_back(atom);
                else if (b_for_superpose) atoms_for_superpose.push_back(atom);

                if (pdb_ostream and Atom::is_match(atom, mask_for_first_frame_output))
                    atoms_for_first_frame_output.push_back(atom);
            });

    allocate_array_memory();

}

void RMSFCal::allocate_array_memory() {
    auto nfit1 = atoms_for_superpose.size();
    auto nfit2 = atoms_for_superpose_and_rmsfcalc.size();
    auto n = atoms_for_rmsfcalc.size();

    x_avg = new double[nfit2 + n];
    y_avg = new double[nfit2 + n];
    z_avg = new double[nfit2 + n];

    x1 = new double[nfit1 + nfit2 + n];
    y1 = new double[nfit1 + nfit2 + n];
    z1 = new double[nfit1 + nfit2 + n];

    x2 = new double[nfit1 + nfit2 + n];
    y2 = new double[nfit1 + nfit2 + n];
    z2 = new double[nfit1 + nfit2 + n];
}

RMSFCal::~RMSFCal() {
    boost::checked_array_delete(x_avg);
    boost::checked_array_delete(y_avg);
    boost::checked_array_delete(z_avg);

    boost::checked_array_delete(x1);
    boost::checked_array_delete(y1);
    boost::checked_array_delete(z1);

    boost::checked_array_delete(x2);
    boost::checked_array_delete(y2);
    boost::checked_array_delete(z2);
}
