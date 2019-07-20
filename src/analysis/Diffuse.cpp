//
// Created by xiamr on 6/14/19.
//
#include <tbb/tbb.h>

#include "Diffuse.hpp"
#include "frame.hpp"

using namespace std;

void Diffuse::process(std::shared_ptr<Frame> &frame) {

    if (first_frame) {
        first_frame = false;
        for (auto &mol : frame->molecule_list) {
            mol->bExculde = true;
        }

        for (auto &atom : group) {
            auto mol = atom->molecule.lock();
            if (mol) {
                mol->bExculde = false;
            }
        }

        for (auto &mol : frame->molecule_list) {
            if (mol->bExculde) continue;
            mol->calc_mass();
            total_mol++;
        }
        cout << "Total molecule number : " << total_mol << endl;
        xcm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);
        ycm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);
        zcm = Eigen::MatrixXd::Zero(total_frame_number, total_mol);

        int mol_index = 0;
        for (auto &mol: frame->molecule_list) {
            if (!mol->bExculde) {
                auto coord = mol->calc_weigh_center(frame);
                xcm(0, mol_index) = get<0>(coord);
                ycm(0, mol_index) = get<1>(coord);
                zcm(0, mol_index) = get<2>(coord);
                mol_index++;
            }
        }
    } else {
        int mol_index = 0;
        for (auto &mol: frame->molecule_list) {
            if (!mol->bExculde) {
                auto coord = mol->calc_weigh_center(frame);
                double xold = xcm(steps - 1, mol_index);
                double yold = ycm(steps - 1, mol_index);
                double zold = zcm(steps - 1, mol_index);
                double xr = get<0>(coord) - xold;
                double yr = get<1>(coord) - yold;
                double zr = get<2>(coord) - zold;
                frame->image(xr, yr, zr);
                xcm(steps, mol_index) = xr + xold;
                ycm(steps, mol_index) = yr + yold;
                zcm(steps, mol_index) = zr + zold;
                mol_index++;
            }
        }
    }
    steps++;

}

void Diffuse::print(std::ostream &os) {

    vector<int> ntime(total_frame_number, 0);
    vector<double> xmsd(total_frame_number, 0.0);
    vector<double> ymsd(total_frame_number, 0.0);
    vector<double> zmsd(total_frame_number, 0.0);

    if (bSerial) {
        if (bTradition) {
            for (int i = 0; i < total_frame_number - 1; i++) {
                for (int j = i + 1; j < total_frame_number; j++) {
                    int m = j - i - 1;
                    ntime[m]++;
                    for (int k = 0; k < total_mol; k++) {
                        double xdiff = xcm(j, k) - xcm(i, k);
                        double ydiff = ycm(j, k) - ycm(i, k);
                        double zdiff = zcm(j, k) - zcm(i, k);
                        xmsd[m] += xdiff * xdiff;
                        ymsd[m] += ydiff * ydiff;
                        zmsd[m] += zdiff * zdiff;

                    }
                }
            }
        } else {
            for (int m = 0; m < total_frame_number - 1; m++) {
                for (int i = 0; i < total_frame_number - 1; i += m + 1) {
                    int j = i + m + 1;
                    if (j < total_frame_number) {
                        ntime[m]++;
                        for (int k = 0; k < total_mol; k++) {
                            double xdiff = xcm(j, k) - xcm(i, k);
                            double ydiff = ycm(j, k) - ycm(i, k);
                            double zdiff = zcm(j, k) - zcm(i, k);
                            xmsd[m] += xdiff * xdiff;
                            ymsd[m] += ydiff * ydiff;
                            zmsd[m] += zdiff * zdiff;
                        }
                    }
                }
            }
        }

    } else {

        class Body {
        public:
            int total_frame_number;
            vector<int> ntime;
            vector<double> xmsd;
            vector<double> ymsd;
            vector<double> zmsd;

            int total_mol;
            Eigen::MatrixXd &xcm;
            Eigen::MatrixXd &ycm;
            Eigen::MatrixXd &zcm;

            Body(int total_frame_number, int total_mol, Eigen::MatrixXd &xcm, Eigen::MatrixXd &ycm,
                 Eigen::MatrixXd &zcm) :
                    total_frame_number(total_frame_number),
                    ntime(total_frame_number - 1), xmsd(total_frame_number - 1),
                    ymsd(total_frame_number - 1), zmsd(total_frame_number - 1),
                    total_mol(total_mol), xcm(xcm), ycm(ycm), zcm(zcm) {}

            Body(const Body &body, tbb::split) :
                    total_frame_number(body.total_frame_number),
                    ntime(body.total_frame_number - 1), xmsd(body.total_frame_number - 1),
                    ymsd(body.total_frame_number - 1), zmsd(body.total_frame_number - 1),
                    total_mol(body.total_mol), xcm(body.xcm), ycm(body.ycm), zcm(body.zcm) {}

            void join(const Body &body) {
                for (int i = 0; i < total_frame_number - 1; i++) {
                    ntime[i] += body.ntime[i];
                    xmsd[i] += body.xmsd[i];
                    ymsd[i] += body.ymsd[i];
                    zmsd[i] += body.zmsd[i];
                }

            }

            void operator()(const tbb::blocked_range<int> &range) {
                for (int i = range.begin(); i != range.end(); i++) {
                    for (int j = i + 1; j < total_frame_number; j++) {
                        int m = j - i - 1;
                        ntime[m]++;
                        for (int k = 0; k < total_mol; k++) {
                            double xdiff = xcm(j, k) - xcm(i, k);
                            double ydiff = ycm(j, k) - ycm(i, k);
                            double zdiff = zcm(j, k) - zcm(i, k);
                            xmsd[m] += xdiff * xdiff;
                            ymsd[m] += ydiff * ydiff;
                            zmsd[m] += zdiff * zdiff;

                        }
                    }
                }
            }
        } body(total_frame_number, total_mol, xcm, ycm, zcm);

        tbb::parallel_reduce(tbb::blocked_range<int>(0, total_frame_number - 1), body, tbb::auto_partitioner());


        ntime = body.ntime;
        xmsd = body.xmsd;
        ymsd = body.ymsd;
        zmsd = body.zmsd;
    }

    const double dunits = 10.0;


    for (int i = 0; i < total_frame_number - 1; i++) {
        double counts = total_mol * ntime[i];
        xmsd[i] /= counts;
        ymsd[i] /= counts;
        zmsd[i] /= counts;
    }

    os << "*********************************************************\n";
    os << description();
    os << "Mean Squared Displacements and Self-Diffusion Constant\n";
    os << "    Time Gap      X MSD       Y MSD       Z MSD       R MSD        Diff Const\n";
    os << "      (ps)       (Ang^2)     (Ang^2)     (Ang^2)     (Ang^2)     (x 10^-5 cm**2/sec)\n";

    for (int i = 0; i < total_frame_number - 1; i++) {
        double delta = time_increment_ps * (i + 1);
        double xvalue = xmsd[i];
        double yvalue = ymsd[i];
        double zvalue = zmsd[i];
        double rvalue = xmsd[i] + ymsd[i] + zmsd[i];
        double dvalue = dunits * rvalue / delta / 6.0;
        os << boost::format("%12.2f%12.2f%12.2f%12.2f%12.2f%12.4f\n") %
              delta % xvalue % yvalue % zvalue % rvalue % dvalue;
    }
    os << "*********************************************************\n";

}

void Diffuse::setParameters(const Atom::Node &mask, double time_increment_ps, int total_frames,
                            const std::string &outfilename) {

    this->ids = mask;

    if (time_increment_ps <= 0) {
        throw runtime_error("`time_increment_ps' must great than zero");
    }
    this->time_increment_ps = time_increment_ps;


    if (total_frames <= 0) {
        throw runtime_error("`total_frames' must great than zero");
    }
    this->total_frame_number = total_frames;

    this->outfilename = outfilename;
    boost::trim(this->outfilename);
    if (this->outfilename.empty()) {
        throw runtime_error("outfilename cannot empty");
    }

    bSerial = false;
    enable_tbb = true;
}

void Diffuse::readInfo() {
    while (true) {
        time_increment_ps = 0.1;
        string input_line = input(" Enter the Time Increment in Picoseconds [0.1]: ");
        boost::trim(input_line);
        if (!input_line.empty()) {
            time_increment_ps = stod(input_line);
            if (time_increment_ps <= 0.0) {
                cout << "error time increment " << time_increment_ps << endl;
                continue;
            }
        }
        break;
    }
    while (true) {
        string input_line = input(" Enter the Total Frame Number: ");
        boost::trim(input_line);
        if (!input_line.empty()) {
            total_frame_number = stoi(input_line);
            if (total_frame_number <= 0) {
                cerr << "error total frame number " << total_frame_number << endl;
                continue;
            }
            break;
        }
    }
    Atom::select1group(ids, "select group :");

    while (true) {
        string input_line = input(" serial ? [T]:");
        boost::trim(input_line);
        if (!input_line.empty()) {
            if (input_line.compare("T") == 0) {
                bSerial = true;
                while (true) {
                    string input_line = input(" trandition ? [T]:");
                    boost::trim(input_line);
                    if (!input_line.empty()) {
                        if (input_line.compare("T") == 0) {
                            bTradition = true;
                            break;
                        } else if (input_line.compare("F") == 0) {
                            bTradition = false;
                            break;
                        }
                    }
                }
                break;
            } else if (input_line.compare("F") == 0) {
                bSerial = false;
                enable_tbb = true;
                break;
            }
        }
    }


}

void Diffuse::processFirstFrame(std::shared_ptr<Frame> &frame) {
    std::for_each(frame->atom_list.begin(), frame->atom_list.end(),
                  [this](shared_ptr<Atom> &atom) {
                      if (Atom::is_match(atom, this->ids)) this->group.insert(atom);
                  });
}

string Diffuse::description() {
    stringstream ss;
    string title_line = "------ " + title() + " ------";
    ss << title_line << "\n";
    ss << " mask              = [ " << ids << " ]\n";
    ss << " time_increment_ps = " << time_increment_ps << " (ps)\n";
    ss << " total_frames      = " << total_frame_number << " (frames)\n";
    ss << " outfilename       = " << outfilename << "\n";
    ss << string(title_line.size(), '-') << '\n';
    return ss.str();
}

