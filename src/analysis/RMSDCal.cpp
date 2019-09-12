//
// Created by xiamr on 6/14/19.
//

#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/checked_delete.hpp>
#include "RMSDCal.hpp"
#include "frame.hpp"

void RMSDCal::process(std::shared_ptr<Frame> &frame) {
    rmsds.push_back(rmsvalue(frame));
}

void RMSDCal::print(std::ostream &os) {
    os << "***************************\n";
    os << "****** RMSD Calculator ****\n";
    os << "# AmberMask for superpose : " << mask_for_superpose << '\n';
    os << "# AmberMask for rms calc  : " << mask_for_rmscalc << '\n';
    os << "***************************\n";
    for (const auto &ele : rmsds | boost::adaptors::indexed(1)) {
        os << ele.index() << "     " << ele.value() << '\n';
    }
    os << "***************************\n";
}

void RMSDCal::readInfo() {
    Atom::select1group(mask_for_superpose, "Please enter atoms for superpose > ");
    Atom::select1group(mask_for_rmscalc, "Please enter atoms for rms calc  > ");
}


double RMSDCal::rmsvalue(std::shared_ptr<Frame> &frame) {

    auto nfit = static_cast<int>(this->atoms_for_superpose.size());
    auto n_rms_calc = nfit + static_cast<int>(this->atoms_for_rmscalc.size());

    if (first_frame) {
        first_frame = false;
        save_frame_coord(x1, y1, z1, frame);
        return 0.0;
    }

    save_frame_coord(x2, y2, z2, frame);

    double mid[3];
    center(nfit, x1, y1, z1, nfit, x2, y2, z2, mid, nfit);

    quatfit(n_rms_calc, x1, y1, z1, n_rms_calc, x2, y2, z2, nfit);

    return rmsfit(x1, y1, z1, x2, y2, z2, n_rms_calc);
}

void RMSDCal::save_frame_coord(double x[], double y[], double z[], const std::shared_ptr<Frame> &frame) const {
    int index = 0;
    bool first_atom = true;
    double first_x, first_y, first_z;
    for (auto &atom : atoms_for_superpose) {
        if (first_atom) {
            first_atom = false;
            first_x = x[index] = atom->x;
            first_y = y[index] = atom->y;
            first_z = z[index] = atom->z;
        } else {
            double xr = atom->x - first_x;
            double yr = atom->y - first_y;
            double zr = atom->z - first_z;
            frame->image(xr, yr, zr);
            x[index] = first_x + xr;
            y[index] = first_y + yr;
            z[index] = first_z + zr;
        }
        index++;
    }
    for (auto &atom : atoms_for_rmscalc) {
        double xr = atom->x - first_x;
        double yr = atom->y - first_y;
        double zr = atom->z - first_z;
        frame->image(xr, yr, zr);
        x[index] = first_x + xr;
        y[index] = first_y + yr;
        z[index] = first_z + zr;
        index++;
    }
}

void RMSDCal::center(int n1, double x1[], double y1[], double z1[],
                     int n2, double x2[], double y2[], double z2[],
                     double mid[], int nfit) {
    mid[0] = mid[1] = mid[2] = 0.0;
    double norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        mid[0] += x2[i];
        mid[1] += y2[i];
        mid[2] += z2[i];
        norm += 1.0;
    }
    mid[0] /= norm;
    mid[1] /= norm;
    mid[2] /= norm;
    for (int i = 0; i < n2; ++i) {
        x2[i] -= mid[0];
        y2[i] -= mid[1];
        z2[i] -= mid[2];
    }

    mid[0] = mid[1] = mid[2] = 0.0;
    norm = 0.0;
    for (int i = 0; i < nfit; ++i) {
        mid[0] += x1[i];
        mid[1] += y1[i];
        mid[2] += z1[i];
        norm += 1.0;
    }
    mid[0] /= norm;
    mid[1] /= norm;
    mid[2] /= norm;
    for (int i = 0; i < n1; ++i) {
        x1[i] -= mid[0];
        y1[i] -= mid[1];
        z1[i] -= mid[2];
    }
}

void RMSDCal::quatfit(int /* n1 */, double x1[], double y1[], double z1[],
                      int n_rms_calc, double x2[], double y2[], double z2[], int nfit) {
    int i;
    //    int i1, i2;
    //    double weigh;
    double xrot, yrot, zrot;
    double xxyx, xxyy, xxyz;
    double xyyx, xyyy, xyyz;
    double xzyx, xzyy, xzyz;
    double q[4], d[4];
    double rot[3][3];
    double c[4][4], v[4][4];

    xxyx = 0.0;
    xxyy = 0.0;
    xxyz = 0.0;
    xyyx = 0.0;
    xyyy = 0.0;
    xyyz = 0.0;
    xzyx = 0.0;
    xzyy = 0.0;
    xzyz = 0.0;

    for (i = 0; i < nfit; ++i) {
        xxyx += x1[i] * x2[i];
        xxyy += y1[i] * x2[i];
        xxyz += z1[i] * x2[i];
        xyyx += x1[i] * y2[i];
        xyyy += y1[i] * y2[i];
        xyyz += z1[i] * y2[i];
        xzyx += x1[i] * z2[i];
        xzyy += y1[i] * z2[i];
        xzyz += z1[i] * z2[i];
    }

    c[0][0] = xxyx + xyyy + xzyz;
    c[0][1] = xzyy - xyyz;
    c[1][1] = xxyx - xyyy - xzyz;
    c[0][2] = xxyz - xzyx;
    c[1][2] = xxyy + xyyx;
    c[2][2] = xyyy - xzyz - xxyx;
    c[0][3] = xyyx - xxyy;
    c[1][3] = xzyx + xxyz;
    c[2][3] = xyyz + xzyy;
    c[3][3] = xzyz - xxyx - xyyy;

    jacobi(4, c, d, v);

    q[0] = v[0][3];
    q[1] = v[1][3];
    q[2] = v[2][3];
    q[3] = v[3][3];

    rot[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    rot[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    rot[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);
    rot[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    rot[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
    rot[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);
    rot[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
    rot[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    rot[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];

    for (i = 0; i < n_rms_calc; ++i) {
        xrot = x2[i] * rot[0][0] + y2[i] * rot[0][1] + z2[i] * rot[0][2];
        yrot = x2[i] * rot[1][0] + y2[i] * rot[1][1] + z2[i] * rot[1][2];
        zrot = x2[i] * rot[2][0] + y2[i] * rot[2][1] + z2[i] * rot[2][2];
        x2[i] = xrot;
        y2[i] = yrot;
        z2[i] = zrot;
    }
}

double RMSDCal::rmsfit(double x1[], double y1[], double z1[],
                       double x2[], double y2[], double z2[], int n_rms_calc) {

    double fit;
    double xr, yr, zr, dist2;
    double norm;

    fit = 0.0;
    norm = 0.0;
    for (int i = 0; i < n_rms_calc; ++i) {
        xr = x1[i] - x2[i];
        yr = y1[i] - y2[i];
        zr = z1[i] - z2[i];
        dist2 = xr * xr + yr * yr + zr * zr;
        norm += 1.0;
        fit += dist2;
    }
    return sqrt(fit / norm);
}

void RMSDCal::jacobi(int n, double a[4][4], double d[], double v[4][4]) {
    // taken from tinker

    int i, j, k;
    int ip, iq;
    int nrot, maxrot;
    double sm, tresh, s, c, t;
    double theta, tau, h, g, p;
    //    double *b; //traditional point
    //    double *z;
    double b[4];
    double z[4];
    //    b = new double[n];
    //    z = new double[n];

    maxrot = 100;
    nrot = 0;
    for (ip = 0; ip < n; ++ip) {
        for (iq = 0; iq < n; ++iq) {
            v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    for (ip = 0; ip < n; ++ip) {
        b[ip] = a[ip][ip];
        d[ip] = b[ip];
        z[ip] = 0.0;
    }

    //  perform the jacobi rotations

    for (i = 0; i < maxrot; ++i) {
        sm = 0.0;
        for (ip = 0; ip < n - 1; ++ip) {
            for (iq = ip + 1; iq < n; ++iq) {
                sm += abs(a[ip][iq]);
            }
        }
        if (sm == 0.0) goto label_10;
        if (i < 3) {
            tresh = 0.2 * sm / (n * n);
        } else {
            tresh = 0.0;
        }
        for (ip = 0; ip < n - 1; ++ip) {
            for (iq = ip + 1; iq < n; ++iq) {
                g = 100.0 * abs(a[ip][iq]);
                if (i > 3 && abs(d[ip]) + g == abs(d[ip]) && abs(d[iq]) + g == abs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (abs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (abs(h) + g == abs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 0; j <= ip - 1; ++j) {
                        g = a[j][ip];
                        h = a[j][iq];
                        a[j][ip] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = ip + 1; j <= iq - 1; ++j) {
                        g = a[ip][j];
                        h = a[j][iq];
                        a[ip][j] = g - s * (h + g * tau);
                        a[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j < n; ++j) {
                        g = a[ip][j];
                        h = a[iq][j];
                        a[ip][j] = g - s * (h + g * tau);
                        a[iq][j] = h + s * (g - h * tau);
                    }
                    for (j = 0; j < n; ++j) {
                        g = v[j][ip];
                        h = v[j][iq];
                        v[j][ip] = g - s * (h + g * tau);
                        v[j][iq] = h + s * (g - h * tau);
                    }
                    ++nrot;
                }
            }
        }
        for (ip = 0; ip < n; ++ip) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    label_10:
    //    delete [] b; b = nullptr;
    //    delete [] z; z = nullptr;

    if (nrot == maxrot)
        std::cerr << " JACOBI  --  Matrix Diagonalization not Converged" << std::endl;

    for (i = 0; i < n - 1; ++i) {
        k = i;
        p = d[i];
        for (j = i + 1; j < n; ++j) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; ++j) {
                p = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }


}

void RMSDCal::processFirstFrame(std::shared_ptr<Frame> &frame) {
    boost::for_each(
            frame->atom_list,
            [this](std::shared_ptr<Atom> &atom) {
                if (Atom::is_match(atom, mask_for_superpose)) {
                    atoms_for_superpose.insert(atom);
                } else if (Atom::is_match(atom, mask_for_rmscalc)) {
                    atoms_for_rmscalc.insert(atom);
                }
            });
    auto n_size = atoms_for_superpose.size() + atoms_for_rmscalc.size();

    x1 = new double[n_size];
    y1 = new double[n_size];
    z1 = new double[n_size];


    x2 = new double[n_size];
    y2 = new double[n_size];
    z2 = new double[n_size];

}

RMSDCal::~RMSDCal() {
    boost::checked_array_delete(x1);
    boost::checked_array_delete(y1);
    boost::checked_array_delete(z1);

    boost::checked_array_delete(x2);
    boost::checked_array_delete(y2);
    boost::checked_array_delete(z2);
}
