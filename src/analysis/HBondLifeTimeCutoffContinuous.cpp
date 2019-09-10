//
// Created by xiamr on 9/10/19.
//

#include "std.hpp"
#include "HBondLifeTimeCutoffContinuous.hpp"

void HBondLifeTimeCutoffContinuous::print(std::ostream &os) {
    auto acf = calculateAcf();
    printData(os, acf, title());
}

std::vector<double> HBondLifeTimeCutoffContinuous::calculateAcf() const {
    auto max_time_grap_frame = ceil(max_time_grap_ps / time_increment_ps);
    std::vector<long> acf(1, 0);
    std::vector<long> ntime(1, 0);
    acf.reserve(max_time_grap_frame + 1);
    ntime.reserve(max_time_grap_frame + 1);

    for (auto &sample : hb_histroy) {
        for (std::size_t i = 0; i < sample->size() - 1; ++i) {
            for (std::size_t j = i; j < std::min<std::size_t>(sample->size(), max_time_grap_frame + 1 + i); ++j) {
                auto n = j - i;
                if (ntime.size() <= n) {
                    ntime.push_back(0);
                    acf.push_back(0);
                }
                ++ntime[n];
                if ((*sample)[i] != 0 and (*sample)[i] == (*sample)[j]) {
                    bool bFailed = false;
                    //Brute-force algorithm
                    //TODO: Optimize this loop
                    for (std::size_t k = i + 1; k < j; ++k) {
                        if ((*sample)[i] != (*sample)[k]) {
                            bFailed = true;
                            break;
                        }
                    }
                    if (!bFailed) {
                        ++acf[n];
                    }
                    ++acf[n];
                }
            }
        }
    }

    std::vector<double> acff(acf.size(), 0);
    acff[0] = double(acf[0]) / ntime[0];
    for (std::size_t i = 1; i < acf.size(); ++i) {
        acff[i] = acf[i] / (ntime[i] * acff[0]);
    }
    acff[0] = 1.0;
    return acff;
}
