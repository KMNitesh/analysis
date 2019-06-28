//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_TRAJCONV_HPP
#define TINKER_TRAJCONV_HPP

#include <memory>
#include <string>

#include "BasicAnalysis.hpp"
#include "gro_writer.hpp"
#include "trr_writer.hpp"
#include "xtc_writer.hpp"
#include "netcdf_writer.hpp"

class Frame;

class Trajconv : public BasicAnalysis {
    enum class PBCType {
        None,
        OneAtom,
        OneMol
    } pbc_type;
    int step = 0;
    std::string grofilename;
    std::string xtcfilename;
    std::string trrfilename;
    std::string mdcrdfilename;
    XTCWriter xtc;
    TRRWriter trr;
    NetCDFWriter mdcrd;

    int num;

    bool enable_xtc = true;
    bool enable_trr = true;
    bool enable_gro = true;
    bool enable_mdcrd = true;
public:
    void process(std::shared_ptr<Frame> &frame) override;

    void print() override;

    void readInfo() override;

    static const std::string title() {
        return "Gromacs XTC & TRR & GRO & NetCDF Output";
    }

    void fastConvertTo(std::string target);

    void doPBC(std::shared_ptr<Frame> &frame) const;
};

#endif //TINKER_TRAJCONV_HPP
