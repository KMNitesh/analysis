

#ifndef TINKER_TRAJCONV_HPP
#define TINKER_TRAJCONV_HPP

#include <memory>
#include <string>

#include <boost/optional.hpp>

#include "ana_module/AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/TrajectoryFormatWriter.hpp"
#include "utils/TrajectoryWriterFactoryImpl.hpp"

class PBCUtils;

class Frame;

class Trajconv : public AbstractAnalysis {
public:
    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void CleanUp();

    void readInfo() override;

    static const std::string title() { return "Gromacs XTC & TRR & GRO & NetCDF & ARC Output"; }

    explicit Trajconv(
        std::shared_ptr<TrajectoryWriterFactoryInterface> factory = std::make_shared<TrajectoryWriterFactoryImpl>());

    void fastConvertTo(std::string target) noexcept(false);

    enum class PBCType { None, OneAtom, OneMol, AtomGroup, AllIntoBox };

    const auto &getWriters() const { return writers; }

    void setParameters(const std::string &out, const std::string &pbc, const AmberMask &pbcmask,
                       const AmberMask &outmask);

private:
    PBCType pbc_type;

    int step = 0;
    AmberMask mask;

    std::shared_ptr<TrajectoryWriterFactoryInterface> factory;
    std::vector<std::pair<std::string, std::shared_ptr<TrajectoryFormatWriter>>> writers;

    std::shared_ptr<PBCUtils> pbc_utils;

    AmberMask mask_for_writetraj;
    std::vector<std::shared_ptr<Atom>> atoms_for_writetraj;

protected:
    void selectPBCMode();

    void inputOutputFiles(std::istream &in = std::cin, std::ostream &out = std::cout);
};

#endif  // TINKER_TRAJCONV_HPP
