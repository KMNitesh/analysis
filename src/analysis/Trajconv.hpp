//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_TRAJCONV_HPP
#define TINKER_TRAJCONV_HPP

#include <memory>
#include <string>

#include "AbstractAnalysis.hpp"
#include "TrajectoryWriterFactoryImpl.hpp"
#include "TrajectoryFormatWriter.hpp"

class PBCUtils;

class Frame;

class Trajconv : public AbstractAnalysis {
public:

    void processFirstFrame(std::shared_ptr<Frame> &) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void CleanUp();

    void readInfo() override;

    static const std::string title() {
        return "Gromacs XTC & TRR & GRO & NetCDF Output";
    }

    explicit Trajconv(std::shared_ptr<TrajectoryWriterFactoryInterface> factory
    = std::make_shared<TrajectoryWriterFactoryImpl>());

    void fastConvertTo(std::string target) noexcept(false);

    enum class PBCType {
        None,
        OneAtom,
        OneMol
    };

    const auto &getWriters() const {
        return writers;
    }

private:

    PBCType pbc_type;

    int step = 0;
    int num;

    std::shared_ptr<TrajectoryWriterFactoryInterface> factory;
    std::vector<std::pair<std::string, std::shared_ptr<TrajectoryFormatWriter>>> writers;

    std::shared_ptr<PBCUtils> pbc_utils;
protected:
    void selectPBCMode();

    void inputOutputFiles(std::istream &in = std::cin, std::ostream &out = std::cout);

    void initPBC(PBCType pbc_mode, int num);


};

#endif //TINKER_TRAJCONV_HPP
