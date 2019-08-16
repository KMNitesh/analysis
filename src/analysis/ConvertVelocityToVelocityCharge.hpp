//
// Created by xiamr on 8/16/19.
//

#ifndef TINKER_CONVERTVELOCITYTOVELOCITYCHARGE_HPP
#define TINKER_CONVERTVELOCITYTOVELOCITYCHARGE_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"
#include "trr_writer.hpp"

class Frame;

class ConvertVelocityToVelocityCharge : public BasicAnalysis {
public:
    ConvertVelocityToVelocityCharge();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static std::string title() { return "Convert velocities to velocities * charge in order to use velacc of gmx"; }

protected:
    std::unique_ptr<TRRWriter> writer;

    std::string trr_vq_outfilename;
};


#endif //TINKER_CONVERTVELOCITYTOVELOCITYCHARGE_HPP
