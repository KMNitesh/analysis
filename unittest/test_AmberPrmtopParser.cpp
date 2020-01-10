
#include <fstream>
#include <gmock/gmock.h>
#include <boost/spirit/include/qi.hpp>
#include "trajectory_reader/PrmtopParser.hpp"


using namespace testing;
using namespace boost::spirit;

TEST(TEST_PrmtopGrammar, Line) {

    std::stringstream ss("%VERSION  VERSION_STAMP = V0001.000  DATE = 08/07/18  23:35:35\n");
    ss.unsetf(std::ios::skipws);

    using namespace qi;
    qi::rule<boost::spirit::istream_iterator, std::string(), decltype(boost::spirit::ascii::space - eol)>
            version = "%VERSION" >> lexeme[+(char_ - eol)] >> eol;
    std::string value;
    boost::spirit::istream_iterator it(ss), end;

    ASSERT_TRUE(phrase_parse(it, end, version, boost::spirit::ascii::space - eol, value) and it == end);
    ASSERT_THAT(value, Eq("VERSION_STAMP = V0001.000  DATE = 08/07/18  23:35:35"));
}

TEST(TEST_PrmtopGrammar, Simple) {

    std::ifstream ifstream("amber_test_topology1.prmtop");

    boost::optional<PrmtopStruct> prmtop = PrmtopParser::parse(ifstream);

    ASSERT_THAT(prmtop.has_value(), Eq(true));
    ASSERT_THAT(prmtop->version, Eq("VERSION_STAMP = V0001.000  DATE = 04/06/18  21:15:52"));
    ASSERT_THAT(prmtop->title, Eq("MOL"));
    ASSERT_THAT(prmtop->pointers, ElementsAre(53, 4, 34, 18, 74, 28, 124, 40, 0, 0,
                                              302, 1, 18, 28, 40, 4, 8, 7, 4, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 53, 0,
                                              0));
    ASSERT_THAT(prmtop->atom_name,
                ElementsAre("S1", "P1", "C1", "H1", "H2", "C2", "H3", "C3", "H4", "H5",
                            "C4", "C5", "H6", "H7", "H8", "C6", "H9", "H10", "H11", "C7",
                            "H12", "H13", "H14", "C8", "H15", "H16", "H17", "S2", "C9",
                            "H18", "H19", "C10", "H20", "C11", "H21", "H22", "C12", "C13",
                            "H23", "H24", "H25", "C14", "H26", "H27", "H28", "C15", "H29",
                            "H30", "H31", "C16", "H32", "H33", "H34"));
    ASSERT_THAT(prmtop->charge,
                Pointwise(DoubleEq(), {
                        -9.62541975E+00, 2.23292420E+00, -3.18292559E+00, 8.78296638E-01, 8.78296638E-01,
                        9.61906017E+00, -8.16577708E-01, -5.62693691E+00, -2.98590608E-01, -2.98590608E-01,
                        1.38517360E+01, -6.13030972E+00, 6.90497614E-01, 6.90497614E-01, 6.90497614E-01,
                        -5.84547695E+00, 8.05298104E-01, 8.05298104E-01, 8.05298104E-01, -6.13030972E+00,
                        6.90497614E-01, 6.90497614E-01, 6.90497614E-01, -6.13030972E+00, 6.90497614E-01,
                        6.90497614E-01, 6.90497614E-01, -9.62541975E+00, -3.18292559E+00, 8.78296638E-01,
                        8.78296638E-01, 9.61906017E+00, -8.16577708E-01, -5.62693691E+00, -2.98590608E-01,
                        -2.98590608E-01, 1.38517360E+01, -6.13030972E+00, 6.90497614E-01, 6.90497614E-01,
                        6.90497614E-01, -5.84547695E+00, 8.05298104E-01, 8.05298104E-01, 8.05298104E-01,
                        -6.13030972E+00, 6.90497614E-01, 6.90497614E-01, 6.90497614E-01, -6.13030972E+00,
                        6.90497614E-01, 6.90497614E-01, 6.90497614E-01
                }));

    ASSERT_THAT(prmtop->atomic_number, ElementsAre(
            16, 15, 6, 1, 1, 6, 1, 6, 1, 1,
            6, 6, 1, 1, 1, 6, 1, 1, 1, 6,
            1, 1, 1, 6, 1, 1, 1, 16, 6, 1,
            1, 6, 1, 6, 1, 1, 6, 6, 1, 1,
            1, 6, 1, 1, 1, 6, 1, 1, 1, 6,
            1, 1, 1
    ));

    ASSERT_THAT(prmtop->mass,
                Pointwise(DoubleEq(), {
                        3.20600000E+01, 3.09700000E+01, 1.20100000E+01, 1.00800000E+00, 1.00800000E+00,
                        1.20100000E+01, 1.00800000E+00, 1.20100000E+01, 1.00800000E+00, 1.00800000E+00,
                        1.20100000E+01, 1.20100000E+01, 1.00800000E+00, 1.00800000E+00, 1.00800000E+00,
                        1.20100000E+01, 1.00800000E+00, 1.00800000E+00, 1.00800000E+00, 1.20100000E+01,
                        1.00800000E+00, 1.00800000E+00, 1.00800000E+00, 1.20100000E+01, 1.00800000E+00,
                        1.00800000E+00, 1.00800000E+00, 3.20600000E+01, 1.20100000E+01, 1.00800000E+00,
                        1.00800000E+00, 1.20100000E+01, 1.00800000E+00, 1.20100000E+01, 1.00800000E+00,
                        1.00800000E+00, 1.20100000E+01, 1.20100000E+01, 1.00800000E+00, 1.00800000E+00,
                        1.00800000E+00, 1.20100000E+01, 1.00800000E+00, 1.00800000E+00, 1.00800000E+00,
                        1.20100000E+01, 1.00800000E+00, 1.00800000E+00, 1.00800000E+00, 1.20100000E+01,
                        1.00800000E+00, 1.00800000E+00, 1.00800000E+00
                }));

    ASSERT_THAT(prmtop->atom_type_index, ElementsAre(
            1, 2, 3, 4, 4, 3, 4, 3, 4, 4,
            3, 3, 4, 4, 4, 3, 4, 4, 4, 3,
            4, 4, 4, 3, 4, 4, 4, 1, 3, 4,
            4, 3, 4, 3, 4, 4, 3, 3, 4, 4,
            4, 3, 4, 4, 4, 3, 4, 4, 4, 3,
            4, 4, 4
    ));

    ASSERT_THAT(prmtop->number_excluded_atoms, ElementsAre(
            10, 15, 17, 7, 6, 14, 8, 19, 6, 5,
            13, 11, 4, 3, 2, 3, 2, 1, 1, 7,
            3, 2, 1, 3, 2, 1, 1, 4, 12, 5,
            4, 12, 8, 19, 6, 5, 13, 11, 4, 3,
            2, 3, 2, 1, 1, 7, 3, 2, 1, 3,
            2, 1, 1
    ));

    ASSERT_THAT(prmtop->nonbonded_parm_index, ElementsAre(
            1, 2, 4, 7, 2, 3, 5, 8, 4, 5,
            6, 9, 7, 8, 9, 10
    ));

    ASSERT_THAT(prmtop->residue_label, ElementsAre("MOL"));

    ASSERT_THAT(prmtop->residue_pointer, ElementsAre(1));

    ASSERT_THAT(prmtop->bond_force_constant,
                Pointwise(DoubleEq(), {2.43900000E+02, 2.43300000E+02, 3.30600000E+02, 3.00900000E+02}));

    ASSERT_THAT(prmtop->bond_equil_value,
                Pointwise(DoubleEq(), {1.93410000E+00, 1.83950000E+00, 1.09690000E+00, 1.53750000E+00}));

    ASSERT_THAT(prmtop->angle_force_constant,
                Pointwise(DoubleEq(), {
                        3.69800000E+01, 3.91300000E+01, 5.33100000E+01, 7.69600000E+01, 3.64600000E+01,
                        4.63400000E+01, 6.28600000E+01, 3.94000000E+01
                }));

    ASSERT_THAT(prmtop->angle_equil_value,
                Pointwise(DoubleEq(), {
                        1.99665752E+00, 1.99194513E+00, 1.89246132E+00, 1.95511867E+00, 1.85004980E+00,
                        1.91637234E+00, 1.94621748E+00, 1.87762601E+00
                }));

    ASSERT_THAT(prmtop->dihedral_force_constant,
                Pointwise(DoubleEq(), {
                        2.22222222E-02, 1.55555556E-01, 1.60000000E-01, 2.00000000E-01, 2.50000000E-01,
                        1.80000000E-01, 1.50000000E-01
                }));

    ASSERT_THAT(prmtop->dihedral_periodicity,
                Pointwise(DoubleEq(), {
                        3.00000000E+00, 3.00000000E+00, 3.00000000E+00, 1.00000000E+00, 2.00000000E+00,
                        3.00000000E+00, 3.00000000E+00
                }));

    ASSERT_THAT(prmtop->dihedral_phase,
                Pointwise(DoubleEq(), {
                        0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 3.14159400E+00, 3.14159400E+00,
                        0.00000000E+00, 0.00000000E+00
                }));
}