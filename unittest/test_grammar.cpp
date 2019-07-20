
// Aim to test the spirit parser of amber style selection rule

#define BOOST_TEST_MODULE "Amber Style Selection Grammer Test"

#include <string>
#include <memory>

#include <boost/test/included/unit_test.hpp>

#include "grammar.hpp"
#include "atom.hpp"

struct grammar_fixture {
    grammar_fixture() : atom(std::make_shared<Atom>()) {
        //BOOST_TEST_MESSAGE( "setup fixture" );
    }

    ~grammar_fixture() {
        //BOOST_TEST_MESSAGE( "teardown fixture" );
    }

    bool result(bool need_suceess = true) {
        auto it = input_string.begin();
        bool status = qi::phrase_parse(it, input_string.end(), grammar, qi::ascii::space, mask);
        bool ret = status && it == input_string.end() && Atom::is_match(atom, mask);
        ret = need_suceess == ret;
        if (!ret) {
            boost::apply_visitor(print(), mask);
            if (!(status and (it == input_string.end()))) {
                std::cout << "error-pos : " << std::endl;
                std::cout << input_string << std::endl;
                for (auto iter = input_string.begin(); iter != it; ++iter) std::cout << " ";
                std::cout << "^" << std::endl;

            }
        }
        return ret;
    }

    Grammar<std::string::iterator, qi::ascii::space_type> grammar;
    Atom::AtomIndenter mask;
    std::shared_ptr<Atom> atom;
    std::string input_string;
};

BOOST_FIXTURE_TEST_SUITE(test_grammar, grammar_fixture)

    BOOST_AUTO_TEST_CASE(test_residue_name1) {
        input_string = ":ASP";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection1");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name2) {
        input_string = ":A*P";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection2");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name3) {
        input_string = ":A?P";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection3");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name4) {
        input_string = ":A=P";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection4");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name5) {
        input_string = ":*SP";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection5");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name6) {
        input_string = ":?SP";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection6");
    }

    BOOST_AUTO_TEST_CASE(test_residue_name7) {
        input_string = ":=SP";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue name selection7");
    }


    BOOST_AUTO_TEST_CASE(test_residue_number1) {
        input_string = ":1";
        atom->residue_name = "ASP";
        atom->residue_num = 1;

        BOOST_CHECK_MESSAGE(result(), "residue number selection1");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number2) {
        input_string = ":1*";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection2");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number3) {
        input_string = ":1?";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection3");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number4) {
        input_string = ":1=";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection4");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number5) {
        input_string = ":*0";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection5");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number6) {
        input_string = ":?0";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection6");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number7) {
        input_string = ":=0";
        atom->residue_name = "ASP";
        atom->residue_num = 10;

        BOOST_CHECK_MESSAGE(result(), "residue number selection7");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number_range) {
        input_string = ":1-10";
        atom->residue_name = "ASP";
        atom->residue_num = 5;

        BOOST_CHECK_MESSAGE(result(), "residue number range selection");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number_range2) {
        input_string = ":1-10#2";
        atom->residue_name = "ASP";
        atom->residue_num = 5;

        BOOST_CHECK_MESSAGE(result(), "residue number range selection with step v2");
    }

    BOOST_AUTO_TEST_CASE(test_residue_number_range3) {
        input_string = ":1-10#3";
        atom->residue_name = "ASP";
        atom->residue_num = 5;

        BOOST_CHECK_MESSAGE(!result(), "residue number range selection step v3");
    }


    BOOST_AUTO_TEST_CASE(test_atom_name) {
        input_string = "@CA";
        atom->atom_name = "CA";

        BOOST_CHECK_MESSAGE(result(), "atom name selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_name2) {
        input_string = "@C";
        atom->atom_name = "CA";

        BOOST_CHECK_MESSAGE(result(false), "atom name selection2");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number) {
        input_string = "@10";
        atom->seq = 10;

        BOOST_CHECK_MESSAGE(result(), "atom number selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number2) {
        input_string = "@1,2,5,9";
        atom->seq = 5;

        BOOST_CHECK_MESSAGE(result(), "atom number selection2");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number_seq_add_range) {
        input_string = "@1,2-5#2,9";
        atom->seq = 5;

        BOOST_CHECK_MESSAGE(result(false), "atom number selection2");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number_range) {
        input_string = "@1-10";
        atom->seq = 5;

        BOOST_CHECK_MESSAGE(result(), "atom number range selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number_range_with_step) {
        input_string = "@1-10#2";
        atom->seq = 5;

        BOOST_CHECK_MESSAGE(result(), "atom number range selection with step");
    }

    BOOST_AUTO_TEST_CASE(test_atom_number_range_with_step2) {
        input_string = "@1-10#3";
        atom->seq = 5;

        BOOST_CHECK_MESSAGE(result(false), "atom number range selection with step2");
    }

    BOOST_AUTO_TEST_CASE(test_atom_type_name) {
        input_string = "@%CA";
        atom->type_name = "CA";

        BOOST_CHECK_MESSAGE(result(), "atom type name selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_type_number) {
        input_string = "@%10";
        atom->typ = 10;

        BOOST_CHECK_MESSAGE(result(), "atom type number selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_type_number_range) {
        input_string = "@%1-10";
        atom->typ = 5;

        BOOST_CHECK_MESSAGE(result(), "atom type number range selection");
    }

    BOOST_AUTO_TEST_CASE(test_atom_type_number_range_with_step) {
        input_string = "@%1-10#4";
        atom->typ = 5;

        BOOST_CHECK_MESSAGE(result(), "atom type number range selection with step");
    }

    BOOST_AUTO_TEST_CASE(test_atom_symbol_name) {
        input_string = "@/C";
        atom->atom_symbol = "C";

        BOOST_CHECK_MESSAGE(result(), "atom symbol name selection");
    }


    BOOST_AUTO_TEST_CASE(test_not) {
        input_string = "!@10";
        atom->seq = 6;

        BOOST_CHECK_MESSAGE(result(), "not selection");
    }

    BOOST_AUTO_TEST_CASE(test_and) {
        input_string = ":1 & @10";
        atom->residue_name = "ASP";
        atom->residue_num = 1;
        atom->seq = 10;

        BOOST_CHECK_MESSAGE(result(), "and selection");
    }


    BOOST_AUTO_TEST_CASE(test_or) {
        input_string = ":1 | @10";
        atom->residue_name = "ASP";
        atom->residue_num = 2;
        atom->seq = 10;

        BOOST_CHECK_MESSAGE(result(), "or selection");
    }

    BOOST_AUTO_TEST_CASE(test_and_plus_or) {
        input_string = ":1 & @C | @%CA";
        atom->residue_name = "ASP";
        atom->residue_num = 2;
        atom->atom_name = "CB";
        atom->type_name = "CA";

        BOOST_CHECK_MESSAGE(result(), "and plus or selection");
    }

    BOOST_AUTO_TEST_CASE(test_parentheses1) {
        input_string = ":1 & ( @C | @%CA )";
        atom->residue_name = "ASP";
        atom->residue_num = 1;
        atom->atom_name = "CB";
        atom->type_name = "CA";

        BOOST_CHECK_MESSAGE(result(), "parentheses selection1");
    }

    BOOST_AUTO_TEST_CASE(test_parentheses2) {
        input_string = ":1-10#1 & ( @C | @%CA )";
        atom->residue_name = "ASP";
        atom->residue_num = 9;
        atom->atom_name = "CB";
        atom->type_name = "CA";

        BOOST_CHECK_MESSAGE(result(), "parentheses selection2");
    }

    BOOST_AUTO_TEST_CASE(zero_step) {
        input_string = ":1-10#0 & ( @C | @%CA )";
        atom->residue_name = "ASP";
        atom->residue_num = 9;
        atom->atom_name = "CB";
        atom->type_name = "CA";

        BOOST_CHECK_MESSAGE(result(false), "zero step");
    }


BOOST_AUTO_TEST_SUITE_END()

