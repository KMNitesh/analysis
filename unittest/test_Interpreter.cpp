//
// Created by xiamr on 7/14/19.
//

#include <gmock/gmock.h>

#include "dsl/Interpreter.hpp"

using namespace std;
using namespace testing;


class InterpreterGrammarTest : public Test {
protected:
    void SetUp() override {

    }

    void pass() {
        it = input_string.begin();
        ASSERT_TRUE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
    }

    InterpreterGrammarT grammar;
    boost::any ast;
    string input_string;
    string::iterator it;
    Interpreter interpreter;
    boost::any result;
};

TEST_F(InterpreterGrammarTest, LogicalOperation) {

    input_string = "ret = (1 == 2.0 || 4 > 3.1 + 0.4 && ( 9 + 2 >= 10.0 && ! false));";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(bool));
    ASSERT_THAT(boost::any_cast<bool>(result), Eq(true));

    ASSERT_THAT(interpreter.getVariables().count("ret"), Eq(1));

    ASSERT_EQ(interpreter.getVariables().at("ret").type(), typeid(bool));
    ASSERT_THAT(boost::any_cast<bool>(interpreter.getVariables().at("ret")), Eq(true));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationIntAdd) {

    input_string = R"(a = 1; b = 2; c = a + b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(int));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(3));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationIntSubtract) {

    input_string = R"(a = 1; b = 2; c = a - b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(int));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(-1));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationIntMultiply) {

    input_string = R"(a = 3; b = 2; c = a * b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(int));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(6));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationIntDivide) {

    input_string = R"(a = 3; b = 2; c = a / b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(int));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationIntMod) {

    input_string = R"(a = 3; b = 2; c = a % b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(int));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
}


TEST_F(InterpreterGrammarTest, ArithmeticOperationDoubleAdd) {

    input_string = R"(a = 1.1; b = 2.4; c = a + b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(double));
    ASSERT_THAT(boost::any_cast<double>(result), DoubleEq(3.5));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationDoubleubtract) {

    input_string = R"(a = 1.1; b = 2.4; c = a - b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(double));
    ASSERT_THAT(boost::any_cast<double>(result), DoubleEq(-1.3));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationDoubleMultiply) {

    input_string = R"(a = 1.2; b = 2.4; c = a * b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(double));
    ASSERT_THAT(boost::any_cast<double>(result), DoubleEq(2.88));
}

TEST_F(InterpreterGrammarTest, ArithmeticOperationDoubleDivide) {

    input_string = R"(a = 1.1; b = 2.4; c = a / b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(double));
    ASSERT_THAT(boost::any_cast<double>(result), DoubleEq(1.1 / 2.4));
}


TEST_F(InterpreterGrammarTest, Bitwise_AND_Operation) {

    input_string = R"( a = [:1]; b = [@O]; c = a & b; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(Atom::Node));
    ASSERT_THAT(boost::any_cast<Atom::Node>(result), Eq(
            Atom::Node(make_shared<Atom::Operator>(
                    Atom::Op::AND,
                    boost::any_cast<Atom::Node>(interpreter.getVariables()["a"]),
                    boost::any_cast<Atom::Node>(interpreter.getVariables()["b"])
            ))));
}

TEST_F(InterpreterGrammarTest, Bitwise_OR_Operation) {
    InterpreterGrammarT grammar;

    input_string = R"( a = [:1]; b = [@O]; c = a | b; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));

    ASSERT_EQ(result.type(), typeid(Atom::Node));
    ASSERT_THAT(boost::any_cast<Atom::Node>(result), Eq(
            Atom::Node(make_shared<Atom::Operator>(
                    Atom::Op::OR,
                    boost::any_cast<Atom::Node>(interpreter.getVariables()["a"]),
                    boost::any_cast<Atom::Node>(interpreter.getVariables()["b"])
            ))));
}

TEST_F(InterpreterGrammarTest, Bitwise_NOT_Operation) {

    input_string = R"( a = [:1]; c = ! a; )";
    pass();
    ASSERT_NO_THROW(result = interpreter.execute(ast));

    ASSERT_EQ(result.type(), typeid(Atom::Node));
    ASSERT_THAT(boost::any_cast<Atom::Node>(result), Eq(
            Atom::Node(make_shared<Atom::Operator>(
                    Atom::Op::NOT,
                    boost::any_cast<Atom::Node>(interpreter.getVariables()["a"])
            ))));
}

TEST_F(InterpreterGrammarTest, IF_ELSE_Operation) {

    input_string = R"( a = 2; if ( a > 1 ) { b = true; } else { b = false; } )";
    pass();
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<bool>(interpreter.getVariables()["b"]), Eq(true));
}

TEST_F(InterpreterGrammarTest, WHILE_Operation) {

    input_string = R"( a = 2; while ( a  < 5 ) { a += 2; } )";
    pass();
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["a"]), Eq(6));
}

TEST_F(InterpreterGrammarTest, DO_WHILE_Operation) {

    input_string = R"( a = 2; do { a += 2; } while (a < 5 ); )";
    pass();
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["a"]), Eq(6));
}

TEST_F(InterpreterGrammarTest, DO_UNTIL_Operation) {

    input_string = R"( a = 2; do { a += 2; } until (a > 5 ); )";
    pass();
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["a"]), Eq(6));
}

TEST_F(InterpreterGrammarTest, FOR_Operation) {

    input_string = R"( for (i = 1; i < 4; ++i ) {} )";
    pass();
    Interpreter interpreter;
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(4));
}

TEST_F(InterpreterGrammarTest, PostIncrement_Operation) {

    input_string = R"( i = 1; i++; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(2));
}

TEST_F(InterpreterGrammarTest, PostDecremnt_Operation) {

    input_string = R"( i = 1; i--; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(0));
}

TEST_F(InterpreterGrammarTest, PreIncrement_Operation) {

    input_string = R"( i = 1; ++i; )";
    pass();
    ASSERT_NO_THROW((interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(2));
}

TEST_F(InterpreterGrammarTest, PreDecremnt_Operation) {

    input_string = R"( i = 1; --i; )";
    pass();
    ASSERT_NO_THROW(interpreter.execute(ast));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(0));
}


TEST_F(InterpreterGrammarTest, SimpleAssign) {

    input_string = R"( i = 1; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(1));
}

TEST_F(InterpreterGrammarTest, ContinuousAssign) {
    input_string = R"( i = j = k = 1; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["j"]), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["k"]), Eq(1));
}


TEST_F(InterpreterGrammarTest, CompoundAssignAdd) {

    input_string = R"( i = 1; i+=1; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(2));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(2));
}


TEST_F(InterpreterGrammarTest, CompoundAssignSubstract) {

    input_string = R"( i = 1; i-=1; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(0));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(0));
}

TEST_F(InterpreterGrammarTest, CompoundAssignMultiply) {

    input_string = R"( i = 1; i *= 2; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(2));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(2));
}

TEST_F(InterpreterGrammarTest, CompoundAssignDivide) {

    input_string = R"( i = 1; i /= 2.0; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<double>(result), DoubleEq(0.5));
    ASSERT_THAT(boost::any_cast<double>(interpreter.getVariables()["i"]), DoubleEq(0.5));
}

TEST_F(InterpreterGrammarTest, CompoundAssignMod) {

    input_string = R"( i = 3; i %= 2; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(1));
}

TEST_F(InterpreterGrammarTest, DoubleNot) {

    input_string = R"( i = true; j = !!i; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<bool>(result), Eq(true));
    ASSERT_THAT(boost::any_cast<bool>(interpreter.getVariables()["i"]), Eq(true));
}

TEST_F(InterpreterGrammarTest, DoublePrefixPlus) {

    input_string = R"( i = 1; j = ++++i; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(3));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(3));
}

TEST_F(InterpreterGrammarTest, DoublePrefixPlusWithSuffixPlus) {

    input_string = R"( i = 1; j = i+++++i; )";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(4));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(3));
}

TEST(SkipperGrammarTest, CommentLine) {
    string input_string = R"(
                        // Comment

                             )";
    auto it = input_string.begin();
    ASSERT_TRUE(qi::phrase_parse(it, input_string.end(), SkipperT(), qi::ascii::space) && it == input_string.end());
}


TEST_F(InterpreterGrammarTest, Comments) {
    input_string = R"(#!/bin/analysis
                       // Comment
                         # Comment
                        /*  Comment */  a = [ @1-1536#3];
                            i = 1; # Comment // Comment
                        /*
                             MultiLine Commment
                             MultiLine Commment
                         */
                             )";

    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<int>(result), Eq(1));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["i"]), Eq(1));
    ASSERT_THAT(boost::any_cast<Atom::Node>(interpreter.getVariables()["a"]), Eq(
            Atom::Node(make_shared<Atom::atom_name_nums>(
                    Atom::select_ranges{Atom::numItemType(1, make_pair<uint>(1536, 3))}))
    ));

}

TEST_F(InterpreterGrammarTest, StringPlus) {
    input_string = R"( a = "aa"; b = "bb"; a += b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<string>(result), StrEq("aabb"));
    ASSERT_THAT(boost::any_cast<string>(interpreter.getVariables()["b"]), Eq("bb"));
}

TEST_F(InterpreterGrammarTest, StringPlusInt) {
    input_string = R"( a = "aa"; b = 10; a += b;)";
    pass();
    ASSERT_NO_THROW((result = interpreter.execute(ast)));
    ASSERT_THAT(boost::any_cast<string>(result), StrEq("aa10"));
    ASSERT_THAT(boost::any_cast<int>(interpreter.getVariables()["b"]), Eq(10));
}

TEST_F(InterpreterGrammarTest, StringPlusDouble) {

    input_string = R"( a = "aa"; b = 10.1; a += b;)";
    pass();
    ASSERT_THROW(interpreter.execute(ast), boost::bad_any_cast);
}

TEST_F(InterpreterGrammarTest, StringPlusBool) {

    input_string = R"( a = "aa"; b = true; a += b;)";
    pass();
    ASSERT_THROW(interpreter.execute(ast), boost::bad_any_cast);
}

TEST_F(InterpreterGrammarTest, KeywordForReserve) {

    input_string = R"( for = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}

TEST_F(InterpreterGrammarTest, KeywordIfReserve) {

    input_string = R"( if = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}

TEST_F(InterpreterGrammarTest, KeywordElseReserve) {

    input_string = R"( else = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}

TEST_F(InterpreterGrammarTest, KeywordDoReserve) {

    input_string = R"( do = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}

TEST_F(InterpreterGrammarTest, KeywordWhileReserve) {

    input_string = R"( while = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}

TEST_F(InterpreterGrammarTest, KeywordUntilReserve) {

    input_string = R"( until = 1; )";
    it = input_string.begin();
    ASSERT_FALSE(qi::phrase_parse(it, input_string.end(), grammar, SkipperT(), ast) && it == input_string.end());
}








