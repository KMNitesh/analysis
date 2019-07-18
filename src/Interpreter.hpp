#include <utility>


//
// Created by xiamr on 7/15/19.
//

#ifndef TINKER_INTERPRETER_HPP
#define TINKER_INTERPRETER_HPP

#ifndef BOOST_SPIRIT_DEBUG_MAIN_HPP
#define BOOST_SPIRIT_DEBUG_MAIN_HPP

#endif

#include <utility>
#include <unordered_map>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/phoenix.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/any.hpp>
#include <boost/type_index.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>

#include <boost/range/combine.hpp>
#include <functional>
#include "ThrowAssert.hpp"
#include "common.hpp"
#include "atom.hpp"
#include "grammar.hpp"
#include "TypeUtility.hpp"


namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;
namespace phoenix = boost::phoenix;

struct Identifer {
    explicit Identifer(std::string name, bool &_pass) : name(std::move(name)) {
        std::vector<std::string> keywords = {"for", "if", "else", "do", "while", "until"};
        if (boost::algorithm::one_of_equal(keywords, this->name)) {
            std::cerr << "keyword (" << this->name << ") is reserved and cannot use for indentifer name\n";
            _pass = false;
        }
    }

    Identifer() = default;

    std::string name;
};


struct Function {
    Function(boost::any name, boost::any arguments);

    boost::any name;
    std::vector<boost::any> arguments;
    std::vector<std::pair<boost::any, boost::any>> optional_arguments;
};

struct ForStmt {
    ForStmt(boost::any expr1, boost::any expr2, boost::any expr3, boost::any stmtBlock)
            : expr1(std::move(expr1)), expr2(std::move(expr2)), expr3(std::move(expr3)),
              stmt_block(std::move(stmtBlock)) {}

    boost::any expr1;
    boost::any expr2;
    boost::any expr3;
    boost::any stmt_block;

};


struct WhileStmt {
    WhileStmt(boost::any condition, boost::any stmtBlock)
            : condition(std::move(condition)), stmt_block(std::move(stmtBlock)) {}

    boost::any condition, stmt_block;
};

struct DoUntilStmt {
    DoUntilStmt(boost::any stmtBlock, boost::any condition)
            : stmt_block(std::move(stmtBlock)), condition(std::move(condition)) {}

    boost::any stmt_block, condition;
};

struct DoWhileStmt {
    DoWhileStmt(boost::any stmtBlock, boost::any condition)
            : stmt_block(std::move(stmtBlock)), condition(std::move(condition)) {}

    boost::any stmt_block, condition;
};


struct IfElseStmt {

    IfElseStmt(boost::any condtion, boost::any ifStmtBlock, boost::any elseStmtBlock)
            : condtion(std::move(condtion)), if_stmt_block(std::move(ifStmtBlock)),
              else_stmt_block(std::move(elseStmtBlock)) {}

    boost::any condtion;
    boost::any if_stmt_block, else_stmt_block;
};


enum class LogicalOp {
    Equal,
    NotEqual,
    And,
    Or,
    LessEqual,
    Less,
    GreatEqual,
    Great,
    Not
};

struct LogicalOperation {
    LogicalOperation(LogicalOp op, boost::any operand1, boost::any operand2 = {}) :
            op(op), operand1(std::move(operand1)), operand2(std::move(operand2)) {}

    LogicalOp op;

    boost::any operand1, operand2;
};


enum class ArithmeticOp {
    Plus,
    Minus,
    Subtract,
    Multiply,
    Divide,
    Mod,
    PostIncrement,
    PostDecrement,
    PreIncrement,
    PreDecrement
};

struct ArithmeticOperation {
    ArithmeticOperation(ArithmeticOp op, boost::any operand1, boost::any operand2 = {}) :
            op(op), operand1(std::move(operand1)), operand2(std::move(operand2)) {}

    ArithmeticOp op;
    boost::any operand1, operand2;
};

enum class BitwiseOp {
    And,
    Or
};

struct BitwiseOperation {
    BitwiseOperation(BitwiseOp op, boost::any operand1, boost::any operand2)
            : op(op), operand1(std::move(operand1)), operand2(std::move(operand2)) {}

    BitwiseOp op;
    boost::any operand1, operand2;
};

struct AssignStmt {
    enum class TYPE {
        Simple, Compound
    };

    AssignStmt(boost::any var, boost::any rhs, TYPE t = TYPE::Simple)
            : var(std::move(var)), rhs(std::move(rhs)), t(t) {}

    boost::any var, rhs;
    TYPE t;
};

/*
 *   This the SDT for subset of C and Python
 */

template<typename Iterator, typename Skipper>
struct InterpreterGrammar : qi::grammar<Iterator, boost::any(), Skipper> {

    Grammar<Iterator, qi::ascii::space_type> maskParser;

    qi::rule<Iterator, Atom::Node(), Skipper> mask;
    qi::rule<Iterator, boost::any(), Skipper> languague, stmt, quoted_string;
    qi::rule<Iterator, boost::any(), Skipper> assign_expr;
    qi::rule<Iterator, boost::any(), qi::locals<std::vector<boost::any>>, Skipper> stmts;
    qi::rule<Iterator, boost::any(), Skipper> parenthese_expr, function_call_expr;
    qi::rule<Iterator, boost::any(), Skipper> logical_or_expr, logical_and_expr,
            logical_equal_not_equal_expr, logical_comp_expr,
            arithmetic_low_expr, arithmetic_high_expr, prefix_unary_op_expr,
            literal, identifer, bitwise_or_expr, bitwise_and_expr, suffix_unary_op_expr;

    qi::rule<Iterator, boost::any(), qi::locals<std::vector<boost::any>>, Skipper> function_call_arguments_expr;

    qi::rule<Iterator, boost::any(), qi::locals<boost::any, boost::any, boost::any, boost::any>, Skipper> for_stmt;

    qi::rule<Iterator, boost::any(), qi::locals<boost::any, boost::any>, Skipper> while_stmt;
    qi::rule<Iterator, boost::any(), qi::locals<boost::any, boost::any>, Skipper> do_while_until_stmt;

    qi::rule<Iterator, boost::any(), qi::locals<boost::any, boost::any, boost::any>, Skipper> if_else_stmt;


    InterpreterGrammar() : InterpreterGrammar::base_type(languague) {
        using qi::int_;
        using qi::double_;
        using qi::as_string;
        using qi::lexeme;
        using qi::alnum;
        using qi::alpha;
        using qi::char_;
        using qi::_val;
        using qi::_1;
        using qi::_a;
        using qi::_b;
        using qi::_c;
        using qi::_d;
        using qi::lit;
        using phoenix::construct;
        using phoenix::push_back;
        using phoenix::at;
        using phoenix::at_c;
        using phoenix::ref;
        using phoenix::bind;
        using qi::bool_;
        using qi::eps;
        using qi::skip;
        using qi::_pass;
        using qi::eoi;

        mask %= '[' >> skip(qi::ascii::space)[maskParser] >> ']';

        quoted_string %= as_string[lexeme['"' >> *(char_ - '"') >> '"']];

        identifer = as_string[lexeme[(alpha | char_("_"))
                >> *(alnum | char_("_"))]][_val = construct<Identifer>(_1, _pass)];

        literal %= lexeme[int_ >> !char_('.')] | double_ | bool_ | mask | quoted_string;

        function_call_arguments_expr = eps[_a = construct<std::vector<boost::any>>()]
                >> -(assign_expr[push_back(_a, _1)] >> *(',' >> assign_expr[push_back(_a, _1)]))
                >> eps[_val = _a];

        parenthese_expr = '(' >> assign_expr[_val = _1] >> ')' | literal[_val = _1];

        function_call_expr = parenthese_expr[_val = _1] | identifer[_val = _1]
                >> *('(' >> function_call_arguments_expr[_val = construct<Function>(_val, _1)] >> ')');

        suffix_unary_op_expr = function_call_expr[_val = _1]
                >> -(lit("++")[_val = construct<ArithmeticOperation>(ArithmeticOp::PostIncrement, _val)] |
                     lit("--")[_val = construct<ArithmeticOperation>(ArithmeticOp::PostDecrement, _val)]);

        prefix_unary_op_expr = '!' >> prefix_unary_op_expr[_val = construct<LogicalOperation>(LogicalOp::Not, _1)] |
                               "++" >> prefix_unary_op_expr[_val = construct<ArithmeticOperation>(
                                       ArithmeticOp::PreIncrement, _1)] |
                               "--" >> prefix_unary_op_expr[_val = construct<ArithmeticOperation>(
                                       ArithmeticOp::PreDecrement, _1)] |
                               '+' >> prefix_unary_op_expr[_val = _1] |
                               '-' >> prefix_unary_op_expr[
                                       _val = construct<ArithmeticOperation>(ArithmeticOp::Minus, _1)] |
                               suffix_unary_op_expr[_val = _1];

        arithmetic_high_expr = prefix_unary_op_expr[_val = _1]
                >> *('*' >> prefix_unary_op_expr[
                        _val = construct<ArithmeticOperation>(ArithmeticOp::Multiply, _val, _1)] |
                     '/' >> prefix_unary_op_expr[
                             _val = construct<ArithmeticOperation>(ArithmeticOp::Divide, _val, _1)] |
                     '%' >> prefix_unary_op_expr[_val = construct<ArithmeticOperation>(ArithmeticOp::Mod, _val, _1)]);


        arithmetic_low_expr = arithmetic_high_expr[_val = _1]
                >> *('+' >> arithmetic_high_expr[_val = construct<ArithmeticOperation>(ArithmeticOp::Plus, _val, _1)] |
                     '-' >> arithmetic_high_expr[
                             _val = construct<ArithmeticOperation>(ArithmeticOp::Subtract, _val, _1)]);

        logical_comp_expr = arithmetic_low_expr[_val = _1]
                >> *("<=" >> arithmetic_low_expr[_val = construct<LogicalOperation>(LogicalOp::LessEqual, _val, _1)] |
                     '<' >> arithmetic_low_expr[_val = construct<LogicalOperation>(LogicalOp::Less, _val, _1)] |
                     ">=" >> arithmetic_low_expr[_val = construct<LogicalOperation>(LogicalOp::GreatEqual, _val, _1)] |
                     '>' >> arithmetic_low_expr[_val = construct<LogicalOperation>(LogicalOp::Great, _val, _1)]);

        logical_equal_not_equal_expr = logical_comp_expr[_val = _1]
                >> *("==" >> logical_comp_expr[_val = construct<LogicalOperation>(LogicalOp::Equal, _val, _1)] |
                     "!=" >> logical_comp_expr[_val = construct<LogicalOperation>(LogicalOp::NotEqual, _val, _1)]);

        bitwise_and_expr = logical_equal_not_equal_expr[_val = _1]
                >> *('&' >> logical_equal_not_equal_expr[_val = construct<BitwiseOperation>(BitwiseOp::And, _val, _1)]);

        bitwise_or_expr = bitwise_and_expr[_val = _1]
                >> *('|' >> bitwise_and_expr[_val = construct<BitwiseOperation>(BitwiseOp::Or, _val, _1)]);

        logical_and_expr = bitwise_or_expr[_val = _1]
                >> *("&&" >> bitwise_or_expr[_val = construct<LogicalOperation>(LogicalOp::And, _val, _1)]);

        logical_or_expr = logical_and_expr[_val = _1]
                >> *("||" >> logical_and_expr[_val = construct<LogicalOperation>(LogicalOp::Or, _val, _1)]);

        assign_expr =
                logical_or_expr[_val = _1]
                        >> -('=' >> assign_expr[_val = construct<AssignStmt>(_val, _1)] |
                             "*=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<ArithmeticOperation>(ArithmeticOp::Multiply, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "/=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<ArithmeticOperation>(ArithmeticOp::Divide, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "%=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<ArithmeticOperation>(ArithmeticOp::Mod, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "+=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<ArithmeticOperation>(ArithmeticOp::Plus, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "-=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<ArithmeticOperation>(ArithmeticOp::Subtract, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "&=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<BitwiseOperation>(BitwiseOp::And, _val, _1),
                                     AssignStmt::TYPE::Compound)] |
                             "|=" >> assign_expr[_val = construct<AssignStmt>(
                                     _val, construct<BitwiseOperation>(BitwiseOp::Or, _val, _1),
                                     AssignStmt::TYPE::Compound)]);

        stmt %= for_stmt | if_else_stmt | while_stmt | do_while_until_stmt | assign_expr >> ';';

        stmts = eps[_a = construct<std::vector<boost::any>>()] >> *(stmt[push_back(_a, _1)]) >> eps[_val = _a];

        for_stmt = lit("for") >> "(" >> assign_expr[_a = _1] >> ';'
                              >> assign_expr[_b = _1] >> ';'
                              >> assign_expr[_c = _1] >> ')' >> '{'
                              >> stmts[_d = _1] >> '}'
                              >> eps[_val = construct<ForStmt>(_a, _b, _c, _d)];

        if_else_stmt = lit("if") >> '(' >> assign_expr[_a = _1] >> ')' >> '{' >> stmts[_b = _1] >> '}' >>
                                 -(lit("else") >> '{' >> stmts[_c = _1] >> '}')
                                 >> eps[_val = construct<IfElseStmt>(_a, _b, _c)];

        while_stmt = lit("while") >> '(' >> assign_expr[_a = _1] >> ')' >> '{' >> stmts[_b = _1] >> '}'
                                  >> eps[_val = construct<WhileStmt>(_a, _b)];

        do_while_until_stmt = lit("do") >> '{' >> stmts[_a = _1] >> '}'
                                        >> (lit("until") >> '(' >> assign_expr[_b = _1] >> ')' >> ';'
                                                         >> eps[_val = construct<DoUntilStmt>(_a, _b)] |
                                            lit("while") >> '(' >> assign_expr[_b = _1] >> ')' >> ';'
                                                         >> eps[_val = construct<DoWhileStmt>(_a, _b)]);

        languague %= stmts >> eoi;
    }
};


class InterpreterException : public std::runtime_error {
public:
    explicit InterpreterException(const char *string) : runtime_error(string) {}

    InterpreterException(const InterpreterException &error) : runtime_error(error.what()) {}

    explicit InterpreterException(const std::string &arg) : runtime_error(arg) {}
};

class Interpreter;

class FunctionObject {
    Interpreter *iterpreter{};
public:
    using function_sig =  std::function<
            boost::any(
                    std::vector<
                            std::tuple<
                                    std::string,
                                    std::function<bool(const boost::any)>,
                                    std::function<std::vector<std::string>()>,
                                    boost::any
                            >
                    > &
            )
    >;

    FunctionObject() = default;

    FunctionObject(Interpreter *iterpreter, std::string name,
                   function_sig f)
            : iterpreter(iterpreter), name(std::move(name)), f(std::move(f)) {}

    std::string name;
    std::vector<std::tuple<std::string, std::function<bool(const boost::any)>,
            std::function<std::vector<std::string>()>, boost::any>> argument_mapping;

    function_sig f;

    boost::any
    invoke(std::vector<boost::any> &arguments, std::vector<std::pair<boost::any, boost::any>> &optional_arguments);

    template<typename... Ts>
    FunctionObject &addArgument(const std::string &arg_name, boost::any defaultValue = {}) {
        if (!defaultValue.empty()) {
            throw_assert(TypeIs<Ts...>()(defaultValue),
                         "Arguemnt(" << arg_name << ") default value type error, function (" << name << ")");
        }
        argument_mapping.emplace_back(arg_name, TypeIs<Ts...>(), TypePrettyNames<Ts...>(), defaultValue);
        return *this;
    }
};


class Interpreter {
public:

    friend class FunctionObject;

    boost::any execute(boost::any ast) {
        return evalRightValue(std::move(ast));
    }

    std::unordered_map<std::string, boost::any> &getVariables() {
        return variables;
    }

    FunctionObject &registerFunction(const std::string &name,
                                     const FunctionObject::function_sig &f) {
        functions[name] = FunctionObject(this, name, f);
        return functions[name];
    }

    Interpreter();

protected:

    std::map<std::string, FunctionObject> functions;

    std::unordered_map<std::string, boost::any> variables;

    boost::any evalFunction(const std::string &name, std::vector<boost::any> arguments,
                            std::vector<std::pair<boost::any, boost::any>> optional_arguments);

    boost::any evalRightValue(boost::any ast);

    boost::any evalLeftValue(boost::any v);

    boost::any evalArithmeticOperation(const ArithmeticOperation &op);

    boost::any evalLogicalOperation(const LogicalOperation &op);

    boost::any evalBitwiseOperation(const BitwiseOperation &op);

    boost::any evalArithmeticAdd(boost::any lhs, boost::any rhs);

    boost::any evalArithmeticSubtract(boost::any lhs, boost::any rhs);

    boost::any evalArithmeticMultiply(boost::any lhs, boost::any rhs);

    boost::any evalArithmeticDivide(boost::any lhs, boost::any rhs);

    boost::any evalArithmeticMinus(boost::any lhs);


    bool evalLogicalOperationEqual(boost::any lhs, boost::any rhs);

    bool evalLogicalOperationLess(boost::any lhs, boost::any rhs);

    template<typename T>
    T ArithmeticOpAdd(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) + boost::any_cast<T>(rhs);
    }

    template<typename T>
    T ArithmeticOpSubtract(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) - boost::any_cast<T>(rhs);
    }

    template<typename T>
    T ArithmeticOpMultiply(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) * boost::any_cast<T>(rhs);
    }

    template<typename T>
    T ArithemeticOpDivide(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) / boost::any_cast<T>(rhs);
    }

    template<typename T>
    bool LogicalOpEqaul(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) == boost::any_cast<T>(rhs);
    }

    template<typename T>
    bool LogicalOpLess(boost::any &lhs, boost::any &rhs) {
        return boost::any_cast<T>(lhs) < boost::any_cast<T>(rhs);
    }
};


template<typename Iterator>
struct SkipperGrammar : qi::grammar<Iterator> {

    qi::rule<Iterator> block_comment, single_line_comment, skipper;

    SkipperGrammar() : SkipperGrammar::base_type(skipper) {
        using namespace qi;
        single_line_comment = (lit("//") | "#") >> *(char_ - eol) >> (eol | eoi);
        block_comment = ("/*" >> *(block_comment | char_ - "*/")) > "*/";

        skipper = ascii::space | single_line_comment | block_comment;
    }
};

using SkipperT = SkipperGrammar<std::string::iterator>;
using InterpreterGrammarT = InterpreterGrammar<std::string::iterator, SkipperT>;

#endif //TINKER_INTERPRETER_HPP
