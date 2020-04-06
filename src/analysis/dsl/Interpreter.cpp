//
// Created by xiamr on 7/15/19.
//

#include "dsl/Interpreter.hpp"

#include <boost/spirit/include/karma.hpp>
#include <unordered_set>
#include "dsl/AmberMask.hpp"

#include "utils/ThrowAssert.hpp"

using namespace std;

Interpreter::Interpreter() {
    registerFunction("sqrt", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            return std::sqrt(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(double)) {
            return std::sqrt(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permmited");
    }).addArgument<double, int>("value");

    registerFunction("log", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            return std::log(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(double)) {
            return std::log(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<double, int>("value");

    registerFunction("exp", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            return std::exp(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(double)) {
            return std::exp(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<double, int>("value");

    registerFunction("abs", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            return std::abs(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(double)) {
            return std::abs(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<double, int>("value");

    registerFunction("pow",
                     [this](auto &args) -> boost::any {
                         auto arg1 = evalRightValue(get<3>(args.at(0)));
                         auto arg2 = evalRightValue(get<3>(args.at(1)));
                         if (arg1.type() == typeid(int)) {
                             if (arg2.type() == typeid(int)) {
                                 return std::pow(boost::any_cast<int>(arg1), boost::any_cast<int>(arg2));
                             } else if (arg2.type() == typeid(double)) {
                                 return std::pow(boost::any_cast<int>(arg1), boost::any_cast<double>(arg1));
                             }
                         } else if (arg1.type() == typeid(double)) {
                             if (arg2.type() == typeid(int)) {
                                 return std::pow(boost::any_cast<double>(arg1), boost::any_cast<int>(arg2));
                             } else if (arg2.type() == typeid(double)) {
                                 return std::pow(boost::any_cast<double>(arg1), boost::any_cast<double>(arg1));
                             }
                         }
                         throw InterpreterException("operation not permitted");
                     })
        .addArgument<double, int>("value1")
        .addArgument<double, int>("value2");

    registerFunction("print", [this](auto &args) -> boost::any {
        boost::any arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            std::cout << boost::any_cast<int>(arg) << '\n';
        } else if (arg.type() == typeid(double)) {
            std::cout << boost::any_cast<double>(arg) << '\n';
        } else if (arg.type() == typeid(bool)) {
            std::cout << boost::any_cast<bool>(arg) << '\n';
        } else if (arg.type() == typeid(std::string)) {
            std::cout << boost::any_cast<std::string>(arg) << '\n';
        } else {
            throw InterpreterException("Type cannot print");
        }
        return {};
    }).addArgument<int, double, bool, string>("value");

    registerFunction("int", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(int)) {
            return arg;
        } else if (arg.type() == typeid(double)) {
            return static_cast<int>(boost::any_cast<double>(arg));
        } else if (arg.type() == typeid(bool)) {
            return static_cast<int>(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<int, double, bool>("value");

    registerFunction("double", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(double)) {
            return arg;
        } else if (arg.type() == typeid(int)) {
            return static_cast<double>(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(bool)) {
            return static_cast<double>(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<double, int, bool>("value");

    registerFunction("bool", [this](auto &args) -> boost::any {
        auto arg = evalRightValue(get<3>(args.at(0)));
        if (arg.type() == typeid(bool)) {
            return arg;
        } else if (arg.type() == typeid(int)) {
            return static_cast<bool>(boost::any_cast<int>(arg));
        } else if (arg.type() == typeid(double)) {
            return static_cast<bool>(boost::any_cast<double>(arg));
        }
        throw InterpreterException("operation not permitted");
    }).addArgument<bool, int, double>("value");
}

boost::any Interpreter::evalRightValue(boost::any ast) {
    auto condtion = [](boost::any c) {
        if (!c.empty()) {
            if (c.type() == typeid(int)) {
                return static_cast<bool>(boost::any_cast<int>(c));
            } else if (c.type() == typeid(double)) {
                return static_cast<bool>(boost::any_cast<double>(c));
            } else if (c.type() == typeid(bool)) {
                return boost::any_cast<bool>(c);
            } else {
                return true;
            }
        }
        return false;
    };
    if (ast.empty()) {
        return {};
    } else if (ast.type() == typeid(IfElseStmt)) {
        auto if_else_stmt = boost::any_cast<IfElseStmt>(ast);
        if (condtion(evalRightValue(if_else_stmt.condtion))) {
            evalRightValue(if_else_stmt.if_stmt_block);
        } else {
            evalRightValue(if_else_stmt.else_stmt_block);
        }
        return {};
    } else if (ast.type() == typeid(ForStmt)) {
        auto for_stmt = boost::any_cast<ForStmt>(ast);
        for (evalRightValue(for_stmt.expr1); condtion(evalRightValue(for_stmt.expr2)); evalRightValue(for_stmt.expr3)) {
            try {
                evalRightValue(for_stmt.stmt_block);
            } catch (BreakStmt &) {
                break;
            } catch (ContinueStmt &) {
                continue;
            }
        }
        return {};
    } else if (ast.type() == typeid(WhileStmt)) {
        auto while_stmt = boost::any_cast<WhileStmt>(ast);
        while (condtion(evalRightValue(while_stmt.condition))) {
            try {
                evalRightValue(while_stmt.stmt_block);
            } catch (BreakStmt &) {
                break;
            } catch (ContinueStmt &) {
                continue;
            }
        }
        return {};
    } else if (ast.type() == typeid(DoUntilStmt)) {
        auto do_until_stmt = boost::any_cast<DoUntilStmt>(ast);
        do {
            try {
                evalRightValue(do_until_stmt.stmt_block);
            } catch (BreakStmt &) {
                break;
            } catch (ContinueStmt &) {
                continue;
            }
        } while (!condtion(evalRightValue(do_until_stmt.condition)));
        return {};
    } else if (ast.type() == typeid(DoWhileStmt)) {
        auto do_while_stmt = boost::any_cast<DoWhileStmt>(ast);
        do {
            try {
                evalRightValue(do_while_stmt.stmt_block);
            } catch (BreakStmt &) {
                break;
            } catch (ContinueStmt &) {
                continue;
            }
        } while (condtion(evalRightValue(do_while_stmt.condition)));
        return {};
    } else if (ast.type() == typeid(BreakStmt)) {
        throw BreakStmt{};
    } else if (ast.type() == typeid(ContinueStmt)) {
        throw ContinueStmt{};
    } else if (ast.type() == typeid(Function)) {
        auto f = boost::any_cast<Function>(ast);
        return evalFunction(boost::any_cast<Identifer>(f.name).name, f.arguments, f.optional_arguments);
    } else if (ast.type() == typeid(std::vector<boost::any>)) {
        boost::any ret;
        auto stmts = boost::any_cast<std::vector<boost::any>>(ast);
        for (const auto &stmt : stmts) {
            ret = evalRightValue(stmt);
        }
        return ret;
    } else if (ast.type() == typeid(Identifer)) {
        return variables.at(boost::any_cast<Identifer>(ast).name);
    } else if (ast.type() == typeid(AssignStmt)) {
        auto stmt = boost::any_cast<AssignStmt>(ast);
        auto ret = evalRightValue(evalRightValue(stmt.rhs));
        variables[boost::any_cast<Identifer>(stmt.var).name] = ret;
        return ret;
    } else if (ast.type() == typeid(LogicalOperation)) {
        return evalLogicalOperation(boost::any_cast<LogicalOperation>(ast));
    } else if (ast.type() == typeid(ArithmeticOperation)) {
        return evalRightValue(evalArithmeticOperation(boost::any_cast<ArithmeticOperation>(ast)));
    } else if (ast.type() == typeid(BitwiseOperation)) {
        return evalBitwiseOperation(boost::any_cast<BitwiseOperation>(ast));
    } else {
        return ast;
    }
}

boost::any Interpreter::evalLeftValue(boost::any ast) {
    if (ast.type() == typeid(ArithmeticOperation)) {
        return evalArithmeticOperation(boost::any_cast<ArithmeticOperation>(ast));
    } else if (ast.type() == typeid(Identifer)) {
        return ast;
    }
    return evalRightValue(ast);
}

boost::any Interpreter::evalArithmeticOperation(const ArithmeticOperation &op) {
    auto lhs = evalLeftValue(op.operand1);
    switch (op.op) {
    case ArithmeticOp::PostIncrement:
        if (lhs.type() == typeid(Identifer)) {
            auto value = variables.at(boost::any_cast<Identifer>(lhs).name);
            if (value.type() == typeid(int)) {
                variables[boost::any_cast<Identifer>(lhs).name] = boost::any{boost::any_cast<int>(value) + 1};
            } else if (value.type() == typeid(double)) {
                variables[boost::any_cast<Identifer>(lhs).name] = boost::any{boost::any_cast<double>(value) + 1};
            } else {
                throw InterpreterException("PostIncremnt for this type is not permitted");
            }
            return value;
        } else {
            throw InterpreterException("not Support");
        }
        break;

    case ArithmeticOp::PostDecrement:
        if (lhs.type() == typeid(Identifer)) {
            auto value = variables.at(boost::any_cast<Identifer>(lhs).name);
            if (value.type() == typeid(int)) {
                variables[boost::any_cast<Identifer>(lhs).name] = boost::any{boost::any_cast<int>(value) - 1};
            } else if (value.type() == typeid(double)) {
                variables[boost::any_cast<Identifer>(lhs).name] = boost::any{boost::any_cast<double>(value) - 1};
            } else {
                throw InterpreterException("PostDecremnt for this type is not permitted");
            }
            return value;
        } else {
            throw InterpreterException("not Support");
        }
        break;

    case ArithmeticOp::PreIncrement:
        if (lhs.type() == typeid(Identifer)) {
            auto value = variables.at(boost::any_cast<Identifer>(lhs).name);
            boost::any ret;
            if (value.type() == typeid(int)) {
                ret = boost::any_cast<int>(value) + 1;
            } else if (value.type() == typeid(double)) {
                ret = boost::any_cast<double>(value) + 1;
            } else {
                throw InterpreterException("PostIncremnt for this type is not permitted");
            }
            variables[boost::any_cast<Identifer>(lhs).name] = ret;
            return lhs;
        } else {
            throw InterpreterException("not Support");
        }
        break;

    case ArithmeticOp::PreDecrement:
        if (lhs.type() == typeid(Identifer)) {
            auto value = variables.at(boost::any_cast<Identifer>(lhs).name);
            boost::any ret;
            if (value.type() == typeid(int)) {
                ret = boost::any_cast<int>(value) - 1;
            } else if (value.type() == typeid(double)) {
                ret = boost::any_cast<double>(value) - 1;
            } else {
                throw InterpreterException("PostIncremnt for this type is not permitted");
            }
            variables[boost::any_cast<Identifer>(lhs).name] = ret;
            return lhs;
        } else {
            throw InterpreterException("not Support ");
        }
        break;
    default:
        break;
    }

    lhs = evalRightValue(lhs);
    auto rhs = evalRightValue(op.operand2);
    switch (op.op) {
    case ArithmeticOp::Plus:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else if (lhs.type() == typeid(std::string)) {
                if (rhs.type() == typeid(int)) {
                    rhs = std::to_string(boost::any_cast<int>(rhs));
                };
            } else if (rhs.type() == typeid(std::string)) {
                if (lhs.type() == typeid(int)) {
                    lhs = std::to_string(boost::any_cast<int>(lhs));
                };
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalArithmeticAdd(lhs, rhs);
        break;
    case ArithmeticOp::Minus:
        return evalArithmeticMinus(lhs);
        break;
    case ArithmeticOp::Multiply:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalArithmeticMultiply(lhs, rhs);
        break;
    case ArithmeticOp::Divide:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalArithmeticDivide(lhs, rhs);
        break;
    case ArithmeticOp::Subtract:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalArithmeticSubtract(lhs, rhs);
        break;

    case ArithmeticOp::Mod:
        if (lhs.type() != rhs.type()) {
            throw InterpreterException("Type not same");
        }
        if (lhs.type() != typeid(int)) {
            throw InterpreterException("Only integer can do %");
        }
        return boost::any_cast<int>(lhs) % boost::any_cast<int>(rhs);
        break;
    default:
        throw InterpreterException("Operation not valid !!");
    }
}

boost::any Interpreter::evalBitwiseOperation(const BitwiseOperation &op) {
    auto lhs = evalRightValue(op.operand1);
    auto rhs = evalRightValue(op.operand2);
    if (lhs.type() == typeid(AmberMask) && rhs.type() == typeid(AmberMask)) {
        switch (op.op) {
        case BitwiseOp::And:
            return AmberMask(std::make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::AND, boost::any_cast<AmberMask>(lhs),
                                                               boost::any_cast<AmberMask>(rhs)));
            break;
        case BitwiseOp::Or:
            return AmberMask(std::make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::OR, boost::any_cast<AmberMask>(lhs),
                                                               boost::any_cast<AmberMask>(rhs)));
            break;
        }
    }
    throw InterpreterException("Bitwise operation only for AmbberMask");
}
namespace {
inline ostream &operator<<(ostream &os, const vector<string> &vect) {
    boost::spirit::karma::generate(boost::spirit::karma::ostream_iterator<char>(os),
                                   "vector(" << boost::spirit::karma::string % ',' << ')', vect);
    return os;
}
} // namespace

boost::any FunctionObject::invoke(vector<boost::any> &arguments,
                                  vector<pair<boost::any, boost::any>> &optional_arguments) {
    vector<tuple<string, function<bool(const boost::any)>, function<vector<string>()>, boost::any>> args(
        argument_mapping);

    if (arguments.size() > args.size()) {
        std::cerr << "function(" << name << ") too may arguments, except " << args.size() << " , actual "
                  << arguments.size() << '\n';
        exit(EXIT_FAILURE);
    }

    for (std::size_t i = 0; i < arguments.size(); i++) {
        std::get<3>(args[i]) = arguments[i];
    }

    std::unordered_set<std::string> s;
    for (const auto &p : optional_arguments) {
        auto arg_name = boost::any_cast<Identifer>(p.first).name;
        bool bFound = false;
        for (auto &arg : args) {
            if (arg_name == std::get<0>(arg)) {
                if (s.count(arg_name)) {
                    std::cerr << "named argument " << arg_name << " occur more than once\n";
                    exit(EXIT_FAILURE);
                } else {
                    s.insert(arg_name);
                }
                std::get<3>(arg) = p.second;
                bFound = true;
                break;
            }
        }
        if (!bFound) {
            std::cerr << boost::format("argment name(%s) not found in function(%s)\n") % arg_name % name;
            exit(EXIT_FAILURE);
        }
    }

    for (auto &arg : args) {
        std::get<3>(arg) = iterpreter->evalRightValue(std::get<3>(arg));
        throw_assert(!std::get<3>(arg).empty(), "argment name (" + std::get<0>(arg) + "do not have value");
        throw_assert(std::get<1>(arg)(std::get<3>(arg)),
                     "argument (" << std::get<0>(arg) << ") type not match , except " << std::get<2>(arg)() << " Actul"
                                  << getPrettyName(std::get<3>(arg)));
    }
    std::cout << std::vector<std::string>();
    return f(args);
}

boost::any Interpreter::evalLogicalOperation(const LogicalOperation &op) {
    auto lhs = evalRightValue(op.operand1);
    auto rhs = evalRightValue(op.operand2);
    switch (op.op) {
    case LogicalOp::Equal:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalLogicalOperationEqual(lhs, rhs);
        break;
    case LogicalOp::NotEqual:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return !evalLogicalOperationEqual(lhs, rhs);
        break;
    case LogicalOp::And:
        if (lhs.type() == typeid(int)) {
            lhs = static_cast<bool>(boost::any_cast<int>(lhs));
        } else if (lhs.type() == typeid(double)) {
            lhs = static_cast<bool>(boost::any_cast<double>(lhs));
        }

        if (rhs.type() == typeid(int)) {
            rhs = static_cast<bool>(boost::any_cast<int>(rhs));
        } else if (lhs.type() == typeid(double)) {
            rhs = static_cast<bool>(boost::any_cast<double>(rhs));
        }

        if (lhs.type() != typeid(bool) || rhs.type() != typeid(bool)) {
            throw InterpreterException("only bool can do &&");
        }
        return boost::any_cast<bool>(lhs) && boost::any_cast<bool>(rhs);
        break;
    case LogicalOp::Or:
        if (lhs.type() == typeid(int)) {
            lhs = static_cast<bool>(boost::any_cast<int>(lhs));
        } else if (lhs.type() == typeid(double)) {
            lhs = static_cast<bool>(boost::any_cast<double>(lhs));
        }

        if (rhs.type() == typeid(int)) {
            rhs = static_cast<bool>(boost::any_cast<int>(rhs));
        } else if (lhs.type() == typeid(double)) {
            rhs = static_cast<bool>(boost::any_cast<double>(rhs));
        }

        if (lhs.type() != typeid(bool) || rhs.type() != typeid(bool)) {
            throw InterpreterException("only bool can do ||");
        }
        return boost::any_cast<bool>(lhs) || boost::any_cast<bool>(rhs);
        break;
    case LogicalOp::LessEqual:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalLogicalOperationLess(lhs, rhs) || evalLogicalOperationEqual(lhs, rhs);
        break;
    case LogicalOp::Less:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return evalLogicalOperationLess(lhs, rhs);
        break;
    case LogicalOp::GreatEqual:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return !evalLogicalOperationLess(lhs, rhs);
        break;
    case LogicalOp::Great:
        if (lhs.type() != rhs.type()) {
            if (lhs.type() == typeid(int) && rhs.type() == typeid(double)) {
                lhs = static_cast<double>(boost::any_cast<int>(lhs));
            } else if (lhs.type() == typeid(double) && rhs.type() == typeid(int)) {
                rhs = static_cast<double>(boost::any_cast<int>(rhs));
            } else {
                throw InterpreterException("Type not same");
            }
        }
        return !(evalLogicalOperationLess(lhs, rhs) || evalLogicalOperationEqual(lhs, rhs));
        break;
    case LogicalOp::Not:
        if (lhs.type() == typeid(int)) {
            lhs = static_cast<bool>(boost::any_cast<int>(lhs));
        } else if (lhs.type() == typeid(double)) {
            lhs = static_cast<bool>(boost::any_cast<double>(lhs));
        }

        if (lhs.type() == typeid(bool)) {
            return !boost::any_cast<bool>(lhs);
        } else if (lhs.type() == typeid(AmberMask)) {
            return AmberMask(std::make_shared<AmberMaskAST::Operator>(AmberMaskAST::Op::NOT, boost::any_cast<AmberMask>(lhs)));
        } else {
            throw InterpreterException(std::string("Logical Not operation not permitted for type : ") +
                                       lhs.type().name());
        }
        break;
    }
    return {};
}

boost::any Interpreter::evalArithmeticAdd(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(bool)) {
        return ArithmeticOpAdd<bool>(lhs, rhs);
    } else if (lhs.type() == typeid(int)) {
        return ArithmeticOpAdd<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return ArithmeticOpAdd<double>(lhs, rhs);
    } else if (lhs.type() == typeid(std::string)) {
        return ArithmeticOpAdd<std::string>(lhs, rhs);
    } else {
        throw InterpreterException("type cannot do add");
    }
}

boost::any Interpreter::evalArithmeticSubtract(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(bool)) {
        return ArithmeticOpSubtract<bool>(lhs, rhs);
    } else if (lhs.type() == typeid(int)) {
        return ArithmeticOpSubtract<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return ArithmeticOpSubtract<double>(lhs, rhs);
    } else {
        throw InterpreterException("type cannot do subtract");
    }
}

boost::any Interpreter::evalArithmeticMultiply(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(bool)) {
        return ArithmeticOpMultiply_bool(lhs, rhs);
    } else if (lhs.type() == typeid(int)) {
        return ArithmeticOpMultiply<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return ArithmeticOpMultiply<double>(lhs, rhs);
    } else {
        throw InterpreterException("type cannot do multiply");
    }
}

boost::any Interpreter::evalArithmeticDivide(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(bool)) {
        return ArithemeticOpDivide<bool>(lhs, rhs);
    } else if (lhs.type() == typeid(int)) {
        return ArithemeticOpDivide<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return ArithemeticOpDivide<double>(lhs, rhs);
    } else {
        throw InterpreterException("type cannot do divide");
    }
}

boost::any Interpreter::evalArithmeticMinus(boost::any lhs) {
    if (lhs.type() == typeid(bool)) {
        return -boost::any_cast<bool>(lhs);
    } else if (lhs.type() == typeid(int)) {
        return -boost::any_cast<int>(lhs);
    } else if (lhs.type() == typeid(double)) {
        return -boost::any_cast<double>(lhs);
    } else {
        throw InterpreterException("type cannot do minus");
    }
}

bool Interpreter::evalLogicalOperationEqual(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(bool)) {
        return LogicalOpEqaul<bool>(lhs, rhs);
    } else if (lhs.type() == typeid(int)) {
        return LogicalOpEqaul<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return LogicalOpEqaul<double>(lhs, rhs);
    } else {
        throw InterpreterException("type not comparable");
    }
}

bool Interpreter::evalLogicalOperationLess(boost::any lhs, boost::any rhs) {
    if (lhs.type() == typeid(int)) {
        return LogicalOpLess<int>(lhs, rhs);
    } else if (lhs.type() == typeid(double)) {
        return LogicalOpLess<double>(lhs, rhs);
    } else {
        throw InterpreterException("type not comparable");
    }
}

boost::any Interpreter::evalFunction(const std::string &name, std::vector<boost::any> arguments,
                                     std::vector<std::pair<boost::any, boost::any>> optional_arguments) {
    try {
        return functions.at(name).invoke(arguments, optional_arguments);
    } catch (std::out_of_range &e) {
        std::cerr << "function " << name << " not found\n";
        exit(EXIT_FAILURE);
    }
}

Function::Function(boost::any name, boost::any arguments) : name(std::move(name)) {
    auto args = boost::any_cast<std::vector<boost::any>>(arguments);

    bool postional = true;

    for (auto &arg : args) {
        if (postional) {
            if (arg.type() == typeid(AssignStmt)) {
                postional = false;
                auto a = boost::any_cast<AssignStmt>(arg);
                if (a.t == AssignStmt::TYPE::Compound) {
                    std::cerr << "Compound assignment expression not allowed in function argment\n";
                    exit(EXIT_FAILURE);
                }
                optional_arguments.emplace_back(a.var, a.rhs);
            } else {
                this->arguments.emplace_back(arg);
            }
        } else {
            if (arg.type() == typeid(AssignStmt)) {
                auto a = boost::any_cast<AssignStmt>(arg);
                if (a.t == AssignStmt::TYPE::Compound) {
                    std::cerr << "Compound assignment expression not allowed in function argment\n";
                    exit(EXIT_FAILURE);
                }
                optional_arguments.emplace_back(a.var, a.rhs);
            } else {
                std::cerr << "Function Argurment order error, named arguments must after position arguments"
                             "function name :"
                          << boost::any_cast<Identifer>(name).name << '\n';
                exit(EXIT_FAILURE);
            }
        }
    }
}
