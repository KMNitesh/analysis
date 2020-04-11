

#include "PrintTopolgy.hpp"
#include "data_structure/atom.hpp"
#include "data_structure/forcefield.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"
#include "dsl/selection_grammar.hpp"
#include "trajectory_reader/trajectoryreader.hpp"
#include "utils/BondEnergyCalculator.hpp"
#include "utils/EnergyCalculator.hpp"
#include "utils/PBCUtils.hpp"
#include <boost/algorithm/cxx11/one_of.hpp>
#include <boost/range/algorithm.hpp>
#include <fnmatch.h>
#include <regex>

namespace {

std::set<std::string> allowed_field_names{"atom_name", "atom_no", "atom_type", "res_name",
                                          "res_no",    "mol_no",  "charge",    "mass"};

class WhereInterpreter {
public:
    virtual bool eval([[maybe_unused]] const std::shared_ptr<Atom> &atom) { return true; };
};

template <typename Func> class NameEqual : public WhereInterpreter {
public:
    NameEqual(std::string name) : name(std::move(name)) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) { return func(atom) == name; }

private:
    std::string name;
    Func func;
};

template <typename Func> class NameLike : public WhereInterpreter {
public:
    NameLike(std::string pattern) : pattern(std::move(pattern)) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) {
        return fnmatch(pattern.c_str(), func(atom).c_str(), 0) == 0;
    }

private:
    std::string pattern;
    Func func;
};

template <typename Func> class NameRegex : public WhereInterpreter {
public:
    NameRegex(const std::string &regex) : regex(regex) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) { return std::regex_match(func(atom), regex); }

private:
    std::regex regex;
    Func func;
};

template <typename Func> class BetweenOp : public WhereInterpreter {
public:
    BetweenOp(double low_bound, double up_bound) : low_bound(low_bound), up_bound(up_bound) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) {
        auto val = func(atom);
        return val >= low_bound and val <= up_bound;
    }

private:
    double low_bound, up_bound;
    Func func;
};

template <typename Func> class InOp : public WhereInterpreter {
public:
    InOp(std::vector<std::string> fs) : fs(std::move(fs)) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) { return boost::algorithm::one_of_equal(fs, func(atom)); }

private:
    std::vector<std::string> fs;
    Func func;
};

template <typename Func, typename T = std::result_of_t<Func(const std::shared_ptr<Atom> &)>>
class NoCompare : public WhereInterpreter {
public:
    NoCompare(LogicalOperator::Op op, T no) : op(op), no(no) {}
    virtual bool eval(const std::shared_ptr<Atom> &atom) {
        switch (op) {
        case LogicalOperator::Op::Less:
            return func(atom) < no;
        case LogicalOperator::Op::LessEqual:
            return func(atom) <= no;
        case LogicalOperator::Op::Equal:
            return func(atom) == no;
        case LogicalOperator::Op::Great:
            return func(atom) > no;
        case LogicalOperator::Op::GreatEqual:
            return func(atom) >= no;
        default:
            throw std::runtime_error("Unknown LogicalOperator::Op");
        }
    }

private:
    LogicalOperator::Op op;
    T no;
    Func func;
};

auto atom_name_getter = [](auto &atom) { return atom->atom_name; };
auto atom_no_getter = [](auto &atom) { return atom->seq; };
auto atom_type_getter = [](auto &atom) { return atom->type_name; };
auto res_name_getter = [](auto &atom) { return atom->residue_name.value(); };
auto res_no_getter = [](auto &atom) { return atom->residue_num.value(); };
auto mol_no_getter = [](auto &atom) { return atom->molecule.lock()->sequence; };
auto charge_getter = [](auto &atom) { return atom->charge.value(); };
auto mass_getter = [](auto &atom) { return atom->mass.value(); };

using AtomNameEqual = NameEqual<decltype(atom_name_getter)>;
using AtomNameLike = NameLike<decltype(atom_name_getter)>;
using AtomNameRegex = NameRegex<decltype(atom_name_getter)>;
using AtomNameInOp = InOp<decltype(atom_name_getter)>;

using AtomNoCompare = NoCompare<decltype(atom_no_getter)>;
using AtomNoBetweenOp = BetweenOp<decltype(atom_no_getter)>;

using AtomTypeEqual = NameEqual<decltype(atom_type_getter)>;
using AtomTypeRegex = NameRegex<decltype(atom_type_getter)>;
using AtomTypeLike = NameLike<decltype(atom_type_getter)>;
using AtomTypeInOp = InOp<decltype(atom_type_getter)>;

using ResidueNameEqual = NameEqual<decltype(res_name_getter)>;
using ResidueNameLike = NameLike<decltype(res_name_getter)>;
using ResidueNameRegex = NameRegex<decltype(res_name_getter)>;
using ResidueNameInOp = InOp<decltype(res_name_getter)>;

using ResidueNoCompare = NoCompare<decltype(res_no_getter)>;
using ResidueNoBetweenOp = BetweenOp<decltype(res_no_getter)>;

using MoleculeNoCompare = NoCompare<decltype(mol_no_getter)>;
using MoleculeNoBetweenOp = BetweenOp<decltype(mol_no_getter)>;

using ChargeCompare = NoCompare<decltype(charge_getter)>;
using ChargeBetweenOp = BetweenOp<decltype(charge_getter)>;

using MassCompare = NoCompare<decltype(mass_getter)>;
using MassBetweenOp = BetweenOp<decltype(mass_getter)>;

class CombineCompare : public WhereInterpreter {
public:
    CombineCompare(CombineOperator::Op op, std::unique_ptr<WhereInterpreter> cond1,
                   std::unique_ptr<WhereInterpreter> cond2 = {})
        : op(op), cond1(std::move(cond1)), cond2(std::move(cond2)) {}

    virtual bool eval(const std::shared_ptr<Atom> &atom) {
        switch (op) {
        case CombineOperator::Op::AND:
            return cond1->eval(atom) and cond2->eval(atom);
        case CombineOperator::Op::OR:
            return cond1->eval(atom) or cond2->eval(atom);
        case CombineOperator::Op::NOT:
            return not cond1->eval(atom);
        default:
            throw std::runtime_error("Unknown CombineOperator::Op");
        }
    }

private:
    CombineOperator::Op op;
    std::unique_ptr<WhereInterpreter> cond1;
    std::unique_ptr<WhereInterpreter> cond2;
};

class WhereClauseChecker : public boost::static_visitor<std::unique_ptr<WhereInterpreter>> {
public:
    std::unique_ptr<WhereInterpreter> operator()(const boost::blank &) const { return {}; }
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<CombineOperator> &condition) const;
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<LogicalOperator> &condition) const;
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<LikeOperator> &condition) const;
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<RegexOperator> &condition) const;
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<BetweenOperator> &condition) const;
    std::unique_ptr<WhereInterpreter> operator()(const std::shared_ptr<InOperator> &condition) const;
};

std::unique_ptr<WhereInterpreter>
WhereClauseChecker::operator()(const std::shared_ptr<CombineOperator> &condition) const {
    return std::make_unique<CombineCompare>(condition->op, boost::apply_visitor(*this, condition->cond1),
                                            boost::apply_visitor(*this, condition->cond2));
}
std::unique_ptr<WhereInterpreter>
WhereClauseChecker::operator()(const std::shared_ptr<LogicalOperator> &condition) const {
    const auto &id = condition->id;
    const auto op = condition->op;
    const auto &val = condition->value;
    if (id == "atom_name" and op == LogicalOperator::Op::Equal and val.type() == typeid(std::string)) {
        return std::make_unique<AtomNameEqual>(boost::any_cast<std::string>(val));
    } else if (id == "atom_no" and val.type() == typeid(int)) {
        return std::make_unique<AtomNoCompare>(op, boost::any_cast<int>(val));
    } else if (id == "atom_type" and op == LogicalOperator::Op::Equal and val.type() == typeid(std::string)) {
        return std::make_unique<AtomTypeEqual>(boost::any_cast<std::string>(val));
    } else if (id == "res_name" and op == LogicalOperator::Op::Equal and val.type() == typeid(std::string)) {
        return std::make_unique<ResidueNameEqual>(boost::any_cast<std::string>(val));
    } else if (id == "res_no" and val.type() == typeid(int)) {
        return std::make_unique<ResidueNoCompare>(op, boost::any_cast<int>(val));
    } else if (id == "mol_no" and val.type() == typeid(int)) {
        return std::make_unique<MoleculeNoCompare>(op, boost::any_cast<int>(val));
    } else if (id == "charge") {
        if (val.type() == typeid(int)) {
            return std::make_unique<ChargeCompare>(op, boost::any_cast<int>(val));
        } else if (val.type() == typeid(double)) {
            return std::make_unique<ChargeCompare>(op, boost::any_cast<double>(val));
        }
    } else if (id == "mass") {
        if (val.type() == typeid(int)) {
            return std::make_unique<MassCompare>(op, boost::any_cast<int>(val));
        } else if (val.type() == typeid(double)) {
            return std::make_unique<MassCompare>(op, boost::any_cast<double>(val));
        }
    }
    throw std::runtime_error("[where clause grammar error]");
}

std::unique_ptr<WhereInterpreter> WhereClauseChecker::operator()(const std::shared_ptr<LikeOperator> &condition) const {
    if (condition->id == "atom_name") {
        return std::make_unique<AtomNameLike>(condition->pattern);
    } else if (condition->id == "atom_type") {
        return std::make_unique<AtomTypeLike>(condition->pattern);
    } else if (condition->id == "res_name") {
        return std::make_unique<ResidueNameLike>(condition->pattern);
    }
    throw std::runtime_error("[where clause grammar error]");
}

std::unique_ptr<WhereInterpreter>
WhereClauseChecker::operator()(const std::shared_ptr<RegexOperator> &condition) const {
    if (condition->id == "atom_name") {
        return std::make_unique<AtomNameRegex>(condition->regex);
    } else if (condition->id == "atom_type") {
        return std::make_unique<AtomTypeRegex>(condition->regex);
    } else if (condition->id == "res_name") {
        return std::make_unique<ResidueNameRegex>(condition->regex);
    }
    throw std::runtime_error("[where clause grammar error]");
}

std::unique_ptr<WhereInterpreter>
WhereClauseChecker::operator()(const std::shared_ptr<BetweenOperator> &condition) const {
    if (condition->id == "atom_no") {
        return std::make_unique<AtomNoBetweenOp>(condition->low_bound, condition->up_bound);
    } else if (condition->id == "res_no") {
        return std::make_unique<ResidueNoBetweenOp>(condition->low_bound, condition->up_bound);
    } else if (condition->id == "mol_no") {
        return std::make_unique<MoleculeNoBetweenOp>(condition->low_bound, condition->up_bound);
    } else if (condition->id == "charge") {
        return std::make_unique<ChargeBetweenOp>(condition->low_bound, condition->up_bound);
    } else if (condition->id == "mass") {
        return std::make_unique<MassBetweenOp>(condition->low_bound, condition->up_bound);
    }
    throw std::runtime_error("[where clause grammar error]");
}

std::unique_ptr<WhereInterpreter> WhereClauseChecker::operator()(const std::shared_ptr<InOperator> &condition) const {
    if (condition->id == "atom_name") {
        return std::make_unique<AtomNameInOp>(condition->fs);
    } else if (condition->id == "atom_type") {
        return std::make_unique<AtomTypeInOp>(condition->fs);
    } else if (condition->id == "res_name") {
        return std::make_unique<ResidueNameInOp>(condition->fs);
    }
    throw std::runtime_error("[where clause grammar error]");
}

template <typename Func> class OrderBySort {
public:
    OrderBySort(Func func, OrderBy::Order order) : func(std::move(func)), order(order) {}
    template <typename T> bool operator()(const T &lhs, const T &rhs) {
        return order == OrderBy::Order::ASC ? std::less<>()(func(lhs), func(rhs)) : std::less<>()(func(rhs), func(lhs));
    }

private:
    Func func;
    OrderBy::Order order;
};

} // namespace

void PrintTopolgy::action(const std::string &topology_filename) {
    TrajectoryReader reader;

    auto t1 = std::chrono::high_resolution_clock::now();

    reader.set_topology(topology_filename);
    auto frame = reader.readTopology();
    if (forcefield.isValid()) {
        forcefield.assign_forcefield(frame);
    }

    enum class Mode { Noop, Mass, Geom } mode = Mode::Noop;

    std::cout
        << std::left << "Total atom numbers : " << std::setw(10) << frame->atom_list.size()
        << "Total molecule numbers : " << std::setw(10) << frame->molecule_list.size() << "Load time : "
        << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count()
        << " ms\n";

    for (;;) {
    beg:
        SelectionAST r;
        selectCentergroup(r, "> ");
        AmberMask ast;

        std::vector<std::shared_ptr<Atom>> selected_atoms;

        if (auto ret = boost::get<MassCenterRuleNode>(&r)) {
            ast = ret->SelectionMask;
            mode = Mode::Mass;
        } else if (auto ret = boost::get<GeomCenterRuleNode>(&r)) {
            ast = ret->SelectionMask;
            mode = Mode::Geom;
        } else if (auto ret = boost::get<NoopRuleNode>(&r)) {
            ast = ret->SelectionMask;
            mode = Mode::Noop;
        } else if (auto ret = boost::get<EDARuleNode>(&r)) {
            EnergyCalculator calculator;
            calculator.setMask(ret->mask1, ret->mask2, frame);
            auto [ele, lj] = calculator.calculate_energy(frame);
            std::cout << "E(ele) = " << ele << " (kcal/mole)\n"
                      << "E(vdW) = " << lj << " (kcal/mole)\n";
            continue;
        } else if (auto ret = boost::get<BondRuleNode>(&r)) {
            std::array atoms{PBCUtils::find_atom(ret->mask1, frame), PBCUtils::find_atom(ret->mask2, frame)};
            boost::sort(atoms);
            if (auto it = frame->f_bond_params.find(atoms); it != std::end(frame->f_bond_params)) {
                auto distance = atom_distance(atoms, frame);
                auto r = distance - it->second.rA;
                auto energy = 0.5 * it->second.krA * r * r;
                std::cout
                    << boost::format(
                           " K = %10g kcal/mol/Ang2    b0 = %10g Ang    b = %10g Ang   E(bond) = %10g kcal/mol\n") %
                           it->second.krA % it->second.rA % distance % energy;
            } else {
                std::cout << "No bond parameter found !!!\n";
            }
            continue;
        } else if (auto ret = boost::get<AngleRuleNode>(&r)) {
            std::array atoms{PBCUtils::find_atom(ret->mask1, frame), PBCUtils::find_atom(ret->mask2, frame),
                             PBCUtils::find_atom(ret->mask3, frame)};

            boost::sort(atoms, [](const auto &lhs, const auto &rhs) { return lhs->seq < rhs->seq; });

            if (auto it = frame->f_angle_params.find(atoms); it != std::end(frame->f_angle_params)) {
                auto angle = atom_angle(atoms, frame);
                auto theta = (angle - it->second.rA) * degree;
                auto energy = 0.5 * it->second.krA * theta * theta;
                std::cout << boost::format(" K = %10g kcal/mol/rad2    theta0 = %10g degree    theta = %10g degree  "
                                           "E(angle) = %10g kcal/mol\n") %
                                 it->second.krA % it->second.rA % angle % energy;
            } else {
                std::cout << "No angle parameter found !!!\n";
            }
            continue;

        } else if (auto ret = boost::get<DihedralRuleNode>(&r)) {
            std::array atoms{PBCUtils::find_atom(ret->mask1, frame), PBCUtils::find_atom(ret->mask2, frame),
                             PBCUtils::find_atom(ret->mask3, frame), PBCUtils::find_atom(ret->mask4, frame)};

            boost::sort(atoms, [](const auto &lhs, const auto &rhs) { return lhs->seq < rhs->seq; });

            if (auto [lower_bound, upper_bound] = frame->f_dihedral_params.equal_range(atoms);
                lower_bound != upper_bound) {
                auto angle = atom_dihedral(atoms, frame);
                double energy{};
                for (auto it = lower_bound; it != upper_bound; ++it) {
                    std::cout << boost::format(" K = %10g kcal/mol    phi0 = %10g degree    multi = %2d\n") %
                                     it->second.cpA % it->second.phiA % it->second.mult;
                    auto cos = std::cos((angle * it->second.mult - it->second.phiA) * degree);
                    energy += it->second.cpA * (1 + cos);
                }
                std::cout << boost::format(" phi = %10g degree  E(dihedral) = %10g kcal/mol\n") % angle % energy;
            } else {
                std::cout << "No dihedral parameter found !!!\n";
            }
            continue;
        } else if (auto ret = boost::get<BondedEnergyRuleNode>(&r)) {
            auto [bond, angle, dihedral, improper] = BondEnergyCalculator::energy(ret->mask, frame);
            std::cout << boost::format("E(bond) = %g Kcal/mol   E(angle) = %g kcal/mol   E(dihedral) = %g kcal/mol "
                                       "E(improper dihedral) = %g kcal/mol   E_total = %g kcal/mol\n") %
                             bond % angle % dihedral % improper % (bond + angle + dihedral + improper);
            continue;
        } else if (auto ret = boost::get<SelectStmt>(&r)) {
            const auto &stmt = *ret;
            if (stmt.fields.empty()) {
                std::cerr << "select clause with empty field names\n";
                goto beg;
            }
            if (auto it = boost::find(stmt.fields, "*"); it != end(stmt.fields) and stmt.fields.size() != 1) {
                std::cerr << "cannot combine * with other field names\n";
                goto beg;
            }

            for (const auto &name : stmt.fields) {
                if (!allowed_field_names.contains(name) and name != "*") {
                    std::cerr << "Unknown field name :" << name << '\n'
                              << "all allowed names : atom_name, atom_no, atom_type, res_name, res_no, mol_no, charge, "
                                 "mass\n";
                    goto beg;
                }
            }
            std::unique_ptr<WhereInterpreter> where_clause;
            if (stmt.where.has_value()) {
                try {
                    where_clause = boost::apply_visitor(WhereClauseChecker(), stmt.where.get());
                } catch (std::runtime_error &e) {
                    std::cerr << e.what();
                    goto beg;
                }
            } else {
                where_clause = std::make_unique<WhereInterpreter>();
            }

            for (auto &atom : frame->atom_list) {
                if (is_match(atom, stmt.mask) and where_clause->eval(atom)) {
                    selected_atoms.push_back(atom);
                }
            }
            if (stmt.orderby.has_value()) {
                const auto &field = stmt.orderby.get().field;
                const auto order = stmt.orderby.get().order;
                if (field == "atom_name") {
                    boost::sort(selected_atoms, OrderBySort(atom_name_getter, order));
                } else if (field == "atom_no") {
                    boost::sort(selected_atoms, OrderBySort(atom_no_getter, order));
                } else if (field == "atom_type") {
                    boost::sort(selected_atoms, OrderBySort(atom_type_getter, order));
                } else if (field == "res_name") {
                    boost::sort(selected_atoms, OrderBySort(res_name_getter, order));
                } else if (field == "res_no") {
                    boost::sort(selected_atoms, OrderBySort(res_no_getter, order));
                } else if (field == "mol_no") {
                    boost::sort(selected_atoms, OrderBySort(mol_no_getter, order));
                } else if (field == "charge") {
                    boost::sort(selected_atoms, OrderBySort(charge_getter, order));
                } else if (field == "mass") {
                    boost::sort(selected_atoms, OrderBySort(mass_getter, order));
                } else {
                    std::cerr << "unknown order by field\n";
                    goto beg;
                }
            }
            if (stmt.limit.has_value()) {
                selected_atoms.resize(stmt.limit.get());
            }
            if (!(stmt.fields.size() == 1 and stmt.fields.front() == "*")) {
                // atom_name, atom_no, atom_type, res_name, res_no, mol_no, charge, mass
                for (auto &f : stmt.fields) {
                    std::cout << boost::format("%15s") % f;
                }
                std::cout << '\n';

                for (auto &atom : selected_atoms) {
                    for (auto &f : stmt.fields) {
                        if (f == "atom_name") {
                            std::cout << boost::format("%15s") % atom_name_getter(atom);
                        } else if (f == "atom_no") {
                            std::cout << boost::format("%15d") % atom_no_getter(atom);
                        } else if (f == "atom_type") {
                            std::cout << boost::format("%15s") % atom_type_getter(atom);
                        } else if (f == "res_name") {
                            std::cout << boost::format("%15s") % res_name_getter(atom);
                        } else if (f == "res_no") {
                            std::cout << boost::format("%15d") % res_no_getter(atom);
                        } else if (f == "mol_no") {
                            std::cout << boost::format("%15d") % mol_no_getter(atom);
                        } else if (f == "charge") {
                            std::cout << boost::format("%15.3f") % charge_getter(atom);
                        } else if (f == "mass") {
                            std::cout << boost::format("%15.3f") % mass_getter(atom);
                        }
                    }
                    std::cout << '\n';
                }
                continue;
            }
        } else if (boost::get<QuitRuleNode>(&r)) {
            break;
        } else if (auto ret = boost::get<HelpRuleNode>(&r)) {
            if (ret->keyword.has_value()) {
                const auto &keyword = ret->keyword.get();
                if (keyword == "com") {
                    std::cout << "syntax : com mask\n"
                              << "list selected atoms' information and center of mass\n";
                } else if (keyword == "geom") {
                    std::cout << "syntax : geom mask\n"
                              << "list selected atoms' information and center of geometry\n";
                } else if (keyword == "eda") {
                    std::cout << "syntax : eda mask1 mask2\n"
                              << "calculate electrostatic and vdW energies between two groups\n";
                } else if (keyword == "bond") {
                    std::cout << "syntax : bond mask1 mask2\n"
                              << " list bond parameters, bond length and energy\n";
                } else if (keyword == "angle") {
                    std::cout << "syntax : angle mask1 mask2\n"
                              << " list angle parameters, angle degree and energy\n";
                } else if (keyword == "dihedral") {
                    std::cout << "syntax : dihedral mask1 mask2 mask3 mask4\n"
                              << " list dihedral angle parameters, dihedral angle degree and energy\n";
                } else if (keyword == "energy") {
                    std::cout << "syntax : energy mask\n"
                              << " calculate bond, angle, dihedral and improper energy\n";
                } else if (keyword == "help") {
                    std::cout << "syntax : help [keyword]\n"
                              << "print help information with topic of sepcific keyword or help menu\n";
                } else if (keyword == "quit") {
                    std::cout << "exit program immediately\n";
                } else if (keyword == "select") {
                    std::cout << "syntax : select field_list from mask [ where clause]\n"
                              << "selectivity show information based on citeria\n";
                } else {
                    std::cout << "unknown keyword `" << keyword << "`\n";
                }
            } else {
                std::cout << "Help : \n"
                          << " 1. com mask\n"
                          << " 2. geom mask\n"
                          << " 3. mask\n"
                          << " 4. eda mask1 mask2\n"
                          << " 5. bond mask1 mask2\n"
                          << " 6. angle mask1 mask2 mask3\n"
                          << " 7. dihedral mask1 mask2 mask3 mask4\n"
                          << " 8. energy mask\n"
                          << " 9. select field1,.. from mask [where clause] [ order by field [ASC|DESC] ]\n"
                          << "10. help [keyword]\n"
                          << "11. quit\n";
            }
            continue;
        }

        std::cout << boost::format("%-6s %-7s %4s %-7s %4s %-6s  %8s %8s  %8s%8s%8s") % "#Atom" % "Name" % "#Res" %
                         "Name" % "#Mol" % "Type" % "Charge" % "Mass" % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
        if (frame->atom_list.front()->lj_param.has_value())
            std::cout << boost::format("%14s%14s%14s%14s") % "LJ c6" % "LJ c12" % "sigma(nm)" % "eps(kJ/mol)";
        std::cout << '\n';

        double weight = 0;
        std::tuple<double, double, double> sum{};
        const boost::format fmt{"%6d %-7s %4s %-7s %4s %-6s  %8s %8s  %8.3f%8.3f%8.3f"};
        if (selected_atoms.empty())
            selected_atoms = PBCUtils::find_atoms(ast, frame);

        for (auto &atom : selected_atoms) {
            std::cout << boost::format(fmt) % atom->seq % atom->atom_name %
                             (atom->residue_num ? std::to_string(atom->residue_num.get()) : "-") %
                             (atom->residue_name ? atom->residue_name.get() : "-") % atom->molecule.lock()->sequence %
                             atom->type_name %
                             (atom->charge ? (boost::format("%8.4f") % atom->charge.get()).str() : "-") %
                             (atom->mass ? (boost::format("%8.4f") % atom->mass.get()).str() : "-") % atom->x %
                             atom->y % atom->z;
            if (atom->lj_param.has_value()) {
                const auto c6 = (*atom->lj_param).c6;
                const auto c12 = (*atom->lj_param).c12;
                const auto &[sigma, epsilon] = atom->get_lj_p();
                std::cout << boost::format("%14e%14e%14e%14e") % c6 % c12 % sigma % epsilon;
            }
            std::cout << '\n';
            switch (mode) {
            case Mode::Mass:
                if (!atom->mass) {
                    std::cerr << "atom mass not available !\n";
                    goto beg;
                }
                sum += atom->getCoordinate() * atom->mass.get();
                weight += atom->mass.get();
                break;
            case Mode::Geom:
                sum += atom->getCoordinate() * atom->mass.get();
                weight++;
                break;
            case Mode::Noop:
                break;
            }
        }
        if (weight != 0.0) {
            switch (mode) {
            case Mode::Mass:
                sum /= weight;
                std::cout << boost::format("Mass Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                std::cout << boost::format("            %8.3f%8.3f%8.3f\n") % std::get<0>(sum) % std::get<1>(sum) %
                                 std::get<2>(sum);
                break;
            case Mode::Geom:
                sum /= weight;
                std::cout << boost::format("Geom Center %8s%8s%8s\n") % "X(Ang)" % "Y(Ang)" % "Z(Ang)";
                std::cout << boost::format("            %8.3f%8.3f%8.3f\n") % std::get<0>(sum) % std::get<1>(sum) %
                                 std::get<2>(sum);
                break;
            case Mode::Noop:
                break;
            }
        }
    }
}
