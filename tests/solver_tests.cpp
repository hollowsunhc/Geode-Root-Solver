#include "geode_root_solver.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Real = geode::HighPrecFloat;
using Poly = geode::Polynomial<Real>;

unsigned precision_from_other_translation_unit();

namespace {

void check(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <class Exception, class Function>
void expect_throw(Function function) {
    try {
        function();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error("Expected an exception.");
}

template <class Coefficient>
void check_roots(const geode::Polynomial<Coefficient>& polynomial, std::vector<Real> expected,
                 const Real& tolerance = Real("1e-30")) {
    const auto roots = geode::RootSolver<>(polynomial).compute_roots(tolerance);
    check(roots.size() == expected.size(), "Wrong number of real roots.");
    check(std::is_sorted(roots.begin(), roots.end()), "Roots are not sorted.");
    for (std::size_t i = 0; i < roots.size(); ++i) {
        check(geode::abs(roots[i] - expected[i]) <= tolerance, "Root error exceeds tolerance.");
    }
}

void examples() {
    check_roots(Poly({0, -1, 0, 1}), {-1, 0, 1});
    const Real sqrt_two("1.41421356237309504880168872420969807856967187537694807317667973799");
    check_roots(Poly({-2, 0, 1}), {-sqrt_two, sqrt_two});
    check_roots(Poly({-6, 11, -6, 1}), {1, 2, 3});
    check_roots(Poly({-5, 2}), {Real("2.5")});
    check_roots(Poly({5}), {});
}

void repeated_roots() {
    check_roots(Poly({1, -2, 1}), {1});
    check_roots(Poly({0, 0, -1, 0, 1}), {-1, 0, 1});
    check_roots(Poly({4, -4, -3, 2, 1}), {-2, 1}); // (x + 2)^2 (x - 1)^2
}

void no_real_roots() {
    check_roots(Poly({1, 0, 1}), {});
    check_roots(Poly({1, 1, 1}), {});
    check_roots(Poly({1, 0, 2, 0, 1}), {});
}

void exact_coefficients() {
    using RationalPoly = geode::Polynomial<geode::mpq>;
    geode::set_global_bits(128);
    check_roots(RationalPoly({geode::mpq(1, 100), geode::mpq(-1, 5), 1}), {Real("0.1")});

    const mpz_class denominator("1" + std::string(100, '0'));
    const geode::mpq constant = 1 + geode::mpq(1, denominator);
    check_roots(RationalPoly({constant, -2, 1}), {}); // (x - 1)^2 + 1e-100
    check_roots(Poly({Real(constant), -2, 1}), {1});  // Stored MPFR input has rounded to (x - 1)^2.

    check_roots(RationalPoly({geode::mpq(mpz_class(-2), mpz_class(-4)), geode::mpq(-2, 2),
                              geode::mpq(0, 3)}),
                {Real("0.5")});
    check_roots(RationalPoly({5}), {});
    expect_throw<std::invalid_argument>([] { geode::RootSolver<> solver(RationalPoly({0})); });
    expect_throw<std::invalid_argument>([] {
        geode::mpq invalid(1);
        mpz_set_ui(invalid.get_den_mpz_t(), 0);
        geode::RootSolver<> solver(RationalPoly({invalid, 1}));
    });
}

void invalid_input() {
    expect_throw<std::invalid_argument>([] { Real value("not a number"); });
    expect_throw<std::invalid_argument>([] { geode::RootSolver<> solver(Poly({0})); });
    expect_throw<std::invalid_argument>([] { geode::RootSolver<> solver(Poly({Real("nan"), 1})); });
    for (const char* tolerance : {"0", "-1", "nan", "inf"}) {
        expect_throw<std::invalid_argument>(
            [&] { geode::RootSolver<>(Poly({-1, 1})).compute_roots(Real(tolerance)); });
    }
}

void geode_coefficients() {
    const geode::GeodeTable one_variable(1, 4);
    check(one_variable.get({0}) == 1, "G[0] must be 1.");
    check(one_variable.get({1}) == 2, "G[1] must be 2.");
    check(one_variable.get({2}) == 5, "G[2] must be 5.");
    check(one_variable.get({4}) == 42, "G[4] must be 42.");
    const geode::GeodeTable two_variables(2, 4);
    check(two_variables.get({1, 1}) == 16, "G[1,1] must be 16.");
    check(two_variables.get({0, 2}) == 12, "G[0,2] must be 12.");
}

void precision() {
    geode::set_global_bits(512);
    const Real high("1.123456789012345678901234567890123456789");
    geode::set_global_bits(64);
    Real assigned;
    assigned = high;
    check(assigned == high, "Copy assignment lost precision.");
    check(-(-high) == high, "Unary minus lost precision.");
    check(geode::abs(-high) == high, "Absolute value lost precision.");
    check(Real(0) + high == high, "Mixed-precision addition lost digits.");
    check(Real(high.to_string()) != high, "The low precision fixture is not sensitive.");
    geode::set_global_bits(512);
    check(Real(high.to_string()) == high, "Decimal output did not round-trip.");
    check(precision_from_other_translation_unit() == 512,
          "Precision differs between translation units.");
    const Real nan("nan");
    check(!(nan <= high) && !(nan >= high), "NaN must be unordered.");
    geode::set_global_bits(256);
}

Poly from_roots(const std::vector<Real>& roots) {
    std::vector<Real> coefficients{1};
    for (const auto& root : roots) {
        std::vector<Real> product(coefficients.size() + 1);
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            product[i] -= root * coefficients[i];
            product[i + 1] += coefficients[i];
        }
        coefficients = std::move(product);
    }
    return Poly(std::move(coefficients));
}

void clusters_and_scales() {
    const Real gap = geode::pow(Real(2), -40L);
    check_roots(from_roots({1, 1 + gap, 1 + 2 * gap}), {1, 1 + gap, 1 + 2 * gap});
    // A loose tolerance must not merge roots into one result.
    check_roots(from_roots({1, 1 + gap}), {1, 1 + gap}, Real("1e-5"));
    const Real small = geode::pow(Real(2), -80L);
    const Real large = geode::pow(Real(2), 80L);
    check_roots(from_roots({small, large}), {small, large});
    for (const char* scale : {"1e-100", "1e100", "-1e-100"}) {
        check_roots(Poly({-2 * Real(scale), 0, Real(scale)}),
                    {-Real("1.414213562373095048801688724209698078569671875"),
                     Real("1.414213562373095048801688724209698078569671875")});
    }
    // The minimum is positive. A small residual is not proof of a root.
    check_roots(Poly({1 + geode::pow(Real(2), -150L), -2, 1}), {});
}

void high_precision_and_limits() {
    geode::set_global_bits(512);
    const Real reference(
        "1."
        "414213562373095048801688724209698078569671875376948073176679737990732478462107038850387534"
        "3276415727350138462309122970249248360558507372126441214970999358314132");
    check_roots(Poly({-2, 0, 1}), {-reference, reference}, Real("1e-120"));
    geode::set_global_bits(24);
    expect_throw<std::runtime_error>(
        [] { geode::RootSolver<>(Poly({-2, 0, 1})).compute_roots(Real("1e-30")); });
    expect_throw<std::runtime_error>(
        [] { geode::RootSolver<>(Poly({-1, 3})).compute_roots(Real("1e-30")); });
    check_roots(Poly({-2, 1}), {2}); // An exactly representable root needs no extra bits.
    geode::set_global_bits(256);
}

void polynomial_arithmetic() {
    const geode::Polynomial<geode::mpq> polynomial({-6, 11, -6, 1, 0, 0});
    check(polynomial.degree() == 3, "Trailing zero coefficients changed the degree.");
    check(polynomial(geode::mpq(4)) == 6, "Horner evaluation is wrong.");
    check(polynomial.derivative()(geode::mpq(4)) == 11, "The derivative is wrong.");
    check(geode::Polynomial<geode::mpq>({}).coeff() == std::vector<geode::mpq>{0},
          "An empty polynomial must represent zero.");
    check(Poly({5}).derivative()(Real(4)) == 0, "A constant has a nonzero derivative.");
    check_roots(Poly({-5, 2, 0, 0}), {Real("2.5")});
}

void number_validation() {
    expect_throw<std::invalid_argument>([] { Real value(static_cast<const char*>(nullptr)); });
    expect_throw<std::invalid_argument>([] { Real value("1", 1); });
    expect_throw<std::invalid_argument>([] { Real value(std::string("1\0junk", 6)); });
    expect_throw<std::out_of_range>([] { geode::set_global_bits(0); });
    check(Real(geode::mpq(2, 4)) == Real("0.5"), "A rational input was not normalized.");
    check(Real("ff", 16) == 255, "Parsing a nondecimal base failed.");
    check(Real(42u) == 42, "An unsigned integer input changed value.");
    check(Real(1) + 0.5 == Real("1.5"), "An implicit double operand lost its fraction.");
    const unsigned max_unsigned = std::numeric_limits<unsigned>::max();
    check(Real(1) + max_unsigned == Real(max_unsigned) + 1,
          "An implicit unsigned operand changed value.");
    const long max_long = std::numeric_limits<long>::max();
    check(Real(1) + max_long == Real(max_long) + 1, "An implicit long operand changed value.");
    const unsigned long max_ulong = std::numeric_limits<unsigned long>::max();
    check(Real(1) + max_ulong == Real(max_ulong) + 1,
          "An implicit unsigned long operand changed value.");
    check(Real(1) + (-9223372036854775807LL - 1) == Real("-9223372036854775807"),
          "An implicit signed 64-bit operand changed value.");
    check(Real(1) + 18446744073709551615ULL == Real("18446744073709551616"),
          "An implicit unsigned 64-bit operand changed value.");
    check(Real(1) + 0.5L == Real("1.5"), "An implicit long double operand lost its fraction.");
    check(geode::pow(Real(2), -3) == Real("0.125"), "A signed integer exponent failed.");
    check(geode::pow(Real(2), 3u) == 8, "An unsigned integer exponent failed.");
    for (const char* coefficient : {"nan", "inf", "-inf"}) {
        expect_throw<std::invalid_argument>(
            [&] { geode::RootSolver<> solver(Poly({1, Real(coefficient)})); });
    }
    expect_throw<std::invalid_argument>([] { geode::Ball ball(Real(0), Real(-1)); });
    expect_throw<std::invalid_argument>([] { geode::Ball ball(Real(0), Real("nan")); });
    const geode::Ball ball(Real(1), Real("0.5"));
    check(ball.contains(Real("1.5")) && !ball.contains(Real("1.5001")),
          "Ball boundary test failed.");
    check(!ball.contains(Real("nan")), "A ball must not contain NaN.");
    geode::set_global_bits(2);
    const Real midpoint(1);
    const Real value("-0.25");
    geode::set_global_bits(256);
    check(geode::Ball(midpoint, Real("1.3")).contains(value),
          "A precise radius was ignored in the membership test.");
}

void math_identities() {
    check(geode::hyper_catalan({}) == 1, "C[] must be 1.");
    check(geode::hyper_catalan({2}) == 2, "C[2] must be 2.");
    check(geode::hyper_catalan({1, 1}) == 5, "C[1,1] must be 5.");
    check(geode::hyper_catalan({1, 0, 2}) == 45, "C[1,0,2] must be 45.");
    const geode::GeodeTable table(4, 8);
    check(table.get({1, 1, 1, 1}) == 12344, "The published 4D coefficient is wrong.");
    check(table.get({2, 2, 2, 2}) == geode::mpq("2408941884"),
          "The second 4D coefficient is wrong.");
    check(table.get({9, 0, 0, 0}) == 0, "A term outside the table limit must be zero.");
    expect_throw<std::invalid_argument>([&] { table.get({1}); });
    check(geode::GeodeTable(0, 10).get({}) == 1, "The empty Geode index is wrong.");
    // Gessel's identity must improve on the Newton step from x=3/2 for x^2-2.
    const Real next = geode::detail::geode_step({-2, 0, 1}, Real("1.5"));
    check(geode::abs(next * next - 2) < Real("1e-6"),
          "Geode refinement did not improve the estimate.");
    // When t2+t3=0, S=1. Dividing S-1 by that sum would be invalid.
    check(geode::detail::geode_step({1, -1, 1, -1}, Real(0)) == 1,
          "Geode refinement failed when its parameters sum to zero.");
}

void higher_degree() {
    std::vector<Real> roots;
    for (int i = -8; i <= 8; ++i)
        roots.emplace_back(i);
    check_roots(from_roots(roots), roots);
    // Complex conjugate factors must not add real roots.
    check_roots(Poly({-6, 11, -12, 12, -6, 1}), {1, 2, 3}); // (x^2+1)(x-1)(x-2)(x-3)
    check_roots(Poly({-1, -1, 0, 0, 0, 1}),
                {Real("1.16730397826141868425604589985484218072056037152548903914008244928")});
}

void factored_polynomials() {
    std::mt19937 generator(20261010);
    for (unsigned trial = 0; trial < 40; ++trial) {
        std::vector<Real> roots;
        for (unsigned i = 0; i < 3 + trial % 10; ++i) {
            roots.push_back(Real(static_cast<int>(generator() % 25) - 12) / Real(4));
        }
        const Poly polynomial = from_roots(roots);
        std::sort(roots.begin(), roots.end());
        roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
        check_roots(polynomial, roots);
    }
}

} // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"examples", examples},
        {"repeated roots", repeated_roots},
        {"no real roots", no_real_roots},
        {"exact coefficients", exact_coefficients},
        {"invalid input", invalid_input},
        {"Geode coefficients", geode_coefficients},
        {"precision", precision},
        {"clusters and scales", clusters_and_scales},
        {"high precision and limits", high_precision_and_limits},
        {"polynomial arithmetic", polynomial_arithmetic},
        {"number validation", number_validation},
        {"mathematical identities", math_identities},
        {"factored polynomials", factored_polynomials},
        {"higher degree", higher_degree}};
    int failures = 0;
    for (const auto& test : tests) {
        geode::set_global_bits(256);
        try {
            test.second();
            std::cout << "PASS: " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL: " << test.first << ": " << error.what() << '\n';
        }
    }
    mpfr_free_cache();
    return failures == 0 ? 0 : 1;
}
