/**
 * @file series.hpp
 * @brief Exact local root series with bounded work and certified truncation errors.
 */
#pragma once

#include "geode_root_solver.hpp"

namespace geode::browser {

constexpr unsigned series_order = 8; ///< Highest total degree displayed in the browser.
constexpr std::size_t rational_characters = 4096; ///< Maximum printed fraction length.
constexpr std::size_t explained_products = 12;    ///< Worked products per order; sums include all.

/// @brief One exact multi-index product and its unreduced factorial quotient.
struct SeriesProduct {
    std::vector<unsigned> exponents; ///< Tuple m2, ..., md.
    mpz_class numerator;             ///< Factorial of sum k*mk.
    mpz_class denominator;           ///< (1+sum (k-1)*mk)! times product mk!.
    mpz_class coefficient;           ///< Integer hyper-Catalan coefficient C_m.
    mpq value;                       ///< C_m times product tk^mk.
};

/// @brief One total-degree layer of the hyper-Catalan series and its partial sum.
struct SeriesStep {
    unsigned order;         ///< Total degree in the parameters t2, ..., td.
    mpq coefficient;        ///< Exact coefficient of z^order in S(z).
    mpq term;               ///< Exact contribution delta times coefficient.
    HighPrecFloat sum;      ///< Rounded value of a + delta times the partial series.
    HighPrecFloat bound;    ///< Upward-rounded absolute error, including rounding of sum.
    unsigned product_count; ///< Total products in this order, including those not displayed.
    std::vector<SeriesProduct> products; ///< Bounded worked examples, never used for the sum.
};

/// @brief A certified root expansion with an optional exact value, or an explanation of failure.
struct RootSeries {
    std::string kind = "unavailable"; ///< "series" or "unavailable".
    std::string reason; ///< Explanation when a series cannot be provided within the limits.
    std::size_t reduced_degree = 0; ///< Degree after removing repeated factors.
    bool has_exact_value = false;   ///< Whether an exact rational root was identified.
    mpq value;                      ///< Exact root, displayed alongside its expansion.
    mpq center;                     ///< Rational expansion point a.
    mpq scale;                      ///< Rational change of variable delta, with x = a + delta*y.
    mpq radius; ///< Certified analytic radius rho > 1 for the auxiliary variable z.
    std::vector<mpq> parameters;   ///< Exact t2, ..., td of 1-y+sum tk*y^k.
    std::vector<SeriesStep> steps; ///< Partial sums through total degree series_order.
};

/**
 * @brief Bound fraction growth before serializing or continuing exact arithmetic.
 * @throws std::length_error If the numerator and denominator exceed the display budget.
 */
inline void check_fraction(const mpq& value) {
    if (mpz_sizeinbase(value.get_num_mpz_t(), 10) + mpz_sizeinbase(value.get_den_mpz_t(), 10) + 2 >
        rational_characters)
        throw std::length_error("The exact series exceeds the browser's fraction-size limit.");
}

/// @brief Count Sturm variations at negative infinity using each leading coefficient.
inline int left_variations(const std::vector<detail::RationalPoly>& sequence) {
    int previous = 0;
    int count = 0;
    for (const auto& polynomial : sequence) {
        int sign = mpq_sgn(polynomial.back().get_mpq_t());
        if (polynomial.size() % 2 == 0)
            sign = -sign;
        if (previous != 0 && previous != sign)
            ++count;
        previous = sign;
    }
    return count;
}

/// @brief Shift a polynomial to a rational center without coefficient rounding.
inline detail::RationalPoly shift(detail::RationalPoly p, const mpq& center) {
    for (std::size_t i = p.size() - 1; i > 0; --i) {
        for (std::size_t j = i - 1; j + 1 < p.size(); ++j) {
            p[j] += center * p[j + 1];
            check_fraction(p[j]);
        }
    }
    return p;
}

/// @brief Count weak compositions without enumerating them (at most C(22,8)).
inline unsigned product_count(std::size_t variables, unsigned order) {
    if (variables == 0)
        return order == 0 ? 1 : 0;
    unsigned count = 1;
    for (unsigned i = 1; i <= order; ++i)
        count = count * (static_cast<unsigned>(variables) + i - 1) / i;
    return count;
}

/**
 * @brief Enumerate only the first twelve products of one order for the walkthrough.
 * @note Exact grouped sums are computed separately, including every multi-index.
 *       Recursion stops at the display cap; it never traverses the omitted tuples.
 */
inline void explain_products(const std::vector<mpq>& parameters, unsigned remaining,
                             std::size_t position, std::vector<unsigned>& exponents,
                             std::vector<SeriesProduct>& products) {
    if (products.size() == explained_products)
        return;
    if (position == parameters.size()) {
        if (remaining != 0)
            return;
        SeriesProduct product;
        product.exponents = exponents;
        product.coefficient = hyper_catalan(exponents).get_num();
        unsigned weighted = 0;
        product.value = product.coefficient;
        for (std::size_t i = 0; i < exponents.size(); ++i) {
            weighted += static_cast<unsigned>(i + 2) * exponents[i];
            for (unsigned power = 0; power < exponents[i]; ++power) {
                product.value *= parameters[i];
                check_fraction(product.value);
            }
        }
        mpz_fac_ui(product.numerator.get_mpz_t(), weighted);
        product.denominator = product.numerator / product.coefficient;
        products.push_back(std::move(product));
        return;
    }
    for (unsigned count = remaining + 1; count > 0 && products.size() < explained_products;) {
        exponents[position] = --count;
        explain_products(parameters, remaining - count, position + 1, exponents, products);
    }
}

/**
 * @brief Group hyper-Catalan monomials by total degree using exact convolutions.
 * @param report Certified normalization and analytic radius, filled with nine partial sums.
 * @note S(z)=1+z*sum tk*S(z)^k gives s_j=sum tk*[z^(j-1)]S(z)^k.
 *       Cauchy's bound on |z|=rho gives a tail <= 2*rho^(-N)/(rho-1).
 *       The displayed error also includes the exact rounding error of the partial sum.
 */
inline void expand(RootSeries& report) {
    const std::size_t degree = report.parameters.size() + 1;
    std::vector<std::vector<mpq>> powers(degree + 1, std::vector<mpq>(series_order + 1));
    for (auto& power : powers)
        power[0] = 1;
    mpq partial = report.center;
    mpq tail = report.parameters.empty() ? mpq(0) : 2 * abs(report.scale) / (report.radius - 1);
    for (unsigned order = 0; order <= series_order; ++order) {
        mpq coefficient = order == 0 ? 1 : 0;
        if (order != 0) {
            for (std::size_t k = 2; k <= degree; ++k) {
                coefficient += report.parameters[k - 2] * powers[k][order - 1];
                check_fraction(coefficient);
            }
            powers[1][order] = coefficient;
            for (std::size_t k = 2; k <= degree; ++k) {
                for (unsigned i = 0; i <= order; ++i) {
                    powers[k][order] += powers[1][i] * powers[k - 1][order - i];
                    check_fraction(powers[k][order]);
                }
            }
            tail /= report.radius;
        }
        const mpq term = report.scale * coefficient;
        check_fraction(term);
        partial += term;
        check_fraction(partial);
        const HighPrecFloat sum(partial);
        const mpq error = tail + abs(detail::rational(sum) - partial);
        HighPrecFloat bound;
        mpfr_set_q(bound.get_mpfr_ptr(), error.get_mpq_t(), MPFR_RNDU);
        if (!bound.is_finite())
            throw std::length_error("The series error bound exceeds the browser's numeric range.");
        SeriesStep step{order, coefficient, term,
                        sum,   bound,       product_count(report.parameters.size(), order),
                        {}};
        std::vector<unsigned> exponents(report.parameters.size());
        explain_products(report.parameters, order, 0, exponents, step.products);
        report.steps.push_back(std::move(step));
    }
    report.kind = "series";
}

/**
 * @brief Construct a certified series for one indexed root of a square-free polynomial.
 * @param p Nonconstant polynomial with exact canonical coefficients and no repeated factors.
 * @param sequence Sturm sequence of p, empty only for a linear polynomial.
 * @param root Computed root used to propose expansion centers, never to certify them.
 * @param index Zero-based index in the ascending list of distinct real roots.
 * @return A certified expansion (including identified exact roots), or an unavailable status.
 *
 * For B=sum k*|tk|*2^(k-1) < 1/4 and rho=1/(2B), the map
 * y -> 1+z*sum tk*y^k is a contraction on |y-1|<=1 for |z|<=rho.
 * Its analytic fixed point has |S(z)|<=2. Exact Sturm counts identify its root
 * branch before any partial sums are reported. At most seven precision levels
 * are tried, with one extra shift when a center is already an exact root.
 */
inline RootSeries root_series(const detail::RationalPoly& p,
                              const std::vector<detail::RationalPoly>& sequence,
                              const HighPrecFloat& root, std::size_t index) {
    RootSeries report;
    report.reduced_degree = p.size() - 1;
    report.reason = "A convergent local series could not be certified within the browser limits.";
    try {
        if (p.size() == 2) {
            report.value = -p[0] / p[1];
            check_fraction(report.value);
            report.has_exact_value = true;
            report.center = report.value - 1;
            check_fraction(report.center);
            report.scale = 1;
            report.radius = 2; // S(z)=1 is entire; any finite radius above one certifies it.
            expand(report);
            return report;
        }
        const int left = left_variations(sequence);
        if (root == 0 && p[0] == 0 &&
            left - detail::variations(sequence, mpq(0)) == static_cast<int>(index + 1)) {
            report.has_exact_value = true;
        }
        const auto bits = static_cast<unsigned>(mpfr_get_prec(root.get_mpfr_srcptr()));
        for (unsigned precision = 8; precision <= bits; precision *= 2) {
            HighPrecFloat coarse(root);
            mpfr_prec_round(coarse.get_mpfr_ptr(), precision, MPFR_RNDD);
            const auto exponent = (coarse == 0 ? 0 : mpfr_get_exp(coarse.get_mpfr_srcptr())) -
                                  static_cast<mpfr_exp_t>(precision) - 1;
            mpq offset(1);
            if (exponent < 0)
                mpq_div_2exp(offset.get_mpq_t(), offset.get_mpq_t(),
                             static_cast<mp_bitcnt_t>(-exponent));
            else
                mpq_mul_2exp(offset.get_mpq_t(), offset.get_mpq_t(),
                             static_cast<mp_bitcnt_t>(exponent));
            report.center = detail::rational(coarse) + offset;
            check_fraction(report.center);
            auto shifted = shift(p, report.center);
            if (shifted[0] == 0) {
                if (left - detail::variations(sequence, report.center) ==
                    static_cast<int>(index + 1)) {
                    report.value = report.center;
                    report.has_exact_value = true;
                }
                report.center += offset;
                check_fraction(report.center);
                shifted = shift(p, report.center);
                if (shifted[0] == 0)
                    continue;
            }
            if (shifted[1] == 0)
                continue;
            report.scale = -shifted[0] / shifted[1];
            check_fraction(report.scale);
            report.parameters.clear();
            mpq power = report.scale;
            mpq derivative_bound = 0;
            for (std::size_t k = 2; k < p.size(); ++k) {
                const mpq parameter = -shifted[k] * power / shifted[1];
                check_fraction(parameter);
                report.parameters.push_back(parameter);
                derivative_bound +=
                    static_cast<unsigned long>(k) * abs(parameter) * (1UL << (k - 1));
                check_fraction(derivative_bound);
                power *= report.scale;
                check_fraction(power);
            }
            if (derivative_bound <= 0 || derivative_bound >= mpq(1, 4))
                continue;
            const mpq end = report.center + 2 * report.scale;
            const mpq lower = std::min(report.center, end);
            const mpq upper = std::max(report.center, end);
            const int lower_count = detail::variations(sequence, lower);
            if (left - lower_count != static_cast<int>(index) ||
                lower_count - detail::variations(sequence, upper) != 1)
                continue;
            report.radius = 1 / (2 * derivative_bound);
            check_fraction(report.radius);
            expand(report);
            return report;
        }
    } catch (const std::length_error& error) {
        report.reason = error.what();
        report.steps.clear();
    }
    return report;
}

} // namespace geode::browser
