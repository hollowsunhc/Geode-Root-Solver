/**
 * @file geode_root_solver.hpp
 * @brief Real polynomial roots with exact isolation and Geode refinement.
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

/// @brief Polynomial arithmetic and real root solving with GMP and MPFR.
namespace geode {

/// @brief Exact rational number type.
using mpq = mpq_class;

/**
 * @brief An MPFR real number with automatic storage management.
 *
 * New values use the current thread's default precision (initially 256 bits).
 * Copies retain the source precision. Binary arithmetic uses the greater operand
 * precision and rounds to nearest, with ties to even. MPFR NaN and infinity are
 * supported here, but RootSolver rejects them as coefficients and tolerances.
 */
class MpfrFloat {
public:
    mpfr_t m_val; ///< MPFR storage; callers must preserve its initialized state.

    /**
     * @brief Set the default precision for new values in this thread.
     * @param bits Number of binary significand bits; existing values do not change.
     * @throws std::out_of_range If bits is outside the MPFR precision range.
     */
    static void set_global_precision_bits(unsigned bits) {
        if (bits < MPFR_PREC_MIN || bits > static_cast<unsigned long>(MPFR_PREC_MAX)) {
            throw std::out_of_range("MPFR precision is out of range.");
        }
        default_precision_ = bits;
    }

    /// @brief Get the default precision, in bits, for this thread.
    static unsigned get_global_precision_bits() { return default_precision_; }

    /// @brief Construct positive zero at the current default precision.
    MpfrFloat() {
        mpfr_init2(m_val, default_precision_);
        mpfr_set_zero(m_val, 1);
    }
    /// @brief Copy a double; digits already lost in the double cannot be restored.
    MpfrFloat(double value) : MpfrFloat() { mpfr_set_d(m_val, value, MPFR_RNDN); }
    /// @brief Convert a long double at the current default precision.
    MpfrFloat(long double value) : MpfrFloat() { mpfr_set_ld(m_val, value, MPFR_RNDN); }
    /// @brief Convert a signed long at the current default precision.
    MpfrFloat(long value) : MpfrFloat() { mpfr_set_si(m_val, value, MPFR_RNDN); }
    /// @brief Convert an unsigned long at the current default precision.
    MpfrFloat(unsigned long value) : MpfrFloat() { mpfr_set_ui(m_val, value, MPFR_RNDN); }
    /// @brief Convert an integer, including implicit integer operands.
    MpfrFloat(int value) : MpfrFloat(static_cast<long>(value)) {}
    /// @brief Convert an unsigned integer at the current default precision.
    MpfrFloat(unsigned value) : MpfrFloat(static_cast<unsigned long>(value)) {}
    /// @brief Convert a signed long long at the current default precision.
    /// @note Decimal text avoids narrowing when long has only 32 bits.
    MpfrFloat(long long value) : MpfrFloat(std::to_string(value)) {}
    /// @brief Convert an unsigned long long at the current default precision.
    MpfrFloat(unsigned long long value) : MpfrFloat(std::to_string(value)) {}

    /**
     * @brief Parse an MPFR number string.
     * @param value Non-null, null-terminated number string.
     * @param base Base from 2 through 62, or 0 for MPFR base detection.
     * @throws std::invalid_argument If the string or base is invalid.
     */
    explicit MpfrFloat(const char* value, int base = 10) : MpfrFloat() {
        if (value == nullptr || (base != 0 && (base < 2 || base > 62)) ||
            mpfr_set_str(m_val, value, base, MPFR_RNDN) != 0) {
            throw std::invalid_argument("Invalid number or base.");
        }
    }
    /// @brief Parse a number string; embedded null characters are invalid.
    /// @copydetails MpfrFloat(const char*, int)
    explicit MpfrFloat(const std::string& value, int base = 10) : MpfrFloat(value.c_str(), base) {
        if (value.find('\0') != std::string::npos) {
            throw std::invalid_argument("A number must not contain a null character.");
        }
    }
    /**
     * @brief Normalize and round an exact rational to the default precision.
     * @throws std::invalid_argument If the denominator is zero.
     */
    explicit MpfrFloat(mpq value) : MpfrFloat() {
        if (value.get_den() == 0) {
            throw std::invalid_argument("A rational denominator must not be zero.");
        }
        value.canonicalize();
        mpfr_set_q(m_val, value.get_mpq_t(), MPFR_RNDN);
    }

    /// @brief Copy both value and precision.
    MpfrFloat(const MpfrFloat& other) {
        mpfr_init2(m_val, mpfr_get_prec(other.m_val));
        mpfr_set(m_val, other.m_val, MPFR_RNDN);
    }
    /// @brief Replace both value and precision with those of the source.
    MpfrFloat& operator=(const MpfrFloat& other) {
        if (this != &other) {
            mpfr_set_prec(m_val, mpfr_get_prec(other.m_val));
            mpfr_set(m_val, other.m_val, MPFR_RNDN);
        }
        return *this;
    }
    /// @brief Transfer storage; the source remains valid with an unspecified value.
    MpfrFloat(MpfrFloat&& other) noexcept {
        mpfr_init2(m_val, MPFR_PREC_MIN);
        mpfr_swap(m_val, other.m_val);
    }
    /// @brief Exchange storage with the source without allocation.
    MpfrFloat& operator=(MpfrFloat&& other) noexcept {
        mpfr_swap(m_val, other.m_val);
        return *this;
    }
    /// @brief Release the MPFR significand.
    ~MpfrFloat() { mpfr_clear(m_val); }

    /// @brief Access writable MPFR storage without transferring ownership.
    mpfr_ptr get_mpfr_ptr() { return m_val; }
    /// @brief Access read-only MPFR storage without transferring ownership.
    mpfr_srcptr get_mpfr_srcptr() const { return m_val; }
    /// @brief Return true for a finite value, including signed zero.
    bool is_finite() const { return mpfr_number_p(m_val) != 0; }

    /**
     * @brief Format a decimal string with enough digits to recover the stored value.
     * @param rounding MPFR rounding mode for decimal conversion.
     * @return Scientific notation, signed zero, "nan", "inf", or "-inf".
     * @note Exact recovery uses round-to-nearest conversion and parsing at the original precision.
     */
    std::string to_string(mpfr_rnd_t rounding = MPFR_RNDN) const {
        if (mpfr_nan_p(m_val))
            return "nan";
        if (mpfr_inf_p(m_val))
            return mpfr_signbit(m_val) ? "-inf" : "inf";
        if (mpfr_zero_p(m_val))
            return mpfr_signbit(m_val) ? "-0" : "0";
        mpfr_exp_t exponent;
        std::unique_ptr<char, decltype(&mpfr_free_str)> digits(
            mpfr_get_str(nullptr, &exponent, 10, 0, m_val, rounding), mpfr_free_str);
        if (!digits)
            throw std::bad_alloc();
        const bool negative = digits.get()[0] == '-';
        return std::string(negative ? "-0." : "0.") + (digits.get() + negative) + "e" +
               std::to_string(exponent);
    }

    /// @brief Negate the value without changing its precision.
    MpfrFloat operator-() const {
        MpfrFloat result(*this);
        mpfr_neg(result.m_val, result.m_val, MPFR_RNDN);
        return result;
    }
    /// @brief Add rhs, using the greater operand precision.
    MpfrFloat& operator+=(const MpfrFloat& rhs) {
        promote(rhs);
        mpfr_add(m_val, m_val, rhs.m_val, MPFR_RNDN);
        return *this;
    }
    /// @brief Subtract rhs, using the greater operand precision.
    MpfrFloat& operator-=(const MpfrFloat& rhs) {
        promote(rhs);
        mpfr_sub(m_val, m_val, rhs.m_val, MPFR_RNDN);
        return *this;
    }
    /// @brief Multiply by rhs, using the greater operand precision.
    MpfrFloat& operator*=(const MpfrFloat& rhs) {
        promote(rhs);
        mpfr_mul(m_val, m_val, rhs.m_val, MPFR_RNDN);
        return *this;
    }
    /// @brief Divide by rhs, using the greater operand precision.
    MpfrFloat& operator/=(const MpfrFloat& rhs) {
        promote(rhs);
        mpfr_div(m_val, m_val, rhs.m_val, MPFR_RNDN);
        return *this;
    }
    /// @brief Return the sum, using the greater operand precision.
    friend MpfrFloat operator+(MpfrFloat lhs, const MpfrFloat& rhs) {
        lhs += rhs;
        return lhs;
    }
    /// @brief Return the difference, using the greater operand precision.
    friend MpfrFloat operator-(MpfrFloat lhs, const MpfrFloat& rhs) {
        lhs -= rhs;
        return lhs;
    }
    /// @brief Return the product, using the greater operand precision.
    friend MpfrFloat operator*(MpfrFloat lhs, const MpfrFloat& rhs) {
        lhs *= rhs;
        return lhs;
    }
    /// @brief Return the quotient, using the greater operand precision.
    friend MpfrFloat operator/(MpfrFloat lhs, const MpfrFloat& rhs) {
        lhs /= rhs;
        return lhs;
    }

    /// @brief Test less than; NaN is unordered.
    friend bool operator<(const MpfrFloat& a, const MpfrFloat& b) {
        return mpfr_less_p(a.m_val, b.m_val) != 0;
    }
    /// @brief Test greater than; NaN is unordered.
    friend bool operator>(const MpfrFloat& a, const MpfrFloat& b) { return b < a; }
    /// @brief Test less than or equal; NaN is unordered.
    friend bool operator<=(const MpfrFloat& a, const MpfrFloat& b) {
        return mpfr_lessequal_p(a.m_val, b.m_val) != 0;
    }
    /// @brief Test greater than or equal; NaN is unordered.
    friend bool operator>=(const MpfrFloat& a, const MpfrFloat& b) { return b <= a; }
    /// @brief Test numeric equality; NaN is unequal to every value.
    friend bool operator==(const MpfrFloat& a, const MpfrFloat& b) {
        return mpfr_equal_p(a.m_val, b.m_val) != 0;
    }
    /// @brief Test numeric inequality; NaN is unequal to every value.
    friend bool operator!=(const MpfrFloat& a, const MpfrFloat& b) { return !(a == b); }

private:
    void promote(const MpfrFloat& other) {
        if (mpfr_get_prec(other.m_val) > mpfr_get_prec(m_val)) {
            mpfr_prec_round(m_val, mpfr_get_prec(other.m_val), MPFR_RNDN);
        }
    }
    inline static thread_local unsigned default_precision_ = 256;
};

/// @brief Floating-point type supported by RootSolver.
using HighPrecFloat = MpfrFloat;
/// @copydoc MpfrFloat::set_global_precision_bits
inline void set_global_bits(unsigned bits) {
    MpfrFloat::set_global_precision_bits(bits);
}
/// @copydoc MpfrFloat::get_global_precision_bits
inline unsigned get_global_bits() {
    return MpfrFloat::get_global_precision_bits();
}

/// @brief Return the absolute value at the input precision.
inline MpfrFloat abs(MpfrFloat value) {
    mpfr_abs(value.m_val, value.m_val, MPFR_RNDN);
    return value;
}
/// @brief Return the cosine of an angle in radians, at the input precision.
inline MpfrFloat cos(MpfrFloat value) {
    mpfr_cos(value.m_val, value.m_val, MPFR_RNDN);
    return value;
}
/// @brief Return the sine of an angle in radians, at the input precision.
inline MpfrFloat sin(MpfrFloat value) {
    mpfr_sin(value.m_val, value.m_val, MPFR_RNDN);
    return value;
}
/// @brief Return the inverse cosine in radians, at the input precision.
inline MpfrFloat acos(MpfrFloat value) {
    mpfr_acos(value.m_val, value.m_val, MPFR_RNDN);
    return value;
}
/// @brief Return the natural logarithm at the input precision.
inline MpfrFloat log(MpfrFloat value) {
    mpfr_log(value.m_val, value.m_val, MPFR_RNDN);
    return value;
}
/// @brief Round toward positive infinity to an integer, at the input precision.
inline MpfrFloat ceil(MpfrFloat value) {
    mpfr_ceil(value.m_val, value.m_val);
    return value;
}
/// @brief Raise a real base to a real power, using the greater operand precision.
/// @note Wrap noninteger and long long exponents in MpfrFloat to select this overload.
inline MpfrFloat pow(MpfrFloat base, const MpfrFloat& exponent) {
    mpfr_prec_round(base.m_val, std::max(mpfr_get_prec(base.m_val), mpfr_get_prec(exponent.m_val)),
                    MPFR_RNDN);
    mpfr_pow(base.m_val, base.m_val, exponent.m_val, MPFR_RNDN);
    return base;
}
/// @brief Raise a base to an unsigned integer power, at the base precision.
inline MpfrFloat pow(MpfrFloat base, unsigned long exponent) {
    mpfr_pow_ui(base.m_val, base.m_val, exponent, MPFR_RNDN);
    return base;
}
/// @brief Raise a base to a signed integer power, at the base precision.
inline MpfrFloat pow(MpfrFloat base, long exponent) {
    mpfr_pow_si(base.m_val, base.m_val, exponent, MPFR_RNDN);
    return base;
}
/// @brief Raise a base to a signed integer power, at the base precision.
inline MpfrFloat pow(MpfrFloat base, int exponent) {
    return pow(std::move(base), static_cast<long>(exponent));
}
/// @brief Raise a base to an unsigned integer power, at the base precision.
inline MpfrFloat pow(MpfrFloat base, unsigned exponent) {
    return pow(std::move(base), static_cast<unsigned long>(exponent));
}
/// @brief Return magnitude * cos(angle), a real projection, with angle in radians.
inline MpfrFloat polar(const MpfrFloat& magnitude, const MpfrFloat& angle) {
    return magnitude * cos(angle);
}
/// @brief Legacy 256-bit pi constant; use acos(HighPrecFloat(-1)) at other precisions.
inline const HighPrecFloat PI = acos(HighPrecFloat(-1));

/// @brief A closed real interval [mid - rad, mid + rad], without interval arithmetic.
struct Ball {
    HighPrecFloat mid; ///< Finite midpoint.
    HighPrecFloat rad; ///< Finite, nonnegative radius.
    /// @brief Construct the interval [0, 0].
    Ball() = default;
    /// @brief Construct an interval from a midpoint and radius.
    /// @throws std::invalid_argument For a nonfinite input or a negative radius.
    Ball(const HighPrecFloat& midpoint, const HighPrecFloat& radius) : mid(midpoint), rad(radius) {
        if (!mid.is_finite() || !rad.is_finite() || rad < 0) {
            throw std::invalid_argument(
                "A ball needs a finite midpoint and a finite, nonnegative radius.");
        }
    }
    /// @brief Test membership, including endpoints; nonfinite values return false.
    bool contains(const HighPrecFloat& value) const {
        if (!value.is_finite() || !mid.is_finite() || !rad.is_finite() || rad < 0)
            return false;
        HighPrecFloat distance(rad);
        mpfr_set_prec(distance.m_val,
                      std::max({mpfr_get_prec(value.m_val), mpfr_get_prec(mid.m_val),
                                mpfr_get_prec(rad.m_val)}));
        mpfr_sub(distance.m_val, (value >= mid ? value : mid).m_val,
                 (value >= mid ? mid : value).m_val, MPFR_RNDU);
        return distance <= rad;
    }
};

/**
 * @brief A dense polynomial with coefficients in ascending order {a0, a1, ..., an}.
 * @tparam R Coefficient type with integer construction and basic arithmetic.
 * @note An empty list represents zero. Trailing zero coefficients are removed.
 */
template <class R>
class Polynomial {
public:
    using value_type = R; ///< Coefficient type.
    /// @brief Store and normalize coefficients in ascending order.
    explicit Polynomial(std::vector<R> coefficients) : coefficients_(std::move(coefficients)) {
        while (coefficients_.size() > 1 && coefficients_.back() == R(0))
            coefficients_.pop_back();
        if (coefficients_.empty())
            coefficients_.push_back(R(0));
    }
    /// @brief Get the degree; zero and constant polynomials have degree zero.
    std::size_t degree() const noexcept { return coefficients_.size() - 1; }
    /// @brief Read the normalized coefficients in ascending order.
    const std::vector<R>& coeff() const noexcept { return coefficients_; }

    /// @brief Evaluate the polynomial at x with Horner arithmetic.
    /// @tparam X Evaluation type that can represent a coefficient.
    template <class X>
    X operator()(const X& x) const {
        X value(coefficients_.back());
        for (std::size_t i = degree(); i > 0; --i) {
            value *= x;
            value += X(coefficients_[i - 1]);
        }
        return value;
    }
    /// @brief Return the formal first derivative; a constant gives zero.
    Polynomial derivative() const {
        if (degree() == 0)
            return Polynomial({R(0)});
        std::vector<R> result(degree());
        for (std::size_t i = 1; i < coefficients_.size(); ++i) {
            result[i - 1] = coefficients_[i] * R(static_cast<unsigned long>(i));
        }
        return Polynomial(std::move(result));
    }

private:
    std::vector<R> coefficients_;
};

/**
 * @brief Compute an exact hyper-Catalan coefficient C[m2, m3, ...].
 * @param index Counts of triangles, quadrilaterals, and subsequent polygon types.
 * @return (sum k*mk)! / ((1 + sum (k-1)*mk)! * product mk!), or 1 for an empty index.
 * @throws std::overflow_error If an index sum exceeds unsigned long capacity.
 * @see https://arxiv.org/abs/2507.04552 (Theorem 3)
 */
inline mpq hyper_catalan(const std::vector<unsigned>& index) {
    unsigned long edges = 0;
    unsigned long faces = 0;
    mpz_class denominator = 1;
    mpz_class factorial;
    for (std::size_t i = 0; i < index.size(); ++i) {
        const auto count = index[i];
        if (i > std::numeric_limits<unsigned long>::max() - 2 ||
            count > (std::numeric_limits<unsigned long>::max() - edges) / (i + 2)) {
            throw std::overflow_error("The coefficient index is too large.");
        }
        edges += static_cast<unsigned long>(i + 2) * count;
        faces += count;
        mpz_fac_ui(factorial.get_mpz_t(), count);
        denominator *= factorial;
    }
    if (faces == 0)
        return 1;
    mpz_fac_ui(factorial.get_mpz_t(), edges - faces + 1);
    denominator *= factorial;
    mpz_fac_ui(factorial.get_mpz_t(), edges);
    mpq result(factorial, denominator);
    result.canonicalize();
    return result;
}

/**
 * @brief Calculate and cache exact Geode coefficients G[m2, m3, ...] on demand.
 * @note G differs from the hyper-Catalan array C; S = 1 + (sum tk) G.
 * @warning Concurrent access to one table requires an external lock.
 * @see https://arxiv.org/abs/2507.04552 (Geode recurrence)
 */
class GeodeTable {
public:
    /// @brief Create an empty cache with the given index dimension and sum limit.
    /// @param variables Number of index components, starting with m2.
    /// @param max_faces Largest permitted sum of index components.
    GeodeTable(unsigned variables, unsigned max_faces)
        : variables_(variables), max_faces_(max_faces) {}

    /**
     * @brief Read a coefficient, calculating missing cache entries as needed.
     * @param index Vector of exactly the dimension supplied to the constructor.
     * @return G[index], or zero if the index sum exceeds max_faces.
     * @throws std::invalid_argument If the index has the wrong dimension.
     * @throws std::overflow_error If an intermediate index exceeds its capacity.
     */
    mpq get(const std::vector<unsigned>& index) const {
        if (index.size() != variables_)
            throw std::invalid_argument("Wrong Geode index dimension.");
        unsigned remaining = max_faces_;
        for (auto count : index) {
            if (count > remaining)
                return 0;
            remaining -= count;
        }
        if (index.empty())
            return 1;
        return coefficient(index);
    }

private:
    struct IndexHash {
        std::size_t operator()(const std::vector<unsigned>& index) const noexcept {
            std::size_t hash = 0;
            for (auto value : index)
                hash = hash * 1315423911u ^ value;
            return hash;
        }
    };
    mpq coefficient(const std::vector<unsigned>& index) const {
        const auto found = cache_.find(index);
        if (found != cache_.end())
            return found->second;
        auto next = index;
        const auto pivot =
            static_cast<std::size_t>(std::max_element(next.begin(), next.end()) - next.begin());
        if (next[pivot] == std::numeric_limits<unsigned>::max()) {
            throw std::overflow_error("The Geode index is too large.");
        }
        ++next[pivot];
        mpq result = hyper_catalan(next);
        // Rubine (2025): C[m + e_p] = G[m] + sum_{j != p} G[m + e_p - e_j].
        for (std::size_t j = 0; j < next.size(); ++j) {
            if (j != pivot && next[j] != 0) {
                --next[j];
                result -= coefficient(next);
                ++next[j];
            }
        }
        cache_.emplace(index, result);
        return result;
    }
    unsigned variables_;
    unsigned max_faces_;
    mutable std::unordered_map<std::vector<unsigned>, mpq, IndexHash> cache_;
};

/// @cond INTERNAL
namespace detail {

using RationalPoly = std::vector<mpq>;

/// @brief Convert a finite MPFR value to its exact binary rational value.
inline mpq rational(const HighPrecFloat& value) {
    if (!value.is_finite())
        throw std::invalid_argument("A coefficient must be finite.");
    mpq result;
    mpfr_get_q(result.get_mpq_t(), value.m_val);
    return result;
}

/// @brief Remove trailing zeros from a nonempty coefficient vector.
inline void trim(RationalPoly& p) {
    while (p.size() > 1 && p.back() == 0)
        p.pop_back();
}

/// @brief Evaluate a nonempty rational polynomial without rounding.
inline mpq evaluate(const RationalPoly& p, const mpq& x) {
    mpq value = p.back();
    for (std::size_t i = p.size() - 1; i > 0; --i)
        value = value * x + p[i - 1];
    return value;
}

/// @brief Scale to primitive integer coefficients, preserving every Sturm sign.
inline void primitive(RationalPoly& p) {
    mpz_class denominator = 1;
    for (const auto& value : p)
        mpz_lcm(denominator.get_mpz_t(), denominator.get_mpz_t(), value.get_den_mpz_t());
    mpz_class divisor = 0;
    for (auto& value : p) {
        value *= denominator;
        mpz_gcd(divisor.get_mpz_t(), divisor.get_mpz_t(), value.get_num_mpz_t());
    }
    if (divisor != 0)
        for (auto& value : p)
            value /= divisor;
}

/// @brief Compute the exact formal derivative of a nonempty polynomial.
inline RationalPoly derivative(const RationalPoly& p) {
    if (p.size() == 1)
        return {0};
    RationalPoly result(p.size() - 1);
    for (std::size_t i = 1; i < p.size(); ++i)
        result[i - 1] = p[i] * static_cast<unsigned long>(i);
    return result;
}

/// @brief Return the quotient and replace the dividend with its exact remainder.
/// @pre divisor is nonzero and both polynomials have no trailing zero coefficients.
inline RationalPoly divide(RationalPoly& remainder, const RationalPoly& divisor) {
    RationalPoly quotient(remainder.size() >= divisor.size() ? remainder.size() - divisor.size() + 1
                                                             : 1);
    while (remainder.size() >= divisor.size() && !(remainder.size() == 1 && remainder[0] == 0)) {
        const auto shift = remainder.size() - divisor.size();
        const mpq factor = remainder.back() / divisor.back();
        quotient[shift] = factor;
        for (std::size_t i = 0; i < divisor.size(); ++i)
            remainder[i + shift] -= factor * divisor[i];
        trim(remainder);
    }
    return quotient;
}

/// @brief Build a Sturm sequence whose last member is a scaled polynomial GCD.
/// @pre p has positive degree and no trailing zero coefficients.
inline std::vector<RationalPoly> sturm_sequence(const RationalPoly& p) {
    std::vector<RationalPoly> sequence{p, derivative(p)};
    primitive(sequence.back());
    while (sequence.back().size() > 1) {
        auto remainder = sequence[sequence.size() - 2];
        divide(remainder, sequence.back());
        if (remainder.size() == 1 && remainder[0] == 0)
            break;
        for (auto& value : remainder)
            value = -value;
        primitive(remainder);
        sequence.push_back(std::move(remainder));
    }
    return sequence;
}

/// @brief Reduce p to its primitive square-free part and return its Sturm sequence.
/// @pre p is nonempty and its coefficients are canonical rationals.
/// @return An empty sequence for a constant or linear result.
inline std::vector<RationalPoly> square_free_sequence(RationalPoly& p) {
    trim(p);
    primitive(p);
    if (p.size() <= 2)
        return {};
    auto sequence = sturm_sequence(p);
    if (sequence.back().size() > 1) {
        auto remainder = p;
        p = divide(remainder, sequence.back());
        primitive(p);
        if (p.size() <= 2)
            return {};
        sequence = sturm_sequence(p);
    }
    return sequence;
}

/// @brief Count sign changes at x, omitting zero members as Sturm's theorem requires.
inline int variations(const std::vector<RationalPoly>& sequence, const mpq& x) {
    int previous = 0;
    int count = 0;
    for (const auto& p : sequence) {
        const mpq value = evaluate(p, x);
        const int sign = mpq_sgn(value.get_mpq_t());
        if (sign != 0) {
            if (previous != 0 && previous != sign)
                ++count;
            previous = sign;
        }
    }
    return count;
}

/**
 * @brief Apply Gessel's Geode identity twice near a Newton estimate.
 * @pre coefficients represents a nonconstant polynomial in ascending order.
 * @note This only proposes a point; refine() checks it with exact arithmetic.
 * @see https://arxiv.org/abs/2507.09405 (Theorem 2.1)
 */
inline HighPrecFloat geode_step(const std::vector<HighPrecFloat>& coefficients,
                                const HighPrecFloat& x) {
    auto shifted = coefficients;
    for (std::size_t i = shifted.size() - 1; i > 0; --i) {
        for (std::size_t j = i - 1; j + 1 < shifted.size(); ++j) {
            shifted[j] += x * shifted[j + 1];
        }
    }
    if (shifted[1] == 0)
        return x;
    const HighPrecFloat delta = -shifted[0] / shifted[1];
    HighPrecFloat power(delta);
    HighPrecFloat total(0);
    std::vector<HighPrecFloat> parameters;
    for (std::size_t k = 2; k < shifted.size(); ++k) {
        parameters.push_back(-shifted[k] * power / shifted[1]);
        total += parameters.back();
        power *= delta;
    }
    HighPrecFloat s(1);
    for (unsigned iteration = 0; iteration < 2; ++iteration) {
        HighPrecFloat geometric_sum(1);
        HighPrecFloat denominator(1);
        for (const auto& parameter : parameters) {
            geometric_sum = 1 + s * geometric_sum;
            denominator -= parameter * geometric_sum;
        }
        if (!denominator.is_finite() || abs(denominator) < HighPrecFloat("0.5"))
            return x + delta;
        s = 1 + total / denominator;
        if (!s.is_finite() || abs(s - 1) > HighPrecFloat("0.5"))
            return x + delta;
    }
    return x + delta * s;
}

/// @brief Round an exact root and check its absolute error with rational arithmetic.
inline HighPrecFloat rounded_root(const mpq& root, const mpq& tolerance) {
    HighPrecFloat result(root);
    if (!result.is_finite() || abs(rational(result) - root) > tolerance) {
        throw std::runtime_error("Root accuracy needs more bits. Increase the precision.");
    }
    return result;
}

/**
 * @brief Refine an isolated simple root with checked Geode steps and bisection.
 * @pre The endpoints are not roots, and the open interval contains exactly one root.
 * @pre tolerance is positive, and p has no repeated factors.
 * @note Each successful return has an exact sign bracket within tolerance of the result.
 */
inline HighPrecFloat refine(const RationalPoly& p, mpq lower, mpq upper, const mpq& tolerance) {
    const mpq lower_value = evaluate(p, lower);
    const int lower_sign = mpq_sgn(lower_value.get_mpq_t());
    std::vector<HighPrecFloat> coefficients;
    coefficients.reserve(p.size());
    for (const auto& coefficient : p)
        coefficients.emplace_back(coefficient);
    HighPrecFloat x(mpq((lower + upper) / 2));
    unsigned slow_steps = 0;
    while (true) {
        mpq position = rational(x);
        if (position <= lower || position >= upper) {
            x = HighPrecFloat(mpq((lower + upper) / 2));
            position = rational(x);
            if (position <= lower || position >= upper) {
                if (std::max(mpq(abs(position - lower)), mpq(abs(position - upper))) <= tolerance)
                    return x;
                throw std::runtime_error("Root accuracy needs more bits. Increase the precision.");
            }
        }
        const mpq value = evaluate(p, position);
        const int sign = mpq_sgn(value.get_mpq_t());
        if (sign == 0)
            return x;
        const mpq old_width = upper - lower;
        if (sign == lower_sign) {
            const mpq probe = std::min(upper, mpq(position + tolerance));
            const mpq probe_value = evaluate(p, probe);
            if (mpq_sgn(probe_value.get_mpq_t()) != lower_sign)
                return x;
            lower = probe;
        } else {
            const mpq probe = std::max(lower, mpq(position - tolerance));
            const mpq probe_value = evaluate(p, probe);
            if (mpq_sgn(probe_value.get_mpq_t()) != sign)
                return x;
            upper = probe;
        }
        slow_steps = (upper - lower > old_width / 2) ? slow_steps + 1 : 0;
        HighPrecFloat next = geode_step(coefficients, x);
        if (slow_steps >= 2 || !next.is_finite() || next == x || rational(next) <= lower ||
            rational(next) >= upper) {
            next = HighPrecFloat(mpq((lower + upper) / 2));
            slow_steps = 0;
        }
        x = std::move(next);
    }
}

} // namespace detail
/// @endcond

/**
 * @brief Find every distinct real root of the stored polynomial.
 * @tparam Float Must be HighPrecFloat; the parameter is retained for source compatibility.
 *
 * Exact square-free reduction and Sturm counts isolate the roots. Geode corrections
 * accelerate refinement; exact sign tests certify the absolute error. Nonreal roots
 * are omitted. Repeated roots appear once. No result is accepted by residual alone.
 */
template <class Float = HighPrecFloat>
class RootSolver {
    static_assert(std::is_same_v<Float, HighPrecFloat>,
                  "RootSolver requires geode::HighPrecFloat.");

public:
    using PolyF = Polynomial<Float>; ///< MPFR input polynomial type.
    /**
     * @brief Copy the stored binary coefficients as exact rationals.
     * @param polynomial Nonzero polynomial with finite real coefficients.
     * @throws std::invalid_argument For a zero polynomial or nonfinite coefficient.
     * @note Set precision before constructing coefficients; earlier rounding is permanent.
     */
    explicit RootSolver(const PolyF& polynomial) {
        for (const auto& coefficient : polynomial.coeff())
            coefficients_.push_back(detail::rational(coefficient));
        normalize();
    }

    /**
     * @brief Copy exact rational coefficients without rounding to MPFR precision.
     * @param polynomial Nonzero polynomial; noncanonical rationals are normalized.
     * @throws std::invalid_argument For a zero polynomial or a zero denominator.
     * @note Precision affects the returned roots, never the input coefficients.
     */
    explicit RootSolver(const Polynomial<mpq>& polynomial) : coefficients_(polynomial.coeff()) {
        normalize();
    }

    /**
     * @brief Compute distinct real roots in ascending order.
     * @param tolerance Finite, positive absolute root error; defaults to 1e-30.
     * @return One value per distinct real root, or an empty vector if none exist.
     * @throws std::invalid_argument If tolerance is not finite and positive.
     * @throws std::runtime_error If the current precision cannot meet the error bound
     *         or distinguish two roots. Increase precision and try again.
     * @note Accuracy applies to exact rational inputs or the stored binary values of MPFR inputs.
     */
    std::vector<Float> compute_roots(Float tolerance = Float("1e-30")) const {
        if (!tolerance.is_finite() || tolerance <= 0) {
            throw std::invalid_argument("Tolerance must be finite and positive.");
        }
        const mpq error = detail::rational(tolerance);
        auto p = coefficients_;
        const auto sequence = detail::square_free_sequence(p);
        if (p.size() == 1)
            return {};
        if (p.size() == 2)
            return {detail::rounded_root(mpq(-p[0] / p[1]), error)};
        mpq bound = 0;
        for (std::size_t i = 0; i + 1 < p.size(); ++i)
            bound = std::max(bound, mpq(abs(p[i] / p.back())));
        mpz_class integer_bound;
        mpz_cdiv_q(integer_bound.get_mpz_t(), bound.get_num_mpz_t(), bound.get_den_mpz_t());
        bound = integer_bound + 1;
        struct Interval {
            mpq lower, upper;
            int left, right;
        };
        std::vector<Interval> pending{{-bound, bound, detail::variations(sequence, -bound),
                                       detail::variations(sequence, bound)}};
        std::vector<Float> roots;
        while (!pending.empty()) {
            auto interval = std::move(pending.back());
            pending.pop_back();
            const int count = interval.left - interval.right;
            if (count == 0)
                continue;
            if (count == 1) {
                if (detail::evaluate(p, interval.upper) == 0) {
                    roots.push_back(detail::rounded_root(interval.upper, error));
                    continue;
                }
                if (detail::evaluate(p, interval.lower) != 0) {
                    roots.push_back(detail::refine(p, interval.lower, interval.upper, error));
                    continue;
                }
            }
            // Sturm counts roots in (lower, upper], including a root at the split only on the left.
            const mpq midpoint = (interval.lower + interval.upper) / 2;
            const int middle = detail::variations(sequence, midpoint);
            pending.push_back({midpoint, interval.upper, middle, interval.right});
            pending.push_back({interval.lower, midpoint, interval.left, middle});
        }
        for (std::size_t i = 1; i < roots.size(); ++i) {
            if (roots[i] <= roots[i - 1]) {
                throw std::runtime_error("Distinct roots need more bits. Increase the precision.");
            }
        }
        return roots;
    }

private:
    /// @brief Validate and normalize the copied input before exact polynomial arithmetic.
    void normalize() {
        for (auto& coefficient : coefficients_) {
            if (coefficient.get_den() == 0) {
                throw std::invalid_argument("A rational denominator must not be zero.");
            }
            coefficient.canonicalize();
        }
        detail::trim(coefficients_);
        if (coefficients_.size() == 1 && coefficients_[0] == 0) {
            throw std::invalid_argument("The zero polynomial has infinitely many roots.");
        }
        detail::primitive(coefficients_);
    }

    detail::RationalPoly coefficients_;
};

} // namespace geode
