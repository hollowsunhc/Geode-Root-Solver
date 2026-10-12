/**
 * @file solver.cpp
 * @brief Bounded C interface for the browser's isolated WebAssembly solver.
 */
#include "geode_root_solver.hpp"

#include <cstddef>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t max_input = 2048;
constexpr std::size_t max_coefficients = 17;
constexpr std::size_t max_token = 64;

/// @brief Recognize ASCII whitespace without locale-dependent parsing.
bool space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

/// @brief Recognize a base-ten digit.
bool digit(char c) {
    return c >= '0' && c <= '9';
}

/// @brief Validate a bounded decimal before handing it to MPFR.
void check_decimal(const std::string& token) {
    if (token.empty() || token.size() > max_token)
        throw std::invalid_argument("Use at most 64 characters per coefficient.");
    std::size_t i = 0;
    if (token[i] == '+' || token[i] == '-')
        ++i;
    unsigned digits = 0;
    while (i < token.size() && digit(token[i])) {
        ++digits;
        ++i;
    }
    if (i < token.size() && token[i] == '.') {
        ++i;
        while (i < token.size() && digit(token[i])) {
            ++digits;
            ++i;
        }
    }
    if (digits == 0)
        throw std::invalid_argument("Coefficients must be decimal numbers, not expressions.");
    if (i < token.size() && (token[i] == 'e' || token[i] == 'E')) {
        ++i;
        if (i < token.size() && (token[i] == '+' || token[i] == '-'))
            ++i;
        const std::size_t start = i;
        unsigned exponent = 0;
        while (i < token.size() && digit(token[i])) {
            exponent = exponent * 10 + static_cast<unsigned>(token[i] - '0');
            if (exponent > 100)
                throw std::invalid_argument("Coefficient exponents must be between -100 and 100.");
            ++i;
        }
        if (i == start)
            throw std::invalid_argument("An exponent needs decimal digits.");
    }
    if (i != token.size())
        throw std::invalid_argument("Coefficients must be decimal numbers, not expressions.");
}

/**
 * @brief Parse a bounded coefficient list without invoking numeric arithmetic.
 * @param text Null-terminated ASCII decimal coefficients, ascending by power.
 * @return Validated decimal strings.
 * @throws std::invalid_argument For malformed or excessive input.
 */
std::vector<std::string> parse(const char* text) {
    if (text == nullptr)
        throw std::invalid_argument("Enter polynomial coefficients.");
    std::size_t length = 0;
    while (length <= max_input && text[length] != '\0')
        ++length;
    if (length > max_input)
        throw std::invalid_argument("Use at most 2048 characters of coefficient input.");
    std::vector<std::string> tokens;
    std::size_t i = 0;
    bool field_has_token = false;
    while (i < length) {
        if (space(text[i])) {
            ++i;
            continue;
        }
        if (text[i] == ',') {
            if (!field_has_token)
                throw std::invalid_argument("A comma must have a coefficient on each side.");
            field_has_token = false;
            ++i;
            continue;
        }
        const std::size_t start = i;
        while (i < length && !space(text[i]) && text[i] != ',')
            ++i;
        if (tokens.size() == max_coefficients)
            throw std::invalid_argument("The browser lab supports degree 16 or lower.");
        const std::string token(text + start, i - start);
        check_decimal(token);
        tokens.push_back(token);
        field_has_token = true;
    }
    if (tokens.empty() || !field_has_token)
        throw std::invalid_argument("Enter coefficients without empty comma fields.");
    return tokens;
}

/// @brief Quote an exception message for a JSON string without reflecting raw input.
std::string quote(const std::string& value) {
    std::string result = "\"";
    for (unsigned char c : value) {
        if (c == '"' || c == '\\')
            result += '\\';
        result += c < 32 ? ' ' : static_cast<char>(c);
    }
    return result + '"';
}

/**
 * @brief Serialize the binary root exactly, retaining its absolute error bound.
 * @param value Finite MPFR root returned by the solver.
 * @return Its exact decimal expansion in scientific notation, or zero.
 * @note A shortest round-trip string can lose large absolute amounts for large
 *       roots. A dyadic rational has a finite exact decimal expansion.
 */
std::string exact_decimal(const geode::HighPrecFloat& value) {
    if (mpfr_zero_p(value.get_mpfr_srcptr()))
        return "0";
    mpz_class significand;
    const mpfr_exp_t exponent = mpfr_get_z_2exp(significand.get_mpz_t(), value.get_mpfr_srcptr());
    const bool negative = significand < 0;
    if (negative)
        significand = -significand;
    unsigned long places = 0;
    if (exponent >= 0) {
        mpz_mul_2exp(significand.get_mpz_t(), significand.get_mpz_t(),
                     static_cast<mp_bitcnt_t>(exponent));
    } else {
        places = static_cast<unsigned long>(-exponent);
        mpz_class power;
        mpz_ui_pow_ui(power.get_mpz_t(), 5, places);
        significand *= power;
    }
    std::string decimal = significand.get_str();
    const long decimal_exponent = static_cast<long>(decimal.size()) - static_cast<long>(places) - 1;
    while (decimal.size() > 1 && decimal.back() == '0')
        decimal.pop_back();
    if (decimal.size() > 1)
        decimal.insert(1, ".");
    return (negative ? "-" : "") + decimal + "e" + std::to_string(decimal_exponent);
}
} // namespace

/**
 * @brief Solve a browser request after enforcing independent workload limits.
 * @param coefficients Ascending decimal coefficients, separated by commas or whitespace.
 * @param bits Working precision: 128, 256, or 512 bits.
 * @param digits Absolute error exponent: 1..30, 1..60, or 1..120 respectively.
 * @return JSON with exact decimal encodings of the computed binary approximations
 *         and request metadata, or an error string.
 *         The pointer remains valid until the next call in this Wasm instance.
 * @note Run in a disposable worker with a host-side deadline. The memory ceiling
 *       is enforced by the link settings; this synchronous function cannot cancel itself.
 */
extern "C" const char* geode_solve(const char* coefficients, int bits, int digits) {
    static std::string response;
    try {
        const int max_digits = bits == 128 ? 30 : bits == 256 ? 60 : bits == 512 ? 120 : 0;
        if (max_digits == 0 || digits < 1 || digits > max_digits)
            throw std::invalid_argument("Choose a supported precision and error bound.");
        const auto tokens = parse(coefficients);
        geode::set_global_bits(static_cast<unsigned>(bits));
        std::vector<geode::HighPrecFloat> values;
        values.reserve(tokens.size());
        for (const auto& token : tokens)
            values.emplace_back(token);
        const geode::Polynomial<geode::HighPrecFloat> polynomial(std::move(values));
        const geode::RootSolver<> solver(polynomial);
        geode::HighPrecFloat tolerance;
        // Round down so the actual tolerance never exceeds the requested decimal bound.
        mpfr_set_str(tolerance.get_mpfr_ptr(), ("1e-" + std::to_string(digits)).c_str(), 10,
                     MPFR_RNDD);
        const auto roots = solver.compute_roots(tolerance);
        response = "{\"roots\":[";
        for (std::size_t i = 0; i < roots.size(); ++i) {
            if (i != 0)
                response += ',';
            response += quote(exact_decimal(roots[i]));
        }
        response += "],\"degree\":" + std::to_string(polynomial.degree()) +
                    ",\"bits\":" + std::to_string(bits) + ",\"digits\":" + std::to_string(digits) +
                    '}';
    } catch (const std::bad_alloc&) {
        response = "{\"error\":\"Memory limit reached. Try a smaller polynomial.\"}";
    } catch (const std::exception& error) {
        response = "{\"error\":" + quote(error.what()) + '}';
    } catch (...) {
        response = "{\"error\":\"The solver could not finish this calculation.\"}";
    }
    return response.c_str();
}
