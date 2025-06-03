#pragma once

// Disable rarely‑used stuff for MSVC speed‑up
#ifdef _MSC_VER
#  define NOMINMAX
#  include <windows.h>
#endif

#include <unordered_map>
#include <vector>
#include <array>
#include <numeric>
#include <execution>
#include <mutex>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>

#include <gmpxx.h>  // For mpq_class (GMP C++ interface, library is LGPL)
#include <mpfr.h>   // For MPFR C API (LGPL)

namespace geode {

class MpfrFloat;

// MPFR Math functions
MpfrFloat abs(const MpfrFloat& a);
MpfrFloat pow(const MpfrFloat& base, const MpfrFloat& exp);
MpfrFloat pow(const MpfrFloat& base, unsigned long int exp_ui);
MpfrFloat pow(const MpfrFloat& base, long int exp_si); 
MpfrFloat cos(const MpfrFloat& a);
MpfrFloat sin(const MpfrFloat& a); 
MpfrFloat acos(const MpfrFloat& a);
MpfrFloat log(const MpfrFloat& a);
MpfrFloat ceil(const MpfrFloat& a);
MpfrFloat polar(const MpfrFloat& mag, const MpfrFloat& ang); // Custom: mag * cos(ang)

// Custom C++ wrapper for mpfr_t to replace mpfr::mpreal
class MpfrFloat {
public:
    mpfr_t m_val; // mpfr_t is an array type: struct __mpfr_struct[1]

private:
    static mpfr_prec_t s_current_prec; // Stores current global precision

public:
    // Global precision controller
    static void set_global_precision_bits(unsigned bits) {
        if (bits < MPFR_PREC_MIN) {
            // MPFR_PREC_MAX is usually very large, check against practical limits if necessary
            throw std::out_of_range("MPFR precision bits too low.");
        }
        s_current_prec = static_cast<mpfr_prec_t>(bits);
        mpfr_set_default_prec(s_current_prec); // Set MPFR's internal default for new mpfr_t's
    }
    static unsigned get_global_precision_bits() {
        return static_cast<unsigned>(s_current_prec);
    }

    // Constructors
    MpfrFloat() { mpfr_init2(m_val, s_current_prec); /* mpfr_init2 also sets to NaN by default; explicitly set to 0 if desired */ mpfr_set_ui(m_val, 0, MPFR_RNDN); }
    explicit MpfrFloat(double d) { mpfr_init2(m_val, s_current_prec); mpfr_set_d(m_val, d, MPFR_RNDN); }
    explicit MpfrFloat(long l) { mpfr_init2(m_val, s_current_prec); mpfr_set_si(m_val, l, MPFR_RNDN); }
    explicit MpfrFloat(unsigned long ul) { mpfr_init2(m_val, s_current_prec); mpfr_set_ui(m_val, ul, MPFR_RNDN); }
    explicit MpfrFloat(const char* s, int base = 10) { mpfr_init2(m_val, s_current_prec); mpfr_set_str(m_val, s, base, MPFR_RNDN); }
    explicit MpfrFloat(const std::string& s, int base = 10) : MpfrFloat(s.c_str(), base) {}
    explicit MpfrFloat(const mpq_class& q) {
        mpfr_init2(m_val, s_current_prec);
        mpfr_set_q(m_val, q.get_mpq_t(), MPFR_RNDN);
    }

    // Copy semantics
    MpfrFloat(const MpfrFloat& other) {
        mpfr_init2(m_val, mpfr_get_prec(other.m_val)); // Initialize with other's precision
        mpfr_set(m_val, other.m_val, MPFR_RNDN);
    }
    MpfrFloat& operator=(const MpfrFloat& other) {
        if (this == &other) return *this;
        // mpfr_set stores result using m_val's precision.
        // If you want to match other's precision: mpfr_set_prec(m_val, mpfr_get_prec(other.m_val));
        mpfr_set(m_val, other.m_val, MPFR_RNDN);
        return *this;
    }

    // Move semantics
    MpfrFloat(MpfrFloat&& other) noexcept {
        mpfr_init2(m_val, mpfr_get_prec(other.m_val)); // Init m_val to a valid state
        mpfr_swap(m_val, other.m_val); // Swap contents
        // other.m_val now holds the freshly initialized state, will be cleared by its destructor.
    }
    MpfrFloat& operator=(MpfrFloat&& other) noexcept {
        if (this == &other) return *this;
        mpfr_swap(m_val, other.m_val); // Both are valid and initialized.
        return *this;
    }

    // Destructor
    ~MpfrFloat() { mpfr_clear(m_val); }

    // Access to raw mpfr_t
    mpfr_ptr get_mpfr_ptr() { return m_val; } 
    mpfr_srcptr get_mpfr_srcptr() const { return m_val; }

    // To string for debugging/output
    std::string to_string(mpfr_rnd_t rnd = MPFR_RNDN) const {
        if (mpfr_nan_p(m_val)) return "nan";
        if (mpfr_inf_p(m_val)) return (mpfr_signbit(m_val) ? "-inf" : "inf");
        
        // Determine a suitable number of digits for n_digits in mpfr_get_str
        // Roughly: bits * log10(2) + a small margin
        size_t n_digits = static_cast<size_t>(mpfr_get_prec(m_val) * std::log10(2.0)) + 2;

        mpfr_exp_t exp;
        char* s = mpfr_get_str(nullptr, &exp, 10, n_digits, m_val, rnd);
        if (s == nullptr) return "[error converting to string]";
        
        std::string res;
        if (s[0] == '-') { res += '-'; res += "0."; res += (s + 1); }
        else if (s[0] == '+') { res += "0."; res += (s + 1); } // Some modes might prepend '+'
        else { res += "0."; res += s; }
        res += 'E';
        res += std::to_string(exp);
        
        mpfr_free_str(s);
        return res;
    }
    
    // Basic Arithmetic
    MpfrFloat operator-() const { MpfrFloat r; mpfr_neg(r.m_val, m_val, MPFR_RNDN); return r; }

    MpfrFloat& operator+=(const MpfrFloat& rhs) { mpfr_add(m_val, m_val, rhs.m_val, MPFR_RNDN); return *this; }
    MpfrFloat& operator-=(const MpfrFloat& rhs) { mpfr_sub(m_val, m_val, rhs.m_val, MPFR_RNDN); return *this; }
    MpfrFloat& operator*=(const MpfrFloat& rhs) { mpfr_mul(m_val, m_val, rhs.m_val, MPFR_RNDN); return *this; }
    MpfrFloat& operator/=(const MpfrFloat& rhs) { mpfr_div(m_val, m_val, rhs.m_val, MPFR_RNDN); return *this; }

    friend MpfrFloat operator+(MpfrFloat lhs, const MpfrFloat& rhs) { lhs += rhs; return lhs; }
    friend MpfrFloat operator-(MpfrFloat lhs, const MpfrFloat& rhs) { lhs -= rhs; return lhs; }
    friend MpfrFloat operator*(MpfrFloat lhs, const MpfrFloat& rhs) { lhs *= rhs; return lhs; }
    friend MpfrFloat operator/(MpfrFloat lhs, const MpfrFloat& rhs) { lhs /= rhs; return lhs; }

    // Comparison operators
    friend bool operator<(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_cmp(a.m_val, b.m_val) < 0; }
    friend bool operator>(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_cmp(a.m_val, b.m_val) > 0; }
    friend bool operator<=(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_cmp(a.m_val, b.m_val) <= 0; }
    friend bool operator>=(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_cmp(a.m_val, b.m_val) >= 0; }
    friend bool operator==(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_equal_p(a.m_val, b.m_val) != 0; }
    friend bool operator!=(const MpfrFloat& a, const MpfrFloat& b) { return mpfr_equal_p(a.m_val, b.m_val) == 0; }

    // Allow R{0} by providing constructor from int, implicitly used for 0L
    MpfrFloat(int i) : MpfrFloat(static_cast<long>(i)) {}
};

// Initialize static member - default precision (e.g., 256 bits).
// Call set_global_precision_bits() early in main() to change.
mpfr_prec_t MpfrFloat::s_current_prec = 256; 

// Definitions for MPFR Math functions using MpfrFloat
inline MpfrFloat abs(const MpfrFloat& a) { MpfrFloat r; mpfr_abs(r.m_val, a.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat pow(const MpfrFloat& base, const MpfrFloat& exp) { MpfrFloat r; mpfr_pow(r.m_val, base.get_mpfr_srcptr(), exp.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat pow(const MpfrFloat& base, unsigned long int exp_ui) { MpfrFloat r; mpfr_pow_ui(r.m_val, base.get_mpfr_srcptr(), exp_ui, MPFR_RNDN); return r; }
inline MpfrFloat pow(const MpfrFloat& base, long int exp_si) { MpfrFloat r; mpfr_pow_si(r.m_val, base.get_mpfr_srcptr(), exp_si, MPFR_RNDN); return r; }
inline MpfrFloat cos(const MpfrFloat& a) { MpfrFloat r; mpfr_cos(r.m_val, a.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat sin(const MpfrFloat& a) { MpfrFloat r; mpfr_sin(r.m_val, a.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat acos(const MpfrFloat& a) { MpfrFloat r; mpfr_acos(r.m_val, a.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat log(const MpfrFloat& a) { MpfrFloat r; mpfr_log(r.m_val, a.get_mpfr_srcptr(), MPFR_RNDN); return r; }
inline MpfrFloat ceil(const MpfrFloat& a) { MpfrFloat r; mpfr_ceil(r.m_val, a.get_mpfr_srcptr()); return r; }
inline MpfrFloat polar(const MpfrFloat& mag, const MpfrFloat& ang) { return mag * cos(ang); }


// -----------------------------------------------------------------------------
// EXACT & INTERVAL ARITHMETIC LAYER
// -----------------------------------------------------------------------------

using mpq   = mpq_class;
using HighPrecFloat = MpfrFloat;

struct Ball {
    HighPrecFloat mid;
    HighPrecFloat rad;
    Ball() = default;
    Ball(const HighPrecFloat& m, const HighPrecFloat& r) : mid(m), rad(r) {
        if (r < HighPrecFloat(0L)) throw std::invalid_argument("negative radius");
    }
    bool contains(const HighPrecFloat& x) const { return abs(x - mid) <= rad; }
};

inline void set_global_bits(unsigned bits) { MpfrFloat::set_global_precision_bits(bits); }
inline unsigned get_global_bits()         { return MpfrFloat::get_global_precision_bits(); }

// High-precision PI
const HighPrecFloat PI = acos(HighPrecFloat(-1.0));


// -----------------------------------------------------------------------------
// POLYNOMIAL (dense) – templated on coefficient type R.
// -----------------------------------------------------------------------------

template <typename R>
class Polynomial {
public:
    using value_type = R; // Used by derivative()
    explicit Polynomial(std::vector<R> c) : coeff_(std::move(c)) { trim(); }

    std::size_t degree() const noexcept { return coeff_.empty() ? 0 : coeff_.size() - 1; }
    const std::vector<R>& coeff() const noexcept { return coeff_; }

    template <typename X> // X is the type for evaluation
    X operator()(X x_val) const {
        if (coeff_.empty()) return X(0L); 
        
        X y_val(0L); 
        if (!coeff_.empty()) { // Should always be true due to trim()
            auto it = coeff_.rbegin();
            y_val = static_cast<X>(*it); // Initialize with the highest degree coefficient
            ++it;
            for (; it != coeff_.rend(); ++it) {
                y_val = y_val * x_val + static_cast<X>(*it);
            }
        }
        return y_val;
    }

    Polynomial derivative() const {
        if (coeff_.size() <= 1) return Polynomial({R(0L)});
        std::vector<R> d_coeffs(coeff_.size() - 1);
        for (std::size_t i = 1; i < coeff_.size(); ++i) {
            d_coeffs[i - 1] = coeff_[i] * R(static_cast<long>(i)); // Fixed
        }
        return Polynomial(std::move(d_coeffs));
    }

private:
    void trim() {
        while (!coeff_.empty() && coeff_.back() == R(0L)) coeff_.pop_back();
        if (coeff_.empty()) coeff_.push_back(R(0L));
    }
    std::vector<R> coeff_;
};

// -----------------------------------------------------------------------------
// GEODE TABLE  – memoised coefficients
// -----------------------------------------------------------------------------
class GeodeTable {
public:
    GeodeTable(unsigned degree_of_P, unsigned Fmax_val) // degree_of_P is d in P(x)=x^d+...
        : deg_param_size_(degree_of_P), Fmax_(Fmax_val) // deg_param_size_ is size of m-vector for t_k
    {
        if (deg_param_size_ == 0 && Fmax_val > 0) {
            // This case might need review if deg_param_size_ can legitimately be 0
            // for some polynomial forms, but typically it's P.degree() >= 1
        }
        // Zero multi‑index → 1
        // For deg_param_size_=0, Index is empty vector. table_[{}] = 1
        table_[Index(deg_param_size_, 0)] = mpq{1};
        // Build breadth‑first by |m|
        for (unsigned F_val = 1; F_val <= Fmax_; ++F_val) {
            build_shell(F_val);
        }
    }

    mpq get(const std::vector<unsigned>& m) const {
        auto it = table_.find(m);
        if (it == table_.end()) return mpq{0}; // Return exact 0 if not found
        return it->second;
    }

private:
    using Index = std::vector<unsigned>;
    struct IndexHash {
        std::size_t operator()(Index const& v) const noexcept {
            std::size_t h = 0;
            for (auto x : v) h = (h * 1315423911u) ^ std::hash<unsigned>{}(x);
            return h;
        }
    };

    void build_shell(unsigned F_val) { // Renamed F to F_val
        std::vector<Index> shell_indices;
        Index current_m(deg_param_size_, 0);
        if (deg_param_size_ > 0) {
             enumerate_indices(shell_indices, current_m, 0, F_val);
        } else if (F_val == 0) { // Only relevant if loop started at F=0
            // For deg_param_size_ = 0, only |m|=0 makes sense.
            // This is handled by constructor. build_shell is called for F_val >= 1.
            // If somehow F_val is 0 here and deg_param_size_ is 0:
            // shell_indices.push_back({}); // G(empty_m with |empty_m|=0) = 1 already set
        } // else if deg_param_size_ == 0 and F_val > 0, no such m, shell_indices remains empty.

        if (shell_indices.empty() && F_val > 0) { // No indices to compute for this shell
            return;
        }

        // Temporary storage for results from this shell's parallel computation
        std::vector<std::pair<Index, mpq>> results_for_this_shell(shell_indices.size());

        // Parallel computation phase:
        // Reads from 'table_' for G values where |m| < F_val (which are stable).
        // Writes results to 'results_for_this_shell', not 'table_'.
        // Using std::transform as it's designed for this output pattern.
        std::transform(std::execution::par, shell_indices.begin(), shell_indices.end(), results_for_this_shell.begin(),
            [&](const Index& idx_m_const) -> std::pair<Index, mpq> { // Use a different name for clarity
            mpq sum_val = 0;
            Index m_minus_ek = idx_m_const; // Create a mutable copy for m - e_k
            for (unsigned k_comp = 0; k_comp < deg_param_size_; ++k_comp) {
                if (idx_m_const[k_comp] > 0) { // Check original idx_m_const[k_comp]
                    m_minus_ek[k_comp] = idx_m_const[k_comp] - 1; // Set specific component
                    // All other components of m_minus_ek are same as idx_m_const
                    // This ensures m_minus_ek is correctly m - e_k from original idx_m_const

                    mpq prev_g = get(m_minus_ek); // 'get' reads stable F-1 shell

                    m_minus_ek[k_comp] = idx_m_const[k_comp]; // Restore for next k_comp or just use idx_m_const to build next m_minus_ek

                    sum_val += mpq(k_comp + 1) * mpq(idx_m_const[k_comp]) * prev_g;
                }
            }
            // F_val will be >= 1 here because the outer loop in constructor starts from F_val=1
            mpq g_result = sum_val / mpq(F_val);
            return {idx_m_const, g_result}; // Return the original index and its computed G value
        });

        // Serial insertion phase:
        // Lock the table and insert all computed results for this shell.
        std::lock_guard<std::mutex> lock(table_mtx_);
        for (const auto& p : results_for_this_shell) {
            table_[p.first] = p.second;
        }
    }

    // Renamed from 'enumerate' to 'enumerate_indices'
    void enumerate_indices(std::vector<Index>& out_indices, Index& current_m_build, unsigned current_pos, unsigned remaining_sum) {
        if (current_pos == deg_param_size_ - 1) { // last component takes the rest
            current_m_build[current_pos] = remaining_sum;
            out_indices.push_back(current_m_build);
            return;
        }
        for (unsigned i = 0; i <= remaining_sum; ++i) {
            current_m_build[current_pos] = i;
            enumerate_indices(out_indices, current_m_build, current_pos + 1, remaining_sum - i);
        }
    }

    unsigned deg_param_size_; // Size of multi-indices m (number of t_k parameters, typically P.degree())
    unsigned Fmax_;
    std::unordered_map<Index, mpq, IndexHash> table_;
    std::mutex table_mtx_;
};


// -----------------------------------------------------------------------------
// ROOT SOLVER 
// -----------------------------------------------------------------------------

template <typename Float = HighPrecFloat>
class RootSolver {
public:
    using PolyF = Polynomial<Float>;
    explicit RootSolver(const PolyF& P_in) : P_(P_in), d_poly_degree_(P_in.degree()) {
        if (d_poly_degree_ == 0 && P_.coeff()[0] == Float(0L)) { /* error */ }
        // If d_poly_degree_ == 0 (constant P(x)=c0), handle separately (no roots unless c0=0)

        // Populate c_coeffs_WR from P_.coeff()
        c_coeffs_WR.resize(d_poly_degree_ + 1);
        if (!P_.coeff().empty()) {
            c_coeffs_WR[0] = P_.coeff()[0]; // c_0 = p_0
        } else { // Should not happen if P_ is well-formed
            c_coeffs_WR[0] = Float(0L);
        }

        if (d_poly_degree_ >= 1) {
            c_coeffs_WR[1] = -P_.coeff()[1]; // c_1 = -p_1
        }
        for (std::size_t k = 2; k <= d_poly_degree_; ++k) {
            c_coeffs_WR[k] = P_.coeff()[k]; // c_k = p_k for k >= 2
        }

        // Check for c_1 == 0 (important for W&R Theorem 4)
        if (d_poly_degree_ >= 1 && c_coeffs_WR[1] == Float(0L)) {
            // W&R Theorem 4 explicitly requires c1 != 0.
            // "The formula ... fails when c1 = 0..." (page 6 of W&R PDF)
            // This case needs a different approach (e.g., change of variables, or if c0=0, x=0 is a root)
            // For now, indicate this solver cannot proceed.
            std::cerr << "Warning: Coefficient c1 (for W&R formula) is zero. Solver may fail or give incorrect results." << std::endl;
            // Perhaps let it proceed if c0 is also 0, then x=0 is a root, factor it out.
        }
        build_WR_t_values();
    }

    std::vector<Float> compute_roots(Float tol = Float("1e-30")) {
        if (d_poly_degree_ == 0) { // Constant P(x) = p0
            if (c_coeffs_WR[0] == Float(0L)) { /* P(x)=0, infinite roots - ill-defined */ }
            return {}; // No roots for non-zero constant
        }

        // Handle linear case c0 - c1 x = 0 => x = c0/c1
        if (d_poly_degree_ == 1) {
            if (c_coeffs_WR[1] == Float(0L)) { // 0*x = c0
                return (c_coeffs_WR[0] == Float(0L)) ? std::vector<Float>{Float(0L)} : std::vector<Float>{}; // Or throw for 0=c0!=0
            }
            return {c_coeffs_WR[0] / c_coeffs_WR[1]};
        }
        
        // W&R requires c1 != 0 for their Theorem 4 series.
        if (c_coeffs_WR[1] == Float(0L)) {
            std::cerr << "Error: W&R method requires non-zero c1 coefficient. Cannot proceed." << std::endl;
            // Could try other methods or transformations if c1 is zero.
            // e.g. if c0=0, then x=0 is a root. Factor out x and solve for remaining poly.
            if(c_coeffs_WR[0] == Float(0L)) {
                // TODO: Implement factorization and recurse. For now, just note x=0.
                std::cout << "Found root x=0 due to c0=0, c1=0. Other roots not found by this method." << std::endl;
                // A more robust solver would factor out x and solve P(x)/x = 0.
            }
            return {}; 
        }


        // Fmax calculation
        unsigned Fmax_val;
        if (tol == Float(0L) || tol < Float(0L)) { 
            Fmax_val = (get_global_bits() / 2) + 4; // Reduced safety margin a bit for faster tests
        } else {
            Float log_tol = geode::log(tol); // Ensure geode::log is defined for MpfrFloat
            Float log_2 = geode::log(Float(2.0)); 
            Float term = geode::ceil(-(log_tol / log_2)); // Ensure geode::ceil is defined
            
            double fmax_double_val;
            mpfr_t temp_mpfr_val; 
            mpfr_init2(temp_mpfr_val, mpfr_get_prec(term.get_mpfr_srcptr()));
            mpfr_set(temp_mpfr_val, term.get_mpfr_srcptr(), MPFR_RNDN);
            fmax_double_val = mpfr_get_d(temp_mpfr_val, MPFR_RNDN);
            mpfr_clear(temp_mpfr_val);
            
            Fmax_val = (fmax_double_val < 0) ? 0 : static_cast<unsigned>(fmax_double_val);
            Fmax_val += 2; // Reduced safety margin
        }
        // Ensure Fmax is not excessively large initially, e.g. cap at 15-20 for testing
        // Fmax_val = std::min(Fmax_val, 10u); // For faster initial tests

        // GeodeTable constructor takes number of t_k variables for size of m-vectors.
        // We have t_2, ..., t_D, so D-1 variables if D >= 2.
        unsigned num_t_vars = (d_poly_degree_ >= 2) ? (d_poly_degree_ - 1) : 0;
        GeodeTable C_table(num_t_vars, Fmax_val); // Computes C_m = G.get(m)

        Float alpha = series_eval_alpha(C_table, Fmax_val);

        Float c0 = c_coeffs_WR[0];
        Float c1 = c_coeffs_WR[1];
        Float x_initial_guess = (c0 / c1) * alpha;

        // For now, we only compute one root using this direct application.
        // Finding other roots would require bootstrapping (as W&R do for cubic)
        // or a more complex way to select branches of alpha.
        std::vector<Float> roots;
        Float root = newton_polish(x_initial_guess, tol);
        roots.push_back(root);

        // W&R use bootstrapping for other roots:
        // If x_approx is one root, consider P(x)/(x-x_approx) or P(x_approx + delta_x)
        // This part is more heuristic for finding *all* roots.
        // The "phase_factor" was an attempt for this, but it's not directly in W&R's Theorem 4.

        // If you want to try the phase factor approach (more speculative with W&R's direct formula):
        // It would apply to x_initial_guess.
        // For example, if x_initial_guess is a principal value, other guesses could be
        // x_initial_guess * omega_k where omega_k are D-th roots of unity.
        // This makes sense if the original polynomial was x^D - C = 0.
        // For general polynomials, it's less clear if this simple phase works directly on x.

        // Let's try to find other roots by using the D-th roots of unity on alpha
        // (as if alpha itself had D branches, and c0/c1 was a common factor)
        // This is an INTERPRETATION, not directly from W&R Thm 4 for multiple roots.

        if (d_poly_degree_ > 1) { // Only try for multiple roots if deg > 1
            for (unsigned j = 1; j < d_poly_degree_; ++j) { // Start from j=1 for other phases
                Float num = Float(2L) * PI * Float(static_cast<long>(j));
                Float den = Float(static_cast<long>(d_poly_degree_));
                Float theta = num / den;
                // This is a guess at how to get other roots: assume alpha has branches
                // or that the overall solution x has branches differing by roots of unity.
                // More robustly, one should use deflation or different x0 for transformations.
                Float phase_factor_real = geode::cos(theta); // std::polar gives complex
                // We are in real arithmetic, so this phase factor only gives 1 or -1 if theta is multiple of PI
                // If we need to explore more, the problem becomes complex.
                // Let's assume for now that different starting Newton points might be generated this way.
                // For real roots, this might not be very effective unless D is even and gives -1.

                Float x_guess_variant = x_initial_guess * phase_factor_real; // Simple scaling
                // Or, more speculatively: (c0/c1) * (alpha * phase_factor_real)
                // If alpha is the thing that branches.
                // This part is highly heuristic with real arithmetic.

                // A common technique if Newton gets stuck or to find other roots
                // is to perturb the initial guess or use a different one.
                // For demo, let's try a very simple perturbation if phase_factor_real is 1
                // This part is experimental.
                if (j > 0) {
                    Float perturbation_factor = (j % 2 == 1) ? Float("0.8") : Float("1.2");
                    if (abs(x_initial_guess) > tol) { // Avoid perturbing zero too much
                        x_guess_variant = x_initial_guess * perturbation_factor * phase_factor_real;
                    } else { // If initial guess is near zero, perturb from a small non-zero
                        x_guess_variant = ( (j%2==0) ? Float("0.1") : Float("-0.1") ) * Float(static_cast<long>(j));
                    }
                }
                
                bool already_found = false;
                for(const auto& r : roots) {
                    if (abs(r - x_guess_variant) < tol * Float("100")) { // Check if too close to an existing root
                        already_found = true;
                        break;
                    }
                }
                if (!already_found) {
                    Float another_root = newton_polish(x_guess_variant, tol);
                    // Add only if significantly different from already found roots
                    bool distinct_enough = true;
                    for(const auto& r_found : roots) {
                        if (abs(another_root - r_found) < tol * Float("1000")) {
                            distinct_enough = false;
                            break;
                        }
                    }
                    if(distinct_enough) roots.push_back(another_root);
                }
                if (roots.size() >= d_poly_degree_) break; // Found enough roots
            }
        }

        return roots;
    }

private:

    void build_WR_t_values() {
        // t_k = c_0^{k-1} c_k / c_1^k for k=2, ..., D (d_poly_degree_)
        // t_values_WR will store t_2, t_3, ...
        if (d_poly_degree_ < 2) { // No t_k needed if degree < 2 (linear or constant)
            return;
        }
        t_values_WR.resize(d_poly_degree_ - 1); // t_2 to t_D (D-1 elements)

        Float c0 = c_coeffs_WR[0];
        Float c1 = c_coeffs_WR[1];

        if (c1 == Float(0L)) {
            // Handle division by zero if c1 is zero.
            // As noted, W&R formula assumes c1 != 0.
            // Fill t_values_WR with NaN or throw, or set a flag.
            for(size_t i=0; i<t_values_WR.size(); ++i) t_values_WR[i] = Float(0L); // Placeholder
            return;
        }

        for (std::size_t k = 2; k <= d_poly_degree_; ++k) {
            // t_k maps to t_values_WR[k-2]
            Float ck = c_coeffs_WR[k];
            Float c0_pow_km1 = pow(c0, static_cast<unsigned long>(k - 1));
            Float c1_pow_k = pow(c1, static_cast<unsigned long>(k));
            t_values_WR[k - 2] = (c0_pow_km1 * ck) / c1_pow_k;
        }
    }

    // series_eval_alpha uses the existing GeodeTable (which computes C_m)
    // and the t_values_WR.
    // It computes S[t_2, ..., t_D] = Sum C_m * product(t_k^m_k)
    Float series_eval_alpha(const GeodeTable& C_table, unsigned Fmax_val) {
        Float alpha_sum = Float(0L); // Start with C_0=0 and add it if |m|=0 is handled by loop
                                     // W&R's S starts with C_0 = 1.
                                     // The GeodeTable for m={} (size 0 vector) gives 1.
        
        // Vector m has D-1 components (m_2, ..., m_D)
        // The size of t_values_WR is D-1.
        // deg_param_size for GeodeTable should be D-1 (number of t_k variables).
        std::vector<unsigned> m_indices(t_values_WR.size(), 0);

        // Sum over |m| = F, from F=0 up to Fmax_val
        for (unsigned F_sum = 0; F_sum <= Fmax_val; ++F_sum) {
            if (t_values_WR.empty() && F_sum > 0) continue; // No t_k for F > 0
            if (!t_values_WR.empty() || F_sum == 0) { // Process F=0 or if t_k exist
                Float sum_for_shell_F = shell_sum_for_F_alpha_recursive(C_table, m_indices, 0, F_sum);
                alpha_sum += sum_for_shell_F;
            }
        }
        return alpha_sum;
    }

    // Recursive helper for Sum C_m * product(t_k^m_k) for a given |m|=F_val
    Float shell_sum_for_F_alpha_recursive(const GeodeTable& C_table, std::vector<unsigned>& m_vec,
                                         unsigned k_idx_component, unsigned remaining_F_val) {
        // k_idx_component is index into m_vec (0 to t_values_WR.size()-1)
        // m_vec[k_idx_component] corresponds to m_{k_idx_component + 2}
        // t_values_WR[k_idx_component] corresponds to t_{k_idx_component + 2}

        if (t_values_WR.empty()) { // Special case for D=0 or D=1 (no t_k)
            return (remaining_F_val == 0) ? mpq_to_float(C_table.get(m_vec)) : Float(0L); // C_table.get({}) is 1
        }

        if (k_idx_component == t_values_WR.size() - 1) { // Last component of m_vec
            m_vec[k_idx_component] = remaining_F_val;

            mpq Cm = C_table.get(m_vec); // Cm = G.get(m)
            if (Cm == 0) return Float(0L);

            Float term_prod_tk_mk = Float(1L);
            for (unsigned i = 0; i < t_values_WR.size(); ++i) {
                if (m_vec[i] > 0) {
                    term_prod_tk_mk *= pow(t_values_WR[i], static_cast<unsigned long>(m_vec[i]));
                }
            }
            return mpq_to_float(Cm) * term_prod_tk_mk;
        }

        Float total_sum_for_level = Float(0L);
        for (unsigned val_for_mk = 0; val_for_mk <= remaining_F_val; ++val_for_mk) {
            m_vec[k_idx_component] = val_for_mk;
            total_sum_for_level += shell_sum_for_F_alpha_recursive(C_table, m_vec, k_idx_component + 1, remaining_F_val - val_for_mk);
        }
        return total_sum_for_level;
    }

    Float newton_polish(Float x_initial, Float tol) {
        const unsigned MAX_IT = 60; // Increased iterations slightly
        Float x_current = x_initial;
        PolyF P_prime = P_.derivative(); 

        for (unsigned i = 0; i < MAX_IT; ++i) {
            Float fx_val = P_(x_current);
            if (abs(fx_val) < tol) return x_current;
            
            Float dfx_val = P_prime(x_current);
            if (dfx_val == Float(0L)) { // Derivative is zero
                // std::cerr << "Warning: Newton step failed due to zero derivative near x = " << x_current.to_string() << std::endl;
                return x_current; // Return current best guess
            }
            x_current -= fx_val / dfx_val;
        }
        // std::cerr << "Warning: Newton's method did not converge to tolerance " << tol.to_string() 
        //           << " after " << MAX_IT << " iterations. Last |f(x)| = " << abs(P_(x_current)).to_string() << std::endl;
        return x_current; 
    }

    static Float mpq_to_float(const mpq& r) {
        return Float(r); // MpfrFloat constructor from mpq_class
    }

    PolyF              P_;                 // The polynomial
    std::size_t        d_poly_degree_;     // Degree of P
    std::vector<Float> c_coeffs_WR; // W&R's c_0, c_1, ..., c_D (c_1 will be negated)
    std::vector<Float> t_values_WR; // W&R's t_2, ..., t_D
};

} // namespace geode
