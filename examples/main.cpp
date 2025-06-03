// main.cpp
#include "geode_root_solver.hpp"
#include <iostream>
#include <vector>
#include <iomanip> // For std::setw (used for alignment)
#include <string>  // For std::string manipulation in print_polynomial

// Helper function to print a polynomial in a more readable format
template <typename R>
void print_polynomial(const geode::Polynomial<R>& p, const std::string& poly_name = "P(x)") {
    std::cout << poly_name << " = ";
    const auto& coeffs = p.coeff();
    if (p.degree() == 0 && coeffs[0] == R(0L)) {
        std::cout << "0" << std::endl;
        return;
    }

    bool first_term = true;
    for (int i = static_cast<int>(p.degree()); i >= 0; --i) {
        R current_coeff = coeffs[static_cast<std::size_t>(i)];
        if (current_coeff == R(0L) && !(i == 0 && first_term)) { // Skip zero terms unless it's the only term P(x)=0
            continue;
        }

        std::string coeff_val_str;
        R abs_coeff = geode::abs(current_coeff);

        // Sign
        if (current_coeff < R(0L)) {
            std::cout << (first_term ? "-" : " - ");
        } else if (!first_term) {
            std::cout << " + ";
        }
        
        first_term = false;

        // Coefficient value (don't print 1 for 1x^n unless it's the constant term 1)
        bool print_coeff_val = true;
        if (abs_coeff == R(1L) && i != 0) {
            print_coeff_val = false;
        }
        
        if (print_coeff_val) {
            std::cout << abs_coeff.to_string();
        }


        // Variable and power
        if (i > 0) {
            std::cout << "x";
            if (i > 1) {
                std::cout << "^" << i;
            }
        }
    }
    if (first_term) { // Happens if polynomial was all zeros (e.g. P(x)=0 after trim)
        std::cout << "0";
    }
    std::cout << std::endl;
}


int main() {
    // 1. Set global precision for MPFR (HighPrecFloat)
    unsigned bits_of_precision = 128; // Example: 128 bits (approx. 38 decimal digits)
    geode::set_global_bits(bits_of_precision);
    std::cout << "Using MPFR global precision: " << geode::get_global_bits() << " bits." << std::endl;

    // Tolerance for Newton's method and for series Fmax calculation
    // A tolerance of 1e-N implies N decimal digits of accuracy desired.
    // For 128 bits ( ~38 decimal digits), a tol of 1e-30 is reasonable.
    geode::HighPrecFloat tol("1e-30"); 
    std::cout << "Using tolerance for root finding: " << tol.to_string() << std::endl;

    // --- Example 1: P(x) = x^3 - x ---
    // Roots: -1, 0, 1
    // Coefficients (a0, a1, a2, a3): [0, -1, 0, 1]
    std::cout << "\n--- Example 1 ---" << std::endl;
    std::vector<geode::HighPrecFloat> coeffs1 = {
        geode::HighPrecFloat("0.0"), 
        geode::HighPrecFloat("-1.0"),
        geode::HighPrecFloat("0.0"),
        geode::HighPrecFloat("1.0")
    };
    geode::Polynomial<geode::HighPrecFloat> poly1(coeffs1);
    print_polynomial(poly1, "P1(x)");
    std::cout << "Degree of P1(x): " << poly1.degree() << std::endl;

    geode::RootSolver<geode::HighPrecFloat> solver1(poly1);
    std::vector<geode::HighPrecFloat> roots1 = solver1.compute_roots(tol);

    std::cout << "Computed roots for P1(x):" << std::endl;
    for (size_t i = 0; i < roots1.size(); ++i) {
        geode::HighPrecFloat p_at_root = poly1(roots1[i]);
        std::cout << "Root " << std::setw(2) << i + 1 << ": " << std::setw(45) << roots1[i].to_string()
                  << "  |P(root)|: " << geode::abs(p_at_root).to_string() << std::endl;
    }

    // --- Example 2: P(x) = x^2 - 2 ---
    // Roots: sqrt(2), -sqrt(2) (approx 1.41421356, -1.41421356)
    // Coefficients (a0, a1, a2): [-2, 0, 1]
    std::cout << "\n--- Example 2 ---" << std::endl;
    std::vector<geode::HighPrecFloat> coeffs2 = {
        geode::HighPrecFloat("-2.0"), 
        geode::HighPrecFloat("0.0"),
        geode::HighPrecFloat("1.0")
    };
    geode::Polynomial<geode::HighPrecFloat> poly2(coeffs2);
    print_polynomial(poly2, "P2(x)");
    std::cout << "Degree of P2(x): " << poly2.degree() << std::endl;

    geode::RootSolver<geode::HighPrecFloat> solver2(poly2);
    std::vector<geode::HighPrecFloat> roots2 = solver2.compute_roots(tol);

    std::cout << "Computed roots for P2(x):" << std::endl;
    for (size_t i = 0; i < roots2.size(); ++i) {
        geode::HighPrecFloat p_at_root = poly2(roots2[i]);
        std::cout << "Root " << std::setw(2) << i + 1 << ": " << std::setw(45) << roots2[i].to_string()
                  << "  |P(root)|: " << geode::abs(p_at_root).to_string() << std::endl;
    }
    geode::HighPrecFloat sqrt2_ref = geode::pow(geode::HighPrecFloat("2.0"), geode::HighPrecFloat("0.5"));
    std::cout << "Reference sqrt(2): " << std::setw(40) << sqrt2_ref.to_string() << std::endl;


    // --- Example 3: P(x) = (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 ---
    // Roots: 1, 2, 3
    // Coefficients (a0, a1, a2, a3): [-6, 11, -6, 1]
    std::cout << "\n--- Example 3 ---" << std::endl;
    std::vector<geode::HighPrecFloat> coeffs3 = {
        geode::HighPrecFloat("-6.0"), 
        geode::HighPrecFloat("11.0"),
        geode::HighPrecFloat("-6.0"),
        geode::HighPrecFloat("1.0")
    };
    geode::Polynomial<geode::HighPrecFloat> poly3(coeffs3);
    print_polynomial(poly3, "P3(x)");

    geode::RootSolver<geode::HighPrecFloat> solver3(poly3);
    std::vector<geode::HighPrecFloat> roots3 = solver3.compute_roots(tol);
    std::cout << "Computed roots for P3(x):" << std::endl;
    for (size_t i = 0; i < roots3.size(); ++i) {
        geode::HighPrecFloat p_at_root = poly3(roots3[i]);
        std::cout << "Root " << std::setw(2) << i + 1 << ": " << std::setw(45) << roots3[i].to_string()
                  << "  |P(root)|: " << geode::abs(p_at_root).to_string() << std::endl;
    }
    
    // --- Example 4: Polynomial of degree 1: P(x) = 2x - 5 ---
    // Root: 2.5
    // Coefficients (a0, a1): [-5, 2]
    std::cout << "\n--- Example 4 ---" << std::endl;
    std::vector<geode::HighPrecFloat> coeffs4 = {
        geode::HighPrecFloat("-5.0"), 
        geode::HighPrecFloat("2.0")
    };
    geode::Polynomial<geode::HighPrecFloat> poly4(coeffs4);
    print_polynomial(poly4, "P4(x)");

    geode::RootSolver<geode::HighPrecFloat> solver4(poly4);
    std::vector<geode::HighPrecFloat> roots4 = solver4.compute_roots(tol);
    std::cout << "Computed roots for P4(x):" << std::endl;
    for (size_t i = 0; i < roots4.size(); ++i) {
        geode::HighPrecFloat p_at_root = poly4(roots4[i]);
        std::cout << "Root " << std::setw(2) << i + 1 << ": " << std::setw(45) << roots4[i].to_string()
                  << "  |P(root)|: " << geode::abs(p_at_root).to_string() << std::endl;
    }


    // --- Example 5: Constant Polynomial P(x) = 5 ---
    // Roots: None
    // Coefficients (a0): [5]
    std::cout << "\n--- Example 5 ---" << std::endl;
    std::vector<geode::HighPrecFloat> coeffs5 = {
        geode::HighPrecFloat("5.0")
    };
    geode::Polynomial<geode::HighPrecFloat> poly5(coeffs5);
    print_polynomial(poly5, "P5(x)");
    try {
        geode::RootSolver<geode::HighPrecFloat> solver5(poly5); // Constructor handles d=0
        std::vector<geode::HighPrecFloat> roots5 = solver5.compute_roots(tol); // Should return empty
        std::cout << "Computed roots for P5(x): (expecting 0 roots)" << std::endl;
        if (roots5.empty()) {
            std::cout << "Correctly found 0 roots." << std::endl;
        } else {
            for (size_t i = 0; i < roots5.size(); ++i) {
                std::cout << "Root " << std::setw(2) << i + 1 << ": " << roots5[i].to_string() << std::endl;
            }
        }
    } catch (const std::invalid_argument& e) {
        std::cout << "Caught expected exception for constant polynomial: " << e.what() << std::endl;
    }

    return 0;
}