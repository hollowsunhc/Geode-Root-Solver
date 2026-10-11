#include "geode_root_solver.hpp"

#include <exception>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    try {
        using Real = geode::HighPrecFloat;
        using Poly = geode::Polynomial<Real>;
        geode::set_global_bits(256);
        const Real tolerance("1e-30");
        const std::vector<std::pair<const char*, Poly>> examples = {
            {"x^3 - x", Poly({0, -1, 0, 1})},
            {"x^2 - 2", Poly({-2, 0, 1})},
            {"(x - 1)(x - 2)(x - 3)", Poly({-6, 11, -6, 1})},
            {"2x - 5", Poly({-5, 2})},
            {"(x - 1)^2", Poly({1, -2, 1})},
            {"x^2 + 1", Poly({1, 0, 1})},
            {"5", Poly({5})}};

        for (const auto& example : examples) {
            const auto roots = geode::RootSolver<>(example.second).compute_roots(tolerance);
            std::cout << example.first << ": " << roots.size() << " distinct real root(s)\n";
            for (const auto& root : roots) {
                std::cout << "  " << root.to_string()
                          << "  |P(x)| = " << geode::abs(example.second(root)).to_string() << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Cannot compute roots: " << error.what() << '\n';
        return 1;
    }
}
