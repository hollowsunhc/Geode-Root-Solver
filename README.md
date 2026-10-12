# Geode Root Solver

A C++17, header-only library that finds **all distinct real roots** of a polynomial.
It uses GMP for exact arithmetic and MPFR for values with adjustable precision.
The solver combines Geode refinement with exact root isolation. It returns roots
in ascending order, or throws an exception if it cannot meet the requested accuracy.

The library handles zero roots, repeated roots, missing linear terms, and
polynomials with no real roots. It does not return complex roots or multiplicities.

## Interactive browser lab

Try the [Geode polynomial lab](https://hollowsunhc.github.io/Geode-Root-Solver/)
to explore examples, plot a polynomial, and find its real roots without installing
anything. The same C++ solver runs as WebAssembly in a background worker. Inputs
and results stay in your browser. The graph uses approximate coordinates; expand
a result to see an exact decimal encoding of the computed MPFR approximation.
The requested absolute error bound determines its accuracy, not the number of
printed digits.

The lab includes examples of crossings, repeated roots, absent real roots, and
nearby roots, with explanations of precision and error bounds. Browser calculations
are limited to degree 16, 512 working bits, five seconds, and 128 MiB of Wasm linear
memory. Cancel stops the worker. See [browser security](docs/WEB_SECURITY.md) for
the full limits and their scope.

### Build and verify the website

The web build is separate from the native CMake build. On Linux or WSL, install
Git, Python 3, Make, a native C/C++ compiler, curl, tar, xz, and Node.js 22 or later.
Activate the pinned [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html):

```sh
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install 4.0.23
./emsdk activate 4.0.23
source ./emsdk_env.sh
```

From this repository, run:

```sh
npm ci
npm run build:web
npm test
npx playwright install --with-deps chromium
npm run test:browser
npm run verify:security
python3 -m http.server 4173 --bind 127.0.0.1 --directory out/web
```

Open `http://127.0.0.1:4173`. The first build compiles pinned GMP and MPFR sources;
later builds reuse them. `GEODE_WEB_CACHE` can move the dependency cache from
`out/web-cache`, and `GEODE_BUILD_JOBS` sets its parallel build count (default 4).
The build checks source archive hashes and writes the complete static site to
`out/web`, including license notices and a downloadable source bundle.

Publish the contents of `out/web` to the root of the `gh-pages` branch and select
that branch in **Settings → Pages → Deploy from a branch**. Keep `.nojekyll` and
the `source` directory. All asset paths are relative, including the worker and
Wasm module. No custom domain, backend, or CI workflow is required.

## Build and test

Use CMake 3.23 or later, Ninja, Git, and a C++17 compiler. The presets use vcpkg
manifest mode. `vcpkg.json` pins the dependency versions to a registry baseline.
Set `VCPKG_ROOT` to the directory that contains `vcpkg` and its `scripts` directory.
CMake installs GMP and MPFR when you configure the project.
The first configure step can be slow because vcpkg builds these libraries from source.
If an older vcpkg checkout cannot find the pinned versions, update that checkout
and run its bootstrap script again before you configure.

### Linux

On Ubuntu or Debian, install the build tools:

```sh
sudo apt-get install build-essential cmake ninja-build git curl zip unzip tar \
    pkg-config autoconf automake autoconf-archive libtool
```

If you do not have vcpkg, install it:

```sh
git clone https://github.com/microsoft/vcpkg.git "$HOME/vcpkg"
"$HOME/vcpkg/bootstrap-vcpkg.sh"
export VCPKG_ROOT="$HOME/vcpkg"
```

From this project directory, run:

```sh
cmake --preset linux-release
cmake --build --preset linux-release
ctest --preset linux-release
./out/build/linux-release/demo
```

Use `linux-debug` in all three commands for a Debug build.

### Windows

Install Visual Studio 2022 Build Tools with **Desktop development with C++**,
including CMake and Ninja. Install Git and
[vcpkg](https://learn.microsoft.com/vcpkg/get_started/get-started).
Open an **x64 Native Tools Command Prompt for VS 2022**. Set `VCPKG_ROOT` to your
vcpkg directory, then run these commands from this project directory:

```bat
set VCPKG_ROOT=C:\path\to\vcpkg
cmake --preset msvc-vcpkg-release
cmake --build --preset msvc-vcpkg-release
ctest --preset msvc-vcpkg-release
out\build\msvc-vcpkg-release\demo.exe
```

Use `msvc-vcpkg-debug` in all three commands for a Debug build. Each preset has a
separate build directory. The vcpkg toolchain selects the dependency triplet;
the project does not contain fixed dependency paths.

### VS Code

Install the recommended **C/C++** and **CMake Tools** extensions. Start VS Code
from a shell that has `VCPKG_ROOT` set. On Windows, use the Native Tools prompt.
On Linux, install GDB if you want to debug the example.

1. Open the project folder.
2. Run **CMake: Select Configure Preset** and select a Debug preset.
3. Run **CMake: Configure**, then **CMake: Build**.
4. Run **CMake: Run Tests**.
5. Select `demo` as the launch target. Select the debugger for your system and press F5.

The repository includes configure, build, and test presets, debugger settings,
and a `.clang-format` file. It does not include CI.

### Use system libraries

You can also use installed GMP and MPFR development packages without vcpkg:

```sh
cmake -S . -B out/build/system -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build out/build/system
ctest --test-dir out/build/system --output-on-failure
```

## Use the library

Coefficients are in **ascending** order. For example, `{-6, 11, -6, 1}` means
`x^3 - 6x^2 + 11x - 6`.

```cpp
#include "geode_root_solver.hpp"
#include <iostream>

int main() {
    using Real = geode::HighPrecFloat;
    geode::set_global_bits(256);

    const geode::Polynomial<Real> polynomial({-6, 11, -6, 1});
    const geode::RootSolver<> solver(polynomial);
    for (const Real& root : solver.compute_roots(Real("1e-30"))) {
        std::cout << root.to_string() << '\n'; // 1, 2, and 3
    }
}
```

To use the CMake target in another project:

```cmake
add_subdirectory(path/to/Geode-Root-Solver)
target_link_libraries(your_program PRIVATE geode::solver)
```

Configure the parent project with the vcpkg toolchain and list `gmp` and `mpfr`
in the parent's manifest, or install their development packages. The target
supplies the include paths, link libraries, and C++17 requirement. Examples and
tests are off by default when this project is a subdirectory.

To use the header directly, add `include` to your include path and link MPFR,
GMPXX, and GMP, in that order.

### API documentation

The header uses Doxygen comments for the public API, including precision rules,
parameters, return values, and exceptions. If Doxygen is installed, run:

```sh
doxygen Doxyfile
```

Open `out/docs/html/index.html` to read the generated documentation. Doxygen is
optional and is not needed to build or use the solver.

## Accuracy and errors

`compute_roots(tolerance)` uses an **absolute root error**, not a residual test.
Each returned value is within `tolerance` of a distinct real root. The default
tolerance is `1e-30`. A small value of `|P(x)|` alone does not prove this accuracy,
especially near a repeated root.

The guarantee applies to the **stored coefficients**. MPFR stores binary numbers.
A decimal string can be rounded when you construct a coefficient. Very close or
repeated roots can change when the coefficients change. Set enough precision
**before** you construct the coefficients. Use decimal strings for values that
would lose digits in a C++ `double` literal.

| Input or condition | Result |
| --- | --- |
| Nonzero constant polynomial | Empty root vector |
| Polynomial with no real roots | Empty root vector |
| Repeated real root | One result for that root |
| Empty coefficient list or all zero coefficients | `std::invalid_argument` from the solver |
| NaN or infinite coefficient | `std::invalid_argument` |
| Zero, negative, NaN, or infinite tolerance | `std::invalid_argument` |
| Too few bits to meet accuracy or distinguish roots | `std::runtime_error`; increase precision |
| Invalid number string or base | `std::invalid_argument` |

The default precision is 256 bits, about 77 decimal digits. `set_global_bits` and
`get_global_bits` control the default **for the current thread**. Existing numbers
retain their precision. Arithmetic uses at least the precision of both operands.
`RootSolver` supports `HighPrecFloat`; it does not support `double` as its template
argument.

`Ball` checks membership in an interval with a midpoint and radius. It is not a
complete interval-arithmetic package. The `PI` constant has 256-bit
precision; use `geode::acos(Real(-1))` to compute pi at the selected precision.

## Mathematics and performance

The root series uses **hyper-Catalan coefficients** `C`, not Geode coefficients
`G`. For an index `m = (m2, m3, ...)`,

```text
C[m] = (sum k*mk)! / ((1 + sum (k-1)*mk)! * product mk!)
S = sum C[m] * product tk^mk
S = 1 + (t2 + t3 + ...) G
```

The solver combines the following methods:

- Rubine's recurrence computes `G` from `C`. `GeodeTable` calculates
  coefficients on demand and caches them. `hyper_catalan` returns
  exact rational values. A table index lists `m2, m3, ...`; its first constructor
  argument is the number of these variables. Do not share one mutable table
  between threads without a lock.
- Gessel's identity evaluates the Geode near a root:
  `G = 1 / (1 - sum tk * (1 + S + ... + S^(k-1)))`.
  The solver shifts the polynomial to the current estimate and applies this
  identity twice to improve the Newton correction. It does not divide by
  `t2 + t3 + ...`, which can be zero.
- Exact polynomial division removes repeated factors. A Sturm sequence counts
  and isolates all real roots. Exact sign tests check each proposed correction.
  Bisection supplies progress when a correction is unsuitable. Thus, a formal
  series outside its convergence region cannot create an unchecked result.

Local shifts use quadratic work in the degree and linear storage. GMP integer
normalization reduces growth during exact polynomial division. Exact isolation
can be expensive for high degrees, large coefficients, or tightly clustered roots.

The tests cover the examples, repeated and absent real roots, close roots, large
scale differences, a degree-17 polynomial, 120-digit accuracy, invalid inputs,
precision changes, header linkage, and published Geode coefficients. Checks run
in both Debug and Release builds.

### References

- Wildberger and Rubine, [A Hyper-Catalan Series Solution to Polynomial Equations,
  and the Geode](https://doi.org/10.1080/00029890.2025.2460966), 2025.
- Rubine, [Hyper-Catalan and Geode Recurrences and Three Conjectures of
  Wildberger](https://arxiv.org/abs/2507.04552), 2025.
- Gessel, [Lattice paths and the Geode](https://arxiv.org/abs/2507.09405), 2025,
  Theorem 2.1.
- Rubine, [Computing the 4D Geode](https://arxiv.org/abs/2512.21785), December 2025.
  Its published small diagonal values provide independent test fixtures. Its
  large diagonal computation is separate from numerical root refinement.

## License

This project uses the [Apache License 2.0](LICENSE). GMP and MPFR have their own license
terms; see the dependency distributions when you distribute linked binaries.
