#!/usr/bin/env bash
## @file Build the browser solver and package the static GitHub Pages site.
set -euo pipefail

geode_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
geode_cache=${GEODE_WEB_CACHE:-"$geode_root/out/web-cache"}
geode_output="$geode_root/out/web"
geode_jobs=${GEODE_BUILD_JOBS:-4}
geode_version=4.0.23

if ! command -v emcc >/dev/null || ! emcc --version | head -n 1 | grep -F " $geode_version " >/dev/null; then
    echo "Activate Emscripten $geode_version with emsdk before building." >&2
    exit 1
fi
mkdir -p "$geode_cache" "$geode_output/source"
geode_cache=$(cd "$geode_cache" && pwd)
geode_prefix="$geode_cache/emscripten-$geode_version"
mkdir -p "$geode_cache/downloads" "$geode_cache/sources" "$geode_prefix"

## @brief Download and verify a pinned upstream source archive before extracting it.
fetch_source() {
    local package=$1 url=$2 checksum=$3 archive="$geode_cache/downloads/$1.tar.xz"
    if [[ ! -f "$archive" && -f "$geode_root/downloads/$package.tar.xz" ]]; then
        cp "$geode_root/downloads/$package.tar.xz" "$archive"
    fi
    if [[ ! -f "$archive" ]]; then
        curl --fail --location --retry 3 "$url" --output "$archive.part"
        mv "$archive.part" "$archive"
    fi
    echo "$checksum  $archive" | sha256sum --check --status
    if [[ ! -d "$geode_cache/sources/$package" ]]; then
        mkdir -p "$geode_cache/sources/$package"
        tar -xJf "$archive" --strip-components=1 -C "$geode_cache/sources/$package"
    fi
}

fetch_source gmp-6.3.0 https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz \
    a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
fetch_source mpfr-4.2.2 https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.2.tar.xz \
    b67ba0383ef7e8a8563734e2e889ef5ec3c3b898a01d00fa0a6869ad81c6ce01

if [[ ! -f "$geode_prefix/lib/libgmpxx.a" ]]; then
    mkdir -p "$geode_cache/build-gmp"
    (
        cd "$geode_cache/build-gmp"
        emconfigure "$geode_cache/sources/gmp-6.3.0/configure" \
            --host=wasm32-unknown-emscripten --build="$("$geode_cache/sources/gmp-6.3.0/config.guess")" \
            --prefix="$geode_prefix" --disable-assembly --enable-cxx --disable-shared --enable-static \
            CFLAGS=-O2 CXXFLAGS='-O2 -fexceptions'
        emmake make -j"$geode_jobs"
        emmake make install
    ) >"$geode_cache/gmp-build.log" 2>&1 || { tail -n 60 "$geode_cache/gmp-build.log"; exit 1; }
fi

if [[ ! -f "$geode_prefix/lib/libmpfr.a" ]]; then
    mkdir -p "$geode_cache/build-mpfr"
    (
        cd "$geode_cache/build-mpfr"
        emconfigure "$geode_cache/sources/mpfr-4.2.2/configure" \
            --host=wasm32-unknown-emscripten --build="$("$geode_cache/sources/mpfr-4.2.2/config.guess")" \
            --prefix="$geode_prefix" --with-gmp="$geode_prefix" \
            --disable-shared --enable-static --disable-thread-safe --disable-decimal-float \
            --disable-float128 CFLAGS=-O2
        emmake make -j"$geode_jobs"
        emmake make install
    ) >"$geode_cache/mpfr-build.log" 2>&1 || { tail -n 60 "$geode_cache/mpfr-build.log"; exit 1; }
fi

em++ "$geode_root/web/solver.cpp" -I"$geode_root/include" -isystem "$geode_prefix/include" \
    "$geode_prefix/lib/libmpfr.a" "$geode_prefix/lib/libgmpxx.a" "$geode_prefix/lib/libgmp.a" \
    -std=c++17 -O2 -fexceptions -Wall -Wextra -Wpedantic \
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createGeode \
    -sENVIRONMENT=web,worker,node -sFILESYSTEM=0 -sDYNAMIC_EXECUTION=0 \
    -sEXPORTED_FUNCTIONS=_geode_solve -sEXPORTED_RUNTIME_METHODS=ccall \
    -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=16777216 -sMAXIMUM_MEMORY=134217728 \
    -sSTACK_SIZE=2097152 -sABORTING_MALLOC=0 -o "$geode_output/solver.mjs"

cp "$geode_root"/web/*.html "$geode_root"/web/*.css "$geode_root"/web/*.js \
    "$geode_root"/web/*.svg "$geode_output/"
touch "$geode_output/.nojekyll"
cp "$geode_root/LICENSE" "$geode_output/LICENSE"
{
    printf 'Geode Root Solver — Apache License 2.0. See LICENSE.\n\n'
    printf 'GMP 6.3.0 and MPFR 4.2.2 are distributed here under LGPLv3-or-later.\n'
    printf 'Complete corresponding application/library sources and build instructions:\n'
    printf 'source/geode-browser-source.tar.gz\n\n'
    cat "$geode_cache/sources/gmp-6.3.0/COPYING.LESSERv3"
    cat "$geode_cache/sources/gmp-6.3.0/COPYINGv3"
    printf '\nEmscripten runtime notices:\n'
    cat "$(dirname "$(command -v emcc)")/LICENSE"
} > "$geode_output/THIRD_PARTY_LICENSES.txt"

# Include the application source as well as exact dependency archives for rebuilding/relinking.
tar -czf "$geode_output/source/geode-browser-source.tar.gz" \
    -C "$geode_root" include web scripts README.md LICENSE package.json package-lock.json \
    CMakeLists.txt CMakePresets.json vcpkg.json Doxyfile examples tests playwright.config.mjs docs/WEB_SECURITY.md \
    -C "$geode_cache" downloads/gmp-6.3.0.tar.xz downloads/mpfr-4.2.2.tar.xz
echo "Static site built in $geode_output"
