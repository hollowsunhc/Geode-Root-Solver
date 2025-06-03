cmake_minimum_required(VERSION 3.22)
project(geode-root-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------------------
# Header-only interface library

add_library(geode_solver INTERFACE)
add_library(geode::solver ALIAS geode_solver)

target_include_directories(geode_solver
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# -------------------------------------------------------------------
# Example program
add_executable(demo examples/main.cpp)

# demo only needs to link against the solver target; the GMP/MPFR
# libraries and include paths propagate automatically.
target_link_libraries(demo PRIVATE geode::solver)

# -------------------------------------------------------------------
# Dependencies

if(UNIX)
    set(VCPKG_TRIPLET_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-linux")
    list(APPEND CMAKE_MODULE_PATH "${VCPKG_TRIPLET_DIR}")

    find_package(PkgConfig REQUIRED)

    # GMP (integer) -------------------------------------------------
    pkg_check_modules(gmp   REQUIRED IMPORTED_TARGET gmp)
    add_library(GMP::gmp   ALIAS PkgConfig::gmp)

    # GMPXX (C++ wrapper) ------------------------------------------
    pkg_check_modules(gmpxx REQUIRED IMPORTED_TARGET gmpxx)
    add_library(GMP::gmpxx ALIAS PkgConfig::gmpxx)

    # MPFR (float / interval) --------------------------------------
    pkg_check_modules(mpfr  REQUIRED IMPORTED_TARGET mpfr)
    add_library(MPFR::mpfr ALIAS PkgConfig::mpfr)

    # For std::execution::par
    find_package(Threads REQUIRED)

endif()

if(WIN32)
    set(VCPKG_TRIPLET_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed/x64-windows")
    set(GMP_INCLUDE_DIR   "${VCPKG_TRIPLET_DIR}/include")
    set(LIB_DIR           "${VCPKG_TRIPLET_DIR}/lib")

    # GMP (integer) -------------------------------------------------
    find_library(GMP_LIBRARY     NAMES gmp     PATHS "${LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
    add_library(GMP::gmp STATIC IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
        IMPORTED_LOCATION            "${GMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
    )

    # GMPXX (C++ wrapper) ------------------------------------------
    find_library(GMPXX_LIBRARY   NAMES gmpxx   PATHS "${LIB_DIR}" NO_DEFAULT_PATH)
    if (GMPXX_LIBRARY)
        add_library(GMP::gmpxx STATIC IMPORTED)
        set_target_properties(GMP::gmpxx PROPERTIES
            IMPORTED_LOCATION            "${GMPXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
        )
    else()
        add_library(GMP::gmpxx ALIAS GMP::gmp)
    endif()

    # MPFR (float / interval) --------------------------------------
    find_library(MPFR_LIBRARY    NAMES mpfr    PATHS "${LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
    add_library(MPFR::mpfr STATIC IMPORTED)
    set_target_properties(MPFR::mpfr PROPERTIES
        IMPORTED_LOCATION            "${MPFR_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"  # mpfr headers live in same include dir
    )
endif()

# The solver relies on GMP / MPFR symbols at link time (mpreal uses them)
target_link_libraries(geode_solver
    INTERFACE
        GMP::gmp
        GMP::gmpxx
        MPFR::mpfr
)
