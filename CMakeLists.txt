cmake_minimum_required(VERSION 3.20)
project(photon_mapping LANGUAGES CXX)

option(BUILD_TESTS "build tests" OFF)

# spdlog
if(NOT TARGET spdlog)
    find_package(spdlog REQUIRED)
endif()

# OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(WARNING "failed to find OpenMP")
endif()

# pm(photon_mapping)
add_library(pm INTERFACE)
target_include_directories(pm INTERFACE "include")
target_compile_features(pm INTERFACE cxx_std_20)
set_target_properties(pm PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(pm INTERFACE spdlog::spdlog)
target_link_libraries(pm INTERFACE OpenMP::OpenMP_CXX)

# examples
add_subdirectory("examples")

# tests
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()