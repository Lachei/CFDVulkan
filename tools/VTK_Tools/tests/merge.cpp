// 010-TestCase.cpp
#include "catch.hpp"
#include <iostream>
#include <filesystem>
#include "../helper.hpp"
#include "fmt/format.h"

TEST_CASE("Extract Rank from file name, manual definition", "[manual][preparation]") {
    // Variable Definition
    int rank;
    std::filesystem::path file_path;
    std::string file_name;
    // Sample Testing
    file_path = std::filesystem::path("cavity_0.1.vtk");
    rank = get_rank(file_path);
    REQUIRE(rank == 0);
    file_path = std::filesystem::path("cavity_1.1.vtk");
    rank = get_rank(file_path);
    REQUIRE(rank == 1);
}

TEST_CASE("Extract Rank from file name, auto definition", "[auto][preparation]") {
    // Variable Definition
    int rank;
    std::filesystem::path file_path;
    std::string file_name;
    // Auto Generate different test case
    int i = GENERATE(range(0, 100));
    file_name = fmt::format("cavity_{0}.1.vtk", i);
    file_path = std::filesystem::path(file_name);
    rank = get_rank(file_path);
    REQUIRE(rank == i);
}

TEST_CASE("Extract Timestep from file name, manual definition", "[manual][preparation]") {
    // Variable Definition
    int ts;
    std::filesystem::path file_path;
    std::string file_name;
    // Sample Testing
    file_path = std::filesystem::path("cavity_0.1.vtk");
    ts = get_timestamp(file_path);
    REQUIRE(ts == 1);
    file_path = std::filesystem::path("cavity_1.10.vtk");
    ts = get_timestamp(file_path);
    REQUIRE(ts == 10);
}

TEST_CASE("Extract Timestep from file name, auto definition", "[auto][preparation]") {
    // Variable Definition
    int ts;
    std::filesystem::path file_path;
    std::string file_name;
    // Auto Generate different test case
    int i = GENERATE(range(0, 100));
    file_name = fmt::format("cavity_0.{0}.vtk", i);
    file_path = std::filesystem::path(file_name);
    ts = get_timestamp(file_path);
    REQUIRE(ts == i);
}

