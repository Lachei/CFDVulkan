//
// Created by Moritz Gnisia on 08.02.20.
//

#ifndef HELPER_H
#define HELPER_H

#include "filesystem"
#include "regex"
#include "FileIterator.hpp"

int get_timestamp(std::filesystem::path file) {
    std::smatch base_match;
    std::string file_name = file.string();
    std::regex_search(file_name, base_match, std::regex("\\b\\d+"));
    return  std::stoi(base_match[0]);
}

int get_rank(std::filesystem::path file) {
    std::smatch base_match;
    std::string file_name = file.string();
    std::regex_search(file_name, base_match, std::regex("\\B\\d+\\."));
    return  std::stoi(base_match[0]);
}

std::vector<int> get_unique(vector_path* files, int func (std::filesystem::path)) {
    std::vector<int> values;
    std::transform(files->begin(), files->end(), std::back_inserter(values), func);
    // Sorting, otherwise the unique array contains the same amount of values
    std::sort(values.begin(), values.end());

    // Get Unique Values
    auto last = std::unique(values.begin(), values.end());

    // Remove duplicates
    values.erase(last, values.end()); 
    return values;
}

std::vector<int> subset(std::vector<int> vec, int n) {
    std::vector<int> subset_vector;
    for (size_t i = 0; i < n; i++)
    {
        subset_vector.push_back(vec[i]);
    }
    return subset_vector;
}

std::string print_nomenclature() {
    return 
    "   Nomenclature \n \n"
    "   Rank <---+ +--> Timestep \n"
    "            | |            \n"
    "            | |            \n"
    "     cavity_1.10.vtk       \n"
    "     -------               \n"
    "        |                  \n"
    "        |                  \n"
    "        |                  \n"
    "     prefix  \n";
}

#endif