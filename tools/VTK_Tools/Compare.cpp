//
// Created by Moritz Gnisia on 30.01.20.
// Useful resources:
// https://github.com/liaoyg/vtkUnstructuredGridReader
// https://lorensen.github.io/VTKExamples/site/Cxx/IO/ReadStructuredGrid/
// mkdir build
// cd build
// cmake ..
// ..

#include <string>
#include <iostream>
#include <filesystem>
#include <vtkType.h>
#include <array>
#include <vtkTuple.h>
#include "Reader.hpp"
#include "FileIterator.hpp"
#include "cxxopts.hpp"
#include <vtkSmartPointer.h>

// #define FAILURE
namespace fs = std::filesystem;


int main(int argc, char *argv[])
{
    // Definition of CLI Options
    cxxopts::Options options("VTK Compare Tool", "Comparing two folders of VTK Files");

    // Inputvariables
    std::string first;
    std::string second;
    std::string origin;
    std::string current;
    
    options.add_options()
        ("o,origin", "Origin folder name of the vtk files e.g. folder/; this is the reference solution", cxxopts::value<std::string>())
        ("c,current", "Folder name which contains the vtk files that should be compared e.g. dest/", cxxopts::value<std::string>())
        ("f,first", "prefix of the vtk files in the origin (first) given folder e.g. cavity", cxxopts::value<std::string>())
        ("s,second", "prefix of the vtk files in the current (second) given folder e.g. cavity", cxxopts::value<std::string>())
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(1);
    }
    if (!result.count("origin")) {
        std::cout << "Please provide the origin location of the vtk files" << std::endl;
        exit(1);
    }
    else {
        origin = result["origin"].as<std::string>();
    }
    if (!result.count("current")) {
        std::cout << "Please provide the current location of the vtk files" << std::endl;
        exit(1);
    }
    else {
        current = result["current"].as<std::string>();
    }
    if (!result.count("first")) {
        std::cout << "Please provide the prefix of the vtk files in the origin folder" << std::endl;
        exit(1);
    }
    else {
        first = result["first"].as<std::string>();
    }
    if (!result.count("second")) {
        std::cout << "Please provide the prefix of the vtk files in the second folder" << std::endl;
        exit(1);
    }
    else {
        second = result["second"].as<std::string>();
    }

    // string array for the true solution
    std::string dirs[2] = {origin, current};
    std::string prefix[2] = {first, second};

    // Creating two directory iterators
    std::filesystem::directory_iterator true_sol_iter(dirs[0]); // Comparison Directory
    std::filesystem::directory_iterator curr_dir_iter(dirs[1]); // Current Solution Directory

    // File Iterator
    FileIterator file_iter_true(true_sol_iter, prefix[0], false);
    FileIterator file_iter_curr(curr_dir_iter, prefix[1], false);

    // Get Iterators
    vector_path* vec_true = file_iter_true.get_dir_vec();
    vector_path* vec_curr = file_iter_curr.get_dir_vec();

    // Comparison of the amount of files in the different directories
    /*bool comp = file_iter_true.compare(*vec_curr);
    if (!comp) {
        std::cout << "The amount of files doesnt match in the directories, the tool will exit now" << std::endl;
        return EXIT_FAILURE;
    }
*/
    bool failure(false);
    // Looping simultaneously over both iterators, only compare the files if there name is identical!
    // https://stackoverflow.com/questions/9602918/c-initialising-2-different-iterators-in-a-for-loop
    // https://stackoverflow.com/questions/30983154/get-an-ordered-list-of-files-in-a-folder
    for (std::pair<vector_path::iterator, vector_path::iterator> files(vec_true->begin(), vec_curr->begin());
     files.first != vec_true->end(); ++files.first, ++files.second) {
                #ifdef LOG
                    std::cout << "File True: " << files.first->string() << "\t" << 
                    "File to compare: " << files.second->string() << std::endl;
                #endif
       
                Reader true_solution = Reader(files.first->string());
                Reader curr_solution = Reader(files.second->string());
                // curr_solution = read_cell_point(b.path());
                vtkIdType cells_c = curr_solution.sg->GetNumberOfCells();
                vtkIdType cells_t = true_solution.sg->GetNumberOfCells();


                // Comparing the amount of grid cells
                if (cells_t != cells_c) {
                    std::cerr << "The current solution " << files.second->string() <<
                        " has not same amount of values for the grid as the true solution." << std::endl;
                    failure = true;
                }

                // Comparing the grid values
                vtkTuple<double, 3> p_t;
                vtkTuple<double, 3> p_c;
                for (vtkIdType j = 0; j < cells_t; ++j) {
                    true_solution.sg->GetPoint(j, p_t.GetData());
                    curr_solution.sg->GetPoint(j, p_c.GetData());
                    if (!( p_t.Compare(p_c, 1e-4))) {
                        std::cerr << "The current solution " << files.second->string() <<
                                  " has not identical values for the grid at position " << j << std::endl;
                        failure = true;
                    }
                }

                // Comparing Scalar Values
                for (int j = 0; j < true_solution.cd->GetNumberOfArrays() ; ++j) {

                    // Scalar Name of Array e.g. pressure, geometry
                    const char* scalar_name = true_solution.cd->GetArrayName(j);
                    #ifdef LOG
                        std::cout << "\t Processing " << scalar_name << std::endl;
                    #endif
                    // Compare if the current solution has the same values
                    if (!curr_solution.cd->HasArray(scalar_name)) {
                        std::cerr << "The current solution " << files.second->string() <<
                                  " doesn't contain the same array as in the true solution" << files.first->string()
                                    << scalar_name << std::endl;
                        failure = true;
                    }

                    // Assigning the values to data
                    vtkSmartPointer<vtkDataArray> arr_t = true_solution.cd->GetScalars(scalar_name);
                    vtkSmartPointer<vtkDataArray> arr_c = curr_solution.cd->GetScalars(scalar_name);

                    // Compare Number of values
                    if (arr_t->GetNumberOfValues() != arr_c->GetNumberOfValues()) {
                        std::cerr   << "The current solution has "  << files.second->string() <<
                                    "not same amount of values as the true solution." << std::endl;
                        failure = true; 
                    }


                    // Compare the values itself
                    vtkTuple<double, 1> p_t;
                    vtkTuple<double, 1> p_c;
                    for (vtkIdType i = 0; i < arr_t->GetNumberOfTuples(); ++i) {
                        arr_t->GetTuple(i, p_t.GetData());
                        arr_c->GetTuple(i, p_c.GetData());
                        if (!( p_t.Compare(p_c, 1e-6))) {
                            std::cerr << "The current solution " << files.second->string() <<
                                      " has not identical values for the grid at position " << i << files.first->string() <<
                                      " for array " << scalar_name << "\nOriginal Solution: "  << p_t[0] << "\nCurrent Solution: " << p_c[0] << std::endl;
                            failure = true;
                        }
                    }
                }
                // Comparing Vector Values
                for (int j = 0; j < true_solution.pd->GetNumberOfArrays() ; ++j) {

                    // Scalar Name of Array e.g. velocity
                    const char* vec_name = true_solution.pd->GetArrayName(j);
                    #ifdef LOG
                        std::cout << "\t Processing " << vec_name << std::endl;
                    #endif
                    // Compare if the current solution has the same values
                    if (!curr_solution.pd->HasArray(vec_name)) {
                        std::cerr << "The current solution " << files.second->string() <<
                                  " doesn't contain the same array as in the true solution"
                                    << vec_name << std::endl;
                        failure = true;
                    }

                    // Assigning the values to data
                    vtkSmartPointer<vtkDataArray> arr_t = true_solution.pd->GetScalars(vec_name);
                    vtkSmartPointer<vtkDataArray> arr_c = curr_solution.pd->GetScalars(vec_name);

                    // Compare Number of values
                    if (arr_t->GetNumberOfValues() != arr_c->GetNumberOfValues()) {
                        std::cerr   << "The current solution has "  << files.second->string() <<
                                    "not same amount of values as the true solution." << std::endl;
                        failure = true;
                    }

                    // Compare the values (vectors) itself
                    vtkTuple<double, 3> p_t;
                    vtkTuple<double, 3> p_c;
                    for (vtkIdType i = 0; i < arr_t->GetNumberOfTuples(); ++i) {
                        arr_t->GetTuple(i, p_t.GetData());
                        arr_c->GetTuple(i, p_c.GetData());
                        if (!( p_t.Compare(p_c, 1e-6))) {
                            std::cerr << "The current solution " << files.second->string() <<
                                      " has not identical values for the grid at position " << i <<
                                      " for array " << vec_name << 
                                      "\nOriginal Solution: "  << p_t[0] << " " << p_t[1] << " " << p_t[2] <<
                                      "\nCurrent Solution: "  << p_c[0] << " " << p_c[1] << " " << p_c[2] << std::endl;
                            failure = true;
                        }
                    }

                }

        }
    if (!failure) {
        std::cout << "No errors detected" << std::endl;
        return EXIT_SUCCESS;
    }
    else {
        std::cerr << "Erros detected" << std::endl;
        return EXIT_FAILURE;
    }
}