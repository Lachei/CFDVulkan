cmake_minimum_required(VERSION 3.10)
include_directories (${VTK_Tools_SOURCE_DIR}) 

# Required files fot testing
add_library(lib_base STATIC main.cpp ${VTK_Tools_SOURCE_DIR}/helper.hpp)
find_library(fmt REQUIRED)

# Test Executables
add_executable(test_merge merge.cpp ) 
target_link_libraries(test_merge lib_base fmt::fmt)

# Testing 
add_test (NAME Test_Merge COMMAND test_merge)