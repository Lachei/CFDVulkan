cmake_minimum_required(VERSION 3.10)
include_directories (${CFDLAB_SOURCE_DIR}) 

# Test Executables
add_executable(test_factorial sample_test.cpp ) 

# if you want to copy config files  Config Files
# configure_file(${CFDLAB_SOURCE_DIR}/source.dat ${CMAKE_BINARY_DIR}/tests/target_example.dat COPYONLY)

# Testing
add_test (NAME example_test COMMAND test_factorial)
# Example for MPI RUN
# add_test (NAME example_test COMMAND mpirun -np 4 test [2x2])
