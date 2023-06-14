
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(VTK_Tools VERSION 1.0)

# C++ Standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# VTK Library
find_package(VTK REQUIRED)
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

# FMT Library
find_package(fmt REQUIRED)

# Build executables
# VTK Compares
add_executable(VTK_Comp Compare.cpp Reader.cpp FileIterator.cpp)
target_link_libraries(VTK_Comp PRIVATE ${VTK_LIBRARIES})

# VTK Merge
add_executable(VTK_Merge Merge.cpp FileIterator.cpp)
target_link_libraries(VTK_Merge PRIVATE ${VTK_LIBRARIES} fmt::fmt)

# Testing
enable_testing ()
add_subdirectory(tests)