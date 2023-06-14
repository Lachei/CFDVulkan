# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.0)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(shader)

# You can find package likes
# find_package(MPI)
# Require a package
find_package(MPI REQUIRED)
find_package(Vulkan REQUIRED)
# Find a package with different components e.g. BOOST
# find_package(Boost COMPONENTS filesystem REQUIRED)

# VTK Library
find_package(VTK REQUIRED)
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

# Filesystem
#set(gpp9 False)
set(gpp9 True)
if(gpp9)
    message("g++ Version is 9 or higher")
    add_compile_definitions(gpp9)
else()
    message("g++ Version is lower than Version 9")
    # Link Filesystem Libraries
    link_libraries(stdc++fs)
endif()

# Creating the executable of our project and the required dependencies
# the executable is called sim
file(GLOB sim_SRC
    "*.h"
    "*.cpp"
	"*.hpp"
)

add_executable(sim ${sim_SRC})
# if you use external libraries you have to link them like 
# target_link_libraries(sim PRIVATE MPI::MPI_CXX Boost::filesystem)
target_link_libraries(sim PRIVATE ${VTK_LIBRARIES} MPI::MPI_CXX Vulkan::Vulkan)

# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
enable_testing ()
add_subdirectory(tests)

file(GLOB CPY
  "simulation_settings/*.dat"
  "simulation_settings/*.pgm"
)
file(COPY ${CPY} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/simulation_settings)

# creating output folder for .vtk images
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)