cmake_minimum_required(VERSION 3.20)
project(saga CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------#
# Set options
option(BUILD_DOXYGEN "Build doxygen documentation" TRUE)

#----------------------------------------------------------------------#
# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------#
# Actual Project
add_subdirectory ("saga")

#----------------------------------------------------------------------#
# Doxygen documentation
if (BUILD_DOXYGEN)
    add_subdirectory ("docs")
endif()
