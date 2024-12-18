cmake_minimum_required(VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------#
# PROJECT DEPENDENCY
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#----------------------------------------------------------------------#
# YAML DEPENDENCY
add_subdirectory("yaml-cpp") #NOTE THIS IS A RELATIVE PATH

#----------------------------------------------------------------------#
# ARGPARSE DEPENDENCY
#set(ARGPARSE_TEST_ENABLE "OFF") #NOTE: THIS IS NEEDED AT LEAST ON MACOS
#add_subdirectory("argparse") #NOTE: THIS IS A RELATIVE PATH

#----------------------------------------------------------------------#
# CERN ROOT DEPENDENCY
find_package(ROOT CONFIG COMPONENTS RIO Net RooStats RooFit ROOTDataFrame)
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------#
# ADD LIBRARIES
file(GLOB SRC_HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB SRC_LIB_FILES ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp)
message(STATUS "SRC_LIB_FILES = ${SRC_LIB_FILES}")#DEBUGGING
message(STATUS "SRC_HEADER_FILES = ${SRC_LIB_FILES}")#DEBUGGING
add_library (saga ${SRC_HEADER_FILES} ${SRC_LIB_FILES}) #NOTE: Need both header and src files here even if empty.
target_include_directories(saga PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

#----------------------------------------------------------------------#
# ADD EXECUTABLES

# Set executables directory
set(EXECUTABLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(getKinBinnedAsym1D ${EXECUTABLES_DIR}/getKinBinnedAsym1D.cpp)
target_link_libraries(getKinBinnedAsym1D PUBLIC ROOT::RIO ROOT::Net ROOT::RooFit ROOT::RooStats ROOT::ROOTDataFrame yaml-cpp)
