# valgrind --tool=callgrind
# kcachegrind callgrind.out.11463 
cmake_minimum_required(VERSION 2.8)
project(glores)
add_definitions(-std=c++11 -pg -fPIC) # -pg for gprof
set(CMAKE_BUILD_TYPE Release)         # RelWithDebug

#set(MATLAB_ADDITIONAL_VERSIONS 
#"R2018b=9.5")

message(STATUS "Added dir of CMake Modules: ${glores_SOURCE_DIR}/cmake_modules")
LIST(APPEND CMAKE_MODULE_PATH "${glores_SOURCE_DIR}/cmake_modules") 

# Library Eigen3
find_package(Eigen3 REQUIRED)
if (${EIGEN3_FOUND})
  include_directories(${EIGEN3_INCLUDE_DIR})
  message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
else(${EIGEN3_FOUND})
  message(WARNING "Cannot find Eigen3 Library")
endif(${EIGEN3_FOUND})

# Boost 
find_package(Boost COMPONENTS system filesystem iostreams  REQUIRED)
if (${Boost_FOUND})
  include_directories(${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY})
endif (${Boost_FOUND})
message(STATUS "Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES} ") 


set(MATLAB_ADDITIONAL_VERSIONS "R2018a=9.4")
find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY )


include_directories(include)

# Library files
add_library(glores
  src/fileutils.cpp
  src/functions.cpp
  src/geometry.cpp
  src/ParamMap.cpp
  src/PointSetGlobalRegistration.cpp
  src/Profiler.cpp
)

# Executables
add_executable(testPointSetGlobalRegistration app/testPointSetGlobalRegistration.cpp)
target_link_libraries(testPointSetGlobalRegistration glores ${Boost_LIBRARIES})

add_executable(testCheapBoundQueue app/testCheapBoundQueue.cpp)
target_link_libraries(testCheapBoundQueue glores ${Boost_LIBRARIES})

add_executable(testGeometryIntersection app/testGeometryIntersection.cpp)
target_link_libraries(testGeometryIntersection glores ${Boost_LIBRARIES})

add_executable(testGeometryBranch app/testGeometryBranch.cpp)
target_link_libraries(testGeometryBranch glores ${Boost_LIBRARIES})



matlab_add_mex(
    NAME glores_mex ${CMAKE_SOURCE_DIR}/app/glores_mex.cpp LINK_TO glores ${Boost_LIBRARIES}
)

# Required to solve a bug on findmatlab in version r2018b
target_compile_options(glores_mex PRIVATE "-fvisibility=default")







