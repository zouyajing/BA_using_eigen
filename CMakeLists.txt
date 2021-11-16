cmake_minimum_required(VERSION 2.8)

project(BA_ceres_eigen)

find_package (Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})


add_executable(BA_ceres BA_ceres.cpp helper.cpp)
target_link_libraries(BA_ceres Eigen3::Eigen ${CERES_LIBRARIES})

add_executable(BA_eigen BA_eigen.cpp helper.cpp)
target_link_libraries(BA_eigen Eigen3::Eigen)

