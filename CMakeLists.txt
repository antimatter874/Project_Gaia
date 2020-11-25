cmake_minimum_required(VERSION 3.17)
project(Project_Gaia)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} main.cpp impl/maths.hpp)

find_package(Threads)
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
