cmake_minimum_required(VERSION 3.15)
project(tket-c-api CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# TKET dependency
find_package(tket REQUIRED)

# Library and headers to be installed
add_library(tket-c-api src/tket-c-api.cpp)
target_include_directories(tket-c-api PUBLIC include)

target_link_libraries(tket-c-api PRIVATE tket::tket)

set_target_properties(tket-c-api PROPERTIES PUBLIC_HEADER "include/tket-c-api.h")
install(TARGETS tket-c-api)
