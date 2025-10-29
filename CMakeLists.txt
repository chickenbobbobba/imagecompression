cmake_minimum_required(VERSION 3.11.0)
project(imagecompression VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
cmake_policy(SET CMP0072 NEW)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MAKE_PROGRAM imagecompression)

add_compile_options(-Wall -Wextra -Wpedantic -O3 -march=native -ffast-math -funsafe-math-optimizations)

file(GLOB_RECURSE APP_SOURCES "${CMAKE_SOURCE_DIR}/src/*")

add_executable(imagecompression ${APP_SOURCES})
target_include_directories(imagecompression PRIVATE ${CMAKE_SOURCE_DIR}/include)
