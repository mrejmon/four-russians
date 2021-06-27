cmake_minimum_required(VERSION 3.5.1)
project(four_russians)

# march=native is set here, so that it can be applied
# to Edlib if Edlib is built from scratch
# march=native should only be set in RELEASE mode
add_compile_options("$<$<CONFIG:RELEASE>:-march=native>")

find_library(EDLIB_LOC edlib)
# if Edlib was not found, check edlib/ folder
if (NOT EDLIB_LOC)
  add_subdirectory(edlib EXCLUDE_FROM_ALL)
  set(EDLIB_LOC "edlib")
endif()

find_library(BENCHMARK_LOC benchmark)
# if Google Benchmark was not found, check benchmark/ folder
if (NOT BENCHMARK_LOC)
  # tells Google Benchmark to not build tests
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "")
  add_subdirectory(benchmark EXCLUDE_FROM_ALL)
  set(BENCHMARK_LOC "benchmark")
endif()

# set the required C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# find and set OpenMp compile options
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set the output directory to bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# add warnings
add_compile_options(-Wall -Wextra -pedantic)
# -03 is set by default by -D CMAKE_BUILD_TYPE=Release

# declare the executables and their libraries

add_executable(my_example src/example.cpp)
set_target_properties(my_example PROPERTIES OUTPUT_NAME "example")

add_executable(my_test src/test.cpp)
target_link_libraries(my_test ${EDLIB_LOC})
set_target_properties(my_test PROPERTIES OUTPUT_NAME "test")

add_executable(my_benchmark src/benchmark.cpp)

# linking with posix threads is recommended by Google Benchmark to avoid errors
target_compile_options(my_benchmark PRIVATE -pthread)

set_target_properties(my_benchmark PROPERTIES OUTPUT_NAME "benchmark")
target_link_libraries(my_benchmark ${EDLIB_LOC} ${BENCHMARK_LOC})
