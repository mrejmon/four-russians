This repository contains a mirror of the code for my bachelor's thesis:
## Multi-threaded implementation of the Four-Russians edit distance algorithm

### Abstract

Edit distance can be computed with the well-known dynamic programming algorithm in O(*n*^2) time, where *n* is the length of the input strings. The Four-Russians algorithm improves this complexity by a factor of (log *n*)^2 by using a lookup table. In this thesis, the algorithm is thoroughly examined and important implementation details are discussed, with special consideration given to parallelizing the algorithm and reducing the size of the lookup table. An implementation in C++ is provided and its performance is compared with a popular edit distance library in several experiments. The results indicate that the algorithm is a viable option in practice, but not optimal.

### What is in the implementation

The implementation is comprised of four files:

* four_russians.h -
The implementation of the Four-Russians algorithm.
There is no corresponding cpp file as templated functions have to be defined in the header (unless explicit instantiation is used).

* example.cpp -
Source code for a simple program that showcases how to use four_russians.h.
It takes two strings as command-line arguments and outputs the edit distance between them.

* test.cpp -
Source code for a simple program that tests the correctness of the algorithm.
It generates a decent amount of random strings and verifies that the edit distance computed is equal to the result from Edlib (a C++ edit distance library).

* benchmark.cpp -
Source code for a program that is capable of running all the benchmarks presented in the thesis.
It uses the Google Benchmark library to attempt to reduce the amount of boilerplate needed.

As you can see, each cpp file is meant to be compiled as a single program.

### How to build it

The two required libraries are:

* Google Benchmark -
https://github.com/google/benchmark -
version 1.3.0 or newer.

* Edlib -
https://github.com/Martinsos/edlib -
version 1.2.4 would be safest, but any version will probably work.

Any compiler supporting C++11 and OpenMP 3.0 (although OpenMP 3.0 is only needed to enable parallelization of loops with unsigned indices, otherwise OpenMP 2.0 will suffice) should be able to compile the implementation and the two required libraries. To be more precise, the implementation is compilable even without OpenMP, but will then lack parallelization (and emit a lot of warnings about unused parameters). However, to make sure all goes smoothly, you should use either GCC 4.8 (or newer) or Clang 3.4 (or newer).

A CMakeLists.txt file is included for your convenience. This is a configuration file for the most popular C++ build tool - CMake (https://cmake.org/). Edlib requires version 2.8 and Google Benchmark 3.5.1, so at least version 3.5.1 is required. The command

`sudo apt install cmake`

(or the equivalent in your package manager) should install a sufficiently new version.

The CMakeLists.txt is configured such so that it tries to find the libraries on your system, and if it does not find them, it expects their source code to be in the relative paths edlib/ and benchmark/ (or more precisely expects their CMakeLists.txt files there) (otherwise it fails). For these reasons, one possible way of building the implementation is the following sequence of commands (assuming you are on an Unix-like system, you have git installed and you are in the directory with the CMakeLists.txt file):

````
git clone https://github.com/google/benchmark.git
git clone https://github.com/Martinsos/edlib.git
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make
````

The three built binaries should then be located in the bin/ directory. The `-D CMAKE_BUILD_TYPE=Release` switch is important, as otherwise the implementation and the libraries will be build in debug mode. If needed, the equivalent of make clean for cmake is deleting the build folder.

### How to run it

#### example

Example usage:

`./example hello world`

Example output:

`The edit distance between "hello" and "world" is 4.`

#### test

Example usage:

`./test`

Example output (might take a few minutes):
```
Test case 1 is done.
Test case 2 is done.
Test case 3 is done.
...
````

#### benchmark

As the binary is build with the Google Benchmark library, it uses its command-line arguments. Most notably "--benchmark_filter=<regex>" runs only the benchmarks which names match the regex. By default the results are printed to the console, but "--benchmark_out=<filename>" can be used to print the output to <filename> in format according to "--benchmark_out_format={json|console|csv}".

Example usage:

`./benchmark --benchmark_filter="ecoli/small_m_and_large_n"`

Example output:
````
2019-05-14 13:36:01
Running ./benchmark
Run on (12 X 3400 MHz CPU s)
Load Average: 0.52, 0.58, 0.59
--------------------------------------------------------------------------------------------
Benchmark                                                  Time             CPU   Iterations
--------------------------------------------------------------------------------------------
ecoli/small_m_and_large_n/wf/4/repeats:5_mean           6157 us         6110 us            5
ecoli/small_m_and_large_n/wf/4/repeats:5_median         6172 us         6138 us            5
ecoli/small_m_and_large_n/wf/4/repeats:5_stddev         25.4 us         62.4 us            5
ecoli/small_m_and_large_n/FR/4/repeats:5_mean           4228 us         4217 us            5
ecoli/small_m_and_large_n/FR/4/repeats:5_median         4218 us         4236 us            5
ecoli/small_m_and_large_n/FR/4/repeats:5_stddev         18.4 us         42.1 us            5
ecoli/small_m_and_large_n/MBV/4/repeats:5_mean         11678 us        11670 us            5
ecoli/small_m_and_large_n/MBV/4/repeats:5_median       11678 us        11719 us            5
ecoli/small_m_and_large_n/MBV/4/repeats:5_stddev        61.9 us          204 us            5
...
````
The iteration count is slightly misleading, as each benchmark is run many times, averaged, and then all that is repeated five times and these aggregates are taken, so the real iteration count is for example 5 times 5000.

The filters which correspond to the benchmarks in the thesis are:
````
./benchmark --benchmark_filter="scaling_of_the_preprocessing_step"
./benchmark --benchmark_filter="ecoli/scaling_of_the_computation_step"
./benchmark --benchmark_filter="war_and_peace/scaling_of_the_computation_step"
./benchmark --benchmark_filter="ecoli/small_m_and_small_n"
./benchmark --benchmark_filter="war_and_peace/small_m_and_small_n"
./benchmark --benchmark_filter="ecoli/small_m_and_large_n"
./benchmark --benchmark_filter="war_and_peace/small_m_and_large_n"
./benchmark --benchmark_filter="ecoli/large_m_and_large_n"
./benchmark --benchmark_filter="war_and_peace/large_m_and_large_n"
````
