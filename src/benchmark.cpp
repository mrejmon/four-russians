#include <algorithm> // std::min
#include <cassert>
#include <cstddef> // size_t
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>
#include <edlib.h>

#include "four_russians.h"

// this file generates all the required benchmarks with the google benchmark library
// BM is used as a shorthand for benchmark

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------

// the results should be the same
// even if this number is changed
static unsigned seed = 18;

// an overexplained function which opens an input file
// and loads it whole into a vector
// returns success/fail bool
bool open_input_file(std::string filepath, std::vector<char> & buffer)
{
  // attempt to read the file
  std::ifstream file(filepath.c_str());
  // if we didn't manage to open it
  // return false
  if (!file.is_open())
  {
    return false;
  }
  // set the file stream to the end
  file.seekg(0, std::ios::end);
  // read the file size
  size_t file_size = file.tellg();
  // resize the buffer
  buffer.resize(file_size);  
  // set the file stream back to beginning
  file.seekg(0, std::ios::beg);
  // read the file
  file.read(buffer.data(), buffer.size());
  // if an unrecoverable error happens during reading
  // we return false
  if (file.bad())
  {
    return false;
  }
  // on the other hand, if the only thing that happened is that
  // we read fewer characters than we expected, try to continue execution
  size_t read_chars = file.gcount();
  if (read_chars != buffer.size())
  {
    std::cerr << "Warning: Did not manage to read the whole file " + filepath << std::endl;
    buffer.resize(read_chars);
  }
  
  return true;
}

// a quick implementation of the naive algorithm
// with column-major order of computation and linear space consumption
long wagner_fischer(const char * U, size_t U_size, const char * V, size_t V_size)
{
  std::vector<long> column(U_size + 1);
  for (size_t i = 0; i < U_size + 1; i++)
  {
    column[i] = i;
  }

  for (size_t j = 1; j < V_size + 1; j++)
  {
    long diagonal = j - 1;
    column[0] = j;
    for (size_t i = 1; i < U_size + 1; i++)
    {
      long tmp = column[i];
      column[i] = std::min(std::min(
        column[i] + 1,
        column[i - 1] + 1),
        diagonal + (U[i - 1] != V[j - 1]));
      diagonal = tmp;
    }
  }
  
  return column[U_size];
}

// -----------------------------------------------------------------------------
// ATOMIC BENCHMARKS
// -----------------------------------------------------------------------------

// runs the benchmark for the naive algorithm for the specified file and sizes of input strings
void BM_wagner_fischer(benchmark::State & state, const std::vector<char> & file, size_t U_size, size_t V_size)
{
  // this check here is to make sure the sizes are appropriate and the file was not tampered with
  if (file.size() < U_size || file.size() < V_size)
  {
    // benchmarking loop is not entered
    state.SkipWithError("Size of one of the input strings is bigger than the size of the input file.");
  }
  // to generate random locations we need to generate random numbers
  const char * U = file.data(), * V = file.data();
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> U_start(0, file.size() - U_size), V_start(0, file.size() - V_size);
  
  // this loop is repeated until its total run time is long enough according to the library
  // depending on the length of one iteration
  // the results are then printed to the console or to a file
  for (auto _ : state)
  {
    // move the ptrs to a random location
    U = file.data() + U_start(gen);
    V = file.data() + V_start(gen);
    // run the alg
    wagner_fischer(U, U_size, V, V_size);
  }
}

// the rest of these benchmarking functions is not commented as much
// since they are pretty much the same

// runs the benchmark for the preprocessing step for the specified tm, tn and number of threads
template <int tm, int tn>
void BM_four_russians_preprocessing_step(benchmark::State & state, int thread_cnt)
{
  FR::Four_russians<tm, tn> four_russians;
  for (auto _ : state)
  {
    four_russians.preprocessing_step(thread_cnt);
  }
}

// runs the benchmark for the computation step for the specified tm, tn, number of threads, file and sizes of input strings
template <int tm, int tn>
void BM_four_russians_computation_step(benchmark::State & state, const std::vector<char> & file, size_t U_size, size_t V_size, int thread_cnt, int chunk_i, int chunk_j)
{
  if (file.size() < U_size || file.size() < V_size)
  {
    state.SkipWithError("Size of one of the input strings is bigger than the size of the input file.");
  }
  const char * U = file.data(), * V = file.data();
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> U_start(0, file.size() - U_size), V_start(0, file.size() - V_size);
  
  FR::Four_russians<tm, tn> four_russians;
  four_russians.preprocessing_step(-1);
  
  for (auto _ : state)
  {
    U = file.data() + U_start(gen);
    V = file.data() + V_start(gen);
    four_russians.computation_step(U, U_size, V, V_size, thread_cnt, chunk_i, chunk_j);
  }
}

// runs the benchmark for edlib (Myers's bit-vector algorithm) for the specified file and sizes of input strings
void BM_edlib(benchmark::State & state, const std::vector<char> & file, size_t U_size, size_t V_size)
{
  if (file.size() < U_size || file.size() < V_size)
  {
    state.SkipWithError("Size of one of the input strings is bigger than the size of the input file.");
  }
  const char * U = file.data(), * V = file.data();
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> U_start(0, file.size() - U_size), V_start(0, file.size() - V_size);
  
  auto edlibConfig = edlibNewAlignConfig(U_size > V_size ? U_size : V_size, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
  for (auto _ : state)
  {
    U = file.data() + U_start(gen);
    V = file.data() + V_start(gen);
    edlibAlign(U, U_size, V, V_size, edlibConfig);
  }
}

// this benchmark tries to approximate the cost of accessing random positions in a file
// it is only used in the benchmark for small m and n, as otherwise it is very insignificant
void BM_overhead(benchmark::State & state, const std::vector<char> & file, size_t U_size, size_t V_size)
{
  if (file.size() < U_size || file.size() < V_size)
  {
    state.SkipWithError("Size of one of the input strings is bigger than the size of the input file.");
  }
  const char * U = file.data(), * V = file.data();
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> U_start(0, file.size() - U_size), V_start(0, file.size() - V_size);
  
  for (auto _ : state)
  {
    U = file.data() + U_start(gen);
    V = file.data() + V_start(gen);
    benchmark::DoNotOptimize(*U);
    benchmark::DoNotOptimize(*V);
  }
}

// -----------------------------------------------------------------------------
// COMBINED BENCHMARKS
// -----------------------------------------------------------------------------

// these following functions use the just defined functions
// and register them with the library for various needed sizes
// basically each function represents one graph in the thesis

void register_BMs_for_scalability_of_the_preprocessing_step()
{
  #define SCALING_OF_THE_PREPROCESSING_STEP(tm, tn) \
  benchmark::RegisterBenchmark((BM_name + "/scaling_of_the_preprocessing_step/"#tm"-"#tn"/" + s).data(), \
  &BM_four_russians_preprocessing_step<tm, tn>, i)-> \
  Repetitions(10)->ReportAggregatesOnly()->Unit(benchmark::kMicrosecond);
  // the parameters which will be used
  const int tm = 3, tn = 4;
  // loop over all threads
  for (int i = 1; i <= 12; i++)
  {
    auto s = std::to_string(i);
    // registers the BM
    benchmark::RegisterBenchmark(
      // name of the BM in the output
      (std::string("scaling_of_the_preprocessing_step/") + s).data(),
      // which BM are we actually running
      &BM_four_russians_preprocessing_step<tm, tn>,
      // its parameters
      i)
      // the benchmark is repeated 10 * (the determined number of iterations) times
      ->Repetitions(10)
      // only show mean, meadian and stddev
      ->ReportAggregatesOnly()
      // unit of measureament
      ->Unit(benchmark::kMicrosecond);
  }
}

// the rest of these functions is again similar,
// except that they are also parameterized by tm, tn and the file

void register_BMs_for_scalability_of_the_computation_step(std::string BM_name, std::vector<char> & file)
{
  size_t input_size = 65536;
  for (int i = 1; i <= 12; i++)
  {
    auto s = std::to_string(i);
    benchmark::RegisterBenchmark(
      (BM_name + "/scaling_of_the_computation_step/1-1/" + s).data(),
      &BM_four_russians_computation_step<1, 1>,
      file, input_size, input_size, i, 64, 64)
      ->Repetitions(i)
      ->ReportAggregatesOnly()
      ->Unit(benchmark::kMillisecond);
      
    benchmark::RegisterBenchmark(
      (BM_name + "/scaling_of_the_computation_step/2-3/" + s).data(),
      &BM_four_russians_computation_step<2, 3>,
      file, input_size, input_size, i, 64, 64)
      ->Repetitions(i)
      ->ReportAggregatesOnly()
      ->Unit(benchmark::kMillisecond);
      
    benchmark::RegisterBenchmark(
      (BM_name + "/scaling_of_the_computation_step/2-4/" + s).data(),
      &BM_four_russians_computation_step<2, 4>,
      file, input_size, input_size, i, 64, 64)
      ->Repetitions(i)
      ->ReportAggregatesOnly()
      ->Unit(benchmark::kMillisecond);
      
    benchmark::RegisterBenchmark(
      (BM_name + "/scaling_of_the_computation_step/3-4/" + s).data(),
      &BM_four_russians_computation_step<3, 4>,
      file, input_size, input_size, i, 64, 64)
      ->Repetitions(i)
      ->ReportAggregatesOnly()
      ->Unit(benchmark::kMillisecond);
  }
}

template <int tm, int tn>
void register_BMs_for_small_m_small_n(std::string BM_name, std::vector<char> & file)
{
  for (int i = 10; i <= 50; i++)
  {
    auto s = std::to_string(i);
    
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_small_n/overhead/" + s).data(),
    &BM_overhead, file, i, i)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kNanosecond);
 
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_small_n/wf/" + s).data(),
    &BM_wagner_fischer,
    file, i, i)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kNanosecond);
 
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_small_n/FR/" + s).data(),
    &BM_four_russians_computation_step<tm, tn>,
    file, i, i, 1, 64, 64)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kNanosecond);
    
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_small_n/MBV/" + s).data(),
    &BM_edlib, file, i, i)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kNanosecond);
  }  
}

template <int tm, int tn>
void register_BMs_for_small_m_large_n(std::string BM_name, std::vector<char> & file)
{
  for (int i = 4; i <= 64; i += 4)
  {
    auto s = std::to_string(i);
    size_t input_size = 1048576;
   
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_large_n/wf/" + s).data(),
    &BM_wagner_fischer,
    file, i, input_size)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
 
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_large_n/FR/" + s).data(),
    &BM_four_russians_computation_step<tm, tn>,
    file, i, input_size, 1, 64, 64)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
    
    benchmark::RegisterBenchmark(
    (BM_name + "/small_m_and_large_n/MBV/" + s).data(),
    &BM_edlib,
    file, i, input_size)
    ->Repetitions(5)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
  }   
}

template <int single_tm, int single_tn, int multi_tm, int multi_tn>
void register_BMs_for_large_m_large_n(std::string BM_name, std::vector<char> & file)
{
  for (int i = 64; i <= 32768; i *= 2)
  {
    auto s = std::to_string(i);
    int reps = i <= 32768 ? 5 : 1;
    
    benchmark::RegisterBenchmark(
    (BM_name + "/large_m_and_large_n/wf/" + s).data(),
    &BM_wagner_fischer,
    file, i, i)
    ->Repetitions(reps)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
 
    benchmark::RegisterBenchmark(
    (BM_name + "/large_m_and_large_n/FR - single-threaded/" + s).data(),
    &BM_four_russians_computation_step<single_tm, single_tn>,
    file, i, i, 1, 64, 64)
    ->Repetitions(reps)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
    
    reps = i <= 32768 ? 10 : 5;
 
    benchmark::RegisterBenchmark(
    (BM_name + "/large_m_and_large_n/FR - multi-threaded/" + s).data(),
    &BM_four_russians_computation_step<multi_tm, multi_tn>,
    file, i, i, -1, 64, 64)
    ->Repetitions(reps)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
    
    benchmark::RegisterBenchmark(
    (BM_name + "/large_m_and_large_n/MBV/" + s).data(),
    &BM_edlib,
    file, i, i)
    ->Repetitions(reps)
    ->ReportAggregatesOnly()
    ->Unit(benchmark::kMicrosecond);
  }
}

// -----------------------------------------------------------------------------
// MAIN
// -----------------------------------------------------------------------------

int main(int argc, char** argv) {
  // filepaths are hardcoded for simplicity
  std::string ecoli_filepath = "../input/ecoli.txt";
  std::string war_and_peace_filepath = "../input/war_and_peace.txt";
  
  // prefixes for the names of benchmarks
  std::string ecoli_BM_name = "ecoli";
  std::string war_and_peace_BM_name = "war_and_peace";
  
  // containers for the files
  std::vector<char> ecoli_buffer;
  std::vector<char> war_and_peace_buffer;
  
  // try to open and read the files
  bool success = open_input_file(ecoli_filepath, ecoli_buffer);
  if (!success)
  {
    std::perror(("Error while reading file " + ecoli_filepath).c_str());
    return -1;
  }
  success = open_input_file(war_and_peace_filepath, war_and_peace_buffer);
  if (!success)
  {
    std::perror(("Error while reading file " + war_and_peace_filepath).c_str());
    return -1;
  }  
  
  // register all the benchmarks
  register_BMs_for_scalability_of_the_preprocessing_step();
  
  register_BMs_for_scalability_of_the_computation_step(ecoli_BM_name, ecoli_buffer);
  register_BMs_for_scalability_of_the_computation_step(war_and_peace_BM_name, war_and_peace_buffer);
  
  register_BMs_for_small_m_small_n<4, 2>(ecoli_BM_name, ecoli_buffer);
  register_BMs_for_small_m_small_n<5, 2>(war_and_peace_BM_name, war_and_peace_buffer);
  
  register_BMs_for_small_m_large_n<4, 3>(ecoli_BM_name, ecoli_buffer);
  register_BMs_for_small_m_large_n<4, 4>(war_and_peace_BM_name, war_and_peace_buffer);
  
  register_BMs_for_large_m_large_n<3, 4, 2, 3>(ecoli_BM_name, ecoli_buffer);
  register_BMs_for_large_m_large_n<3, 4, 4, 3>(war_and_peace_BM_name, war_and_peace_buffer);
  
  // library handles the rest
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
