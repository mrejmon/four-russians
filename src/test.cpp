#include <climits>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <edlib.h>

#include "four_russians.h"

// this program tests the implementation with a sequence of basic blackbox tests
// unit tests were unfortunately scrapped as the inner implementation changed quite a bit
// as time went on

// this function generates a random string of a certain length over a certain alphabet
std::string generate_random_string(unsigned seed, const std::string & alphabet, size_t length)
{
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(0, alphabet.size() - 1);
  std::string output(length, '0');
  for (size_t i = 0; i < length; ++i)
  {
    output[i] = alphabet[dis(gen)];
  }
  return output;
}

// this function checks that the result of Four-Russians is correct for certain pairs of lengths and counts of threads
template <int tm, int tn>
void test_blackbox(unsigned seed, const std::string & alphabet,
  const std::vector<std::pair<size_t, size_t>> & length_pairs, const std::vector<int> & thread_cnts)
{
  // prepare Four-Russians
  FR::Four_russians<tm, tn> four_russians;
  four_russians.preprocessing_step(-1);
  for (auto length_pair : length_pairs)
  {
    // generate input strings and compute the result with Edlib
    std::string U = generate_random_string(seed++, alphabet, length_pair.first);
    std::string V = generate_random_string(seed++, alphabet, length_pair.second);
    EdlibAlignResult edlib_result = edlibAlign(U.data(), U.size(), V.data(), V.size(),
      edlibDefaultAlignConfig());
    for (auto thread_cnt : thread_cnts)
    {
      // compute the result with the implementation
      long fr_result = four_russians.computation_step(U.data(), U.size(), V.data(), V.size(), thread_cnt);
      // and check they are the same
      if (edlib_result.editDistance != fr_result)
      {
        std::cout << "TEST FAILED!\n"
        << edlib_result.editDistance << " != " << fr_result << '\n'
        << "Parameters: tm = " << tm << "| tn = " << tn << "| U length = " << length_pair.first
        << "| V length = " << length_pair.second << "| thread_cnt = " << thread_cnt << '\n'
        << "U = " << U << '\n'
        << "V = " << V << '\n' 
        << std::endl;
      }
    }
  }
}

int main()
{
  // tests should of course succeed even if this number is changed
  unsigned seed = 18;
   
  // prepare a small, normal and maximal alphabet
  std::string alphabet_small = "01";
  std::string alphabet_normal = "abcdefghijklmnopqrstuvwxyz";
  std::string alphabet_max;
  for (char c = CHAR_MIN; c < CHAR_MAX - 1; c++)
  {
    alphabet_max.push_back(c);
  }
  // prepare a lot of pairs of small lengths
  std::vector<std::pair<size_t, size_t>> small_lengths;
  size_t start = 1, end = 150;
  for (size_t U_length = start; U_length < end; ++U_length)
  {
    for (size_t V_length = start; V_length < end; ++V_length)
    {
      small_lengths.push_back(std::make_pair(U_length, V_length));
    }
  }
    
  // prepare few pairs of big lengths
  std::vector<std::pair<size_t, size_t>> big_lengths;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> dis(5000, 10000);
  size_t amount = 100;
  for (size_t i = 0; i < amount; ++i)
  {
    big_lengths.push_back(std::make_pair(dis(gen), dis(gen)));
  }
  
  // tests will be run with 1, 2 or many threads
  std::vector<int> thread_cnts{ 1, 2, -1 };
     
  // test cases for tm = tn = 1
  test_blackbox<1, 1>(dis(gen), alphabet_normal, small_lengths, thread_cnts);
  std::cout << "Test case 1 done." << std::endl;
  
  test_blackbox<1, 1>(dis(gen), alphabet_normal, big_lengths, thread_cnts);
  std::cout << "Test case 2 done." << std::endl;

  // test cases for tn > tm
  test_blackbox<2, 5>(dis(gen), alphabet_normal, small_lengths, thread_cnts);
  std::cout << "Test case 3 done." << std::endl;
  
  test_blackbox<2, 5>(dis(gen), alphabet_normal, big_lengths, thread_cnts);
  std::cout << "Test case 4 done." << std::endl;

  // test cases for tm = tn
  test_blackbox<5, 2>(dis(gen), alphabet_normal, small_lengths, thread_cnts);
  std::cout << "Test case 5 done." << std::endl;
  
  test_blackbox<5, 2>(dis(gen), alphabet_normal, big_lengths, thread_cnts);
  std::cout << "Test case 6 done." << std::endl;

  // test cases for a small alphabet
  test_blackbox<3, 3>(dis(gen), alphabet_small, small_lengths, thread_cnts);
  std::cout << "Test case 7 done." << std::endl;
    
  test_blackbox<3, 3>(dis(gen), alphabet_small, big_lengths, thread_cnts);
  std::cout << "Test case 8 done." << std::endl;

  // test cases for a maximal alphabet
  test_blackbox<3, 4>(dis(gen), alphabet_max, small_lengths, thread_cnts);
  std::cout << "Test case 9 done." << std::endl;
  
  test_blackbox<3, 4>(dis(gen), alphabet_max, big_lengths, thread_cnts);
  std::cout << "Test case 10 done." << std::endl;
}
