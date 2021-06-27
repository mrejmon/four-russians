#pragma once

#include <array>
#include <cstddef> // size_t
#include <cstdint> // uint8_t etc.
#include <type_traits> // std::conditional_t
#include <vector>

// ifdef enables conditional compilation without OpenMP.
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace FR
{
  // These functions do not necessarily make sense as a part of the Four_Russians class,
  // but to not expose them to the user they are put into an anonymous namespace.
  namespace
  {
    // These need to be forward declared here as they are used in the class definition.
    constexpr int floor_log2(int x);
    constexpr int ceil_log2(int x);
    // This templated struct is the smallest unsigned integer type
    // capable of holding "bits" bits.
    template <int bits>
    struct uintxx_t;
  } // end of anonymous namespace

  // The algorithm is implemented as a single class "Four_Russians"
  // to hide the implementation of the LUT from the user,
  // with two non-type template parameters -- integers tm and tn.
  template <int tm, int tn>
  class Four_russians
  {
  // The entire public interface is three functions.
  public:
    // Reports the size of the LUT table in bytes in compile time according to tm and tn.
    constexpr uint64_t get_LUT_size();
    
    // Allocates the table and computes all values in it.
    // thread_cnt determines the number of threads used by OpenMP (default is 1).
    // If thread_cnt is nonpositive, maximum threads is used.
    void preprocessing_step(int thread_cnt = 1);
    
    // Computes the edit distance between U_size characters of U and V_size characters of V.
    // thread_cnt determines the number of threads used by OpenMP (default is 1).
    // If thread_cnt is nonpositive, maximum threads is used.
    // chunk_i and chunk_j are the height and width (in number of blocks) of the chunks.
    // If chunk_i or chunk_j is nonpositive, their value is set to the default (which is 64).
    // Returns -1 if the preprocessing step has not yet been run.
    long computation_step(const char * U, size_t U_size, const char * V, size_t V_size,
      int thread_cnt = 1, int chunk_i = 64, int chunk_j = 64) const;

  private:
    // Checking template parameters.
    // Values <1 don't make sense.
    // Maximum is 20 to prevent overflow when calculating size of the LUT
    // (practical maximum is lower anyway).
    static_assert(tm >= 1 && tm <= 20,
      "ERROR: Four_russians - Parameter \"tm\" must be a number between 1 and 20 (inclusive).");
    static_assert(tn >= 1 && tn <= 20,
      "ERROR: Four_russians - Parameter \"tn\" must be a number between 1 and 20 (inclusive).");

    // These functions calculate the number of combinations/values of vectors Y, X, A or B.
    // They are mostly used for calculating the size of the LUT.
    static constexpr uint64_t y_vals();
    static constexpr uint64_t x_vals();
    template <int t>
    static constexpr uint64_t ab_vals();

    // LUT_size = number of combinations of vectors Y, X, A, and B multiplied between itself.
    constexpr static uint64_t LUT_size = y_vals() * x_vals() * ab_vals<tm>() * ab_vals<tn>();

    // Vectors A and B have their own types to optimize memory usage as much as possible.
    using A = typename uintxx_t<ceil_log2(ab_vals<tm>())>::type;
    using B = typename uintxx_t<ceil_log2(ab_vals<tn>())>::type;
    
    // Simple struct representing one element of the LUT.
    struct LUT_elem
    {
      A a;
      B b;
    };

    // The LUT itself.
    std::vector<LUT_elem> LUT;
    
    // This bool is used in the computation step to check whether
    // the preprocessing step was run.
    bool preprocessing_complete = false;

    // These functions take care of all the encoding and decoding of
    // vectors Y, X, A, and B from bits to chars/integers and vice versa.
    template <int t>
    static void decode_ab_vals(size_t & input, std::array<int, t + 1> & output);
    template <int t>
    static void decode_ab_steps(size_t input, std::array<int, t> & output);
    static void decode_x(size_t & input, std::array<char, tm> & output);
    static void decode_y(size_t & input, std::array<char, tn> & output);
    template <int t>
    static size_t encode_ab(std::array<int, t + 1> & input);
    static size_t encode_y(const char * input, std::array<char, 256> & alphabet);
    static void encode_x(const char * input, const std::array<char, 256> & alphabet, size_t & output);

    // Private functions for the preprocessing_step.
    // They are described more in the definition.
    static LUT_elem compute_block(std::array<char, tn> & Y, std::array<char, tm> & X,
      std::array<int, tm + 1> & A, std::array<int, tn + 1> & B, size_t YXAB);

    // Private functions for the first step of the computation step -- computing matrix.
    // They are described more in the definition.
    void compute_matrix(std::vector<A> & S, std::vector<B> & T,
      const char * U, const char * V, int thread_cnt, int chunk_i, int chunk_j) const;
    void compute_chunk(std::vector<A> & S, std::vector<B> & T,
      const char * U, const char * V, size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;

    // Private functions for the second step of the computation step - computing edit distance.
    // They are described more in the definition.
    static long compute_edit_distance(std::vector<A> & S, std::vector<B> & T,
      const char * U, size_t U_size, const char * V, size_t V_size);
    template <int divisible_t, int not_divisible_t, typename T>
    static long compute_missing(const T & divisible_vec, const char * divisible_str,
      const char * not_divisible_str, size_t not_divisible_str_size, std::vector<long> & output);
    static void wagner_fischer_step(long prev_corner, long next_corner, char corner_char,
      const char * chars, std::vector<long> & values);
    static long wagner_fischer(const char * U, size_t U_size, const char * V, size_t V_size);
  }; // end of class Four_Russians

  namespace
  {
    // Compile time computation of floor(log2(x)) and ceil(log2(x)).
    // Used for computing the number of bits required by A and B.
    constexpr int floor_log2(int x)
    {
      return x == 1 ? 0 : 1 + floor_log2(x >> 1);
    }

    constexpr int ceil_log2(int x)
    {
      return x == 1 ? 0 : 1 + floor_log2(x - 1);
    }
    
    // The correct type alias is determined by the number of bits.
    // This is slightly more generic than needed, as if A or B need more than 16 bits,
    // the index to the LUT is very probably already too big.
    template <int bits>
    struct uintxx_t
    {
      using type =
        // The syntax of std::conditional is
        // std::conditional<condition, true, false>
        // and ::type is equal to either true or false.
        typename std::conditional<(bits <= 8), uint8_t,
        typename std::conditional<(bits <= 16), uint16_t,
        typename std::conditional<(bits <= 32), uint32_t,
        uint64_t>::type>::type>::type;
    };
    
    // There also a few more functions outside of the class.
    
    // A small helper function that returns the minimum out of three arguments.
    // Used in the preprocessing step when computing blocks 
    // and also in the end when computing missing rows and columns.
    template <typename T> T min3(T a, T b, T c)
    {
      return b < c ? (a < b ? a : b) : (a < c ? a : c);
    }
    
    // Compile time power and factorial. Used for computing number of values
    // vectors Y, X, A, and B can have.
    constexpr uint64_t pow(int num, int exp)
    {
      return exp == 0 ? 1 : num * pow(num, exp - 1);
    }

    constexpr uint64_t factorial(int n)
    {
      return n == 0 ? 1 : n * factorial(n - 1);
    }
  } // end of anonymous namespace
    
  // -------------------------------------------------------------------
  // FUNCTIONS RELATED TO COMPUTING LUT SIZE
  // -------------------------------------------------------------------
  
  // Public function.
  // Reports the size of the LUT table in bytes in compile time according to tm and tn.
  template <int tm, int tn>
  constexpr uint64_t Four_russians<tm, tn>::get_LUT_size()
  {
    return LUT_size * sizeof(LUT_elem);
  }

  // Returns the number of values/combinations of vector Y.
  template <int tm, int tn>
  constexpr uint64_t Four_russians<tm, tn>::y_vals()
  {
    return factorial(tn);
  }

  // Returns the number of values/combinations of vector X.
  template <int tm, int tn>
  constexpr uint64_t Four_russians<tm, tn>::x_vals()
  {
    return pow(tn + 1, tm);
  }

  // Returns the number of combinations of vector A or B.
  // Template parameter should be either tm for A or tn for B.
  // Also used to compute the default value in vectors S and T.
  template <int tm, int tn>
  template <int t> constexpr uint64_t Four_russians<tm, tn>::ab_vals()
  {
    return pow(3, t);
  }
    
  // -------------------------------------------------------------------
  // FUNCTIONS RELATED TO DECODING AND ENCODING
  // -------------------------------------------------------------------

  // Removes t steps from input and stores them in output as values, skipping first elem.
  // Used in the preprocessing step for decoding index. The input is taken as a reference and modified,
  // so that it can be sent immediately to the next decoding function.
  // For example:
  // input = -1/+1/0 (encoded)
  // output = [3, x, x, x] (x means not important)
  // >>>
  // output = [3, 2, 3, 3]
  // The values are encoded in reverse order, than presented here (that is 0/+1/-1).
  // Since we are removing steps from the end, we need the first element to be in the last place.
  template <int tm, int tn>
  template <int t>
  void Four_russians<tm, tn>::decode_ab_vals(size_t & input, std::array<int, t + 1> & output)
  {
    for (int j = 1; j < t + 1; ++j)
    {
      // - 1 is here because we encode -1 as 0, 0 as 1 and 1 as 2,
      // since modulo 3 returns 0, 1 and 2.
      output[j] = (input % 3) - 1 + output[j - 1];
      input /= 3;
    }   
  }
  
  // Decodes t steps from an integer representing either vector A (t=tm) or B(t=tn),
  // and stores them in output as steps. Used in the computation step (second part)
  // for adding up all the values in a row (or column).
  template <int tm, int tn>
  template <int t>
  void Four_russians<tm, tn>::decode_ab_steps(size_t input, std::array<int, t> & output)
  {
    for (int i = 0; i < t; ++i)
    {
      int val = input % 3;
      output[i] = val - 1;
      input /= 3;
    }
  }

  // Removes tm chars from input and stores them in output. Used in the preprocessing step for decoding index.
  // The input is taken as a reference and modified, so that it can be sent immediately to the next decoding function.
  // The array is filled backward, for consistency with Y.
  template <int tm, int tn>
  void Four_russians<tm, tn>::decode_x(size_t & input, std::array<char, tm> & output)
  {
    for (int i = tm; i-- > 0; )
    {
      output[i] = input % (tn + 1);
      input /= (tn + 1);
    }
  }

  // Removes tn chars from input and stores them in output. Used in the preprocessing step for decoding index.
  // The input is taken as a reference and modified, so that it can be sent immediately to the next decoding function.
  // The array is filled backward, since encoding Y from the front is much easier than from the back
  // (that is the first element in the first place, last in the last place, as we remove elements with modulo from the back).
  template <int tm, int tn>
  void Four_russians<tm, tn>::decode_y(size_t & input, std::array<char, tn> & output)
  {
    for (int i = tn; i > 1; --i)
    {
      // + 1 is here since we encode 1 as 0, 2 as 1, 3 as 2 etc.
      // For the reason see the encode_y() function.
      output[i - 1] = (input % i) + 1;
      input /= i;
    }
    output[0] = 1;
  }
  
  // Encodes t elements of input as an integer representing vector A (t=tm) or B(t=tn).
  // Encodes the elements as steps, not as values.
  // Used in the preprocessing step to encode A and B before storing them in the LUT.
  template <int tm, int tn>
  template <int t>
  size_t Four_russians<tm, tn>::encode_ab(std::array<int, t + 1> & input)
  {
    size_t result = 0;
    for (int i = t + 1; i-- > 1; )
    {
      result *= 3;
      // + 1 is needed since we encode -1 as 0, 0 as 1 and 1 as 2.
      result += (input[i] - input[i - 1]) + 1;
    }
    return result;
  }

  // Encodes tn characters of input as an integer representing vector Y.
  // For example: bcacb -> 01210. Alphabet is indexed by the characters and in it their encoded
  // values are saved. The alphabet is the used in encode_x().
  template <int tm, int tn>
  size_t Four_russians<tm, tn>::encode_y(const char * input, std::array<char, 256> & alphabet)
  {
    size_t result = 0;
    // We know that the first elem is always 1 (and also that tn >=1).
    // The static_cast might be needed here, since char might be signed
    // (unless they are possibly converted to size_t by the [] operator of std::array).
    alphabet[static_cast<unsigned char>(input[0])] = 1;
    // next_val represents the largest yet unassigned value.
    unsigned char next_val = 2;

    for (int i = 1; i < tn; ++i)
    {
      // "Shift" result to the left.
      result *= (i + 1);
      // Read the next character.
      unsigned char val = input[i];

      // If we didn't see it yet, assign it the largest yet unassigned value,
      // otherwise just use the value we already assigned to the character.
      if (alphabet[val] == 0)
      {
        alphabet[val] = next_val++;
      }
      // 0 is saved for meaning that the character was not in Y,
      // but we are only saving the characters which are in Y,
      // so we encode 1 as 0, 2 as 1, 3 as 2 etc. We could also use a different
      // value instead of 0 as sentinel value, for example 100, as we know there
      // will never be more than 20 unique characters (since max(tn) == 20).
      result += alphabet[val] - 1;
    }
    return result;
  }

  // Encodes tm characters of input as an integer representing substrings X using
  // the alphabet prepared in encode_y().
  template <int tm, int tn>
  void Four_russians<tm, tn>::encode_x(const char * input, const std::array<char, 256> & alphabet, size_t & output)
  {
    for (int i = 0; i < tm; ++i)
    {
      output *= tn + 1;
      // static_cast is needed here, since char might be signed.
      output += alphabet[static_cast<unsigned char>(input[i])];
    }
  }
  
  // -------------------------------------------------------------------
  // FUNCTIONS RELATED TO THE PREPROCESSING STEP
  // -------------------------------------------------------------------

  // Public function.
  // Allocates the table and computes all values in it.
  // thread_cnt determines the number of threads used by OpenMP (default is 1).
  // If thread_cnt is nonpositive, maximum threads is used.
  template <int tm, int tn>
  void Four_russians<tm, tn>::preprocessing_step(int thread_cnt)
  {
    // Preprocessing_step() has to be run only once.
    if (preprocessing_complete)
    {
      return;
    }
    
  // ifdef for conditional compilation without OpenMP.
  #ifdef _OPENMP
    // Nonpositive thread_cnt means use max threads.
    if (thread_cnt < 1)
    {
      thread_cnt = omp_get_max_threads();
    }
  #else
    thread_cnt = 1;
  #endif

    // Allocate the lookup table.
    LUT.resize(LUT_size);

    // Construct some arrays here so they don't have to
    // be reconstructed in each iteration.
    std::array<char, tn> Y;
    std::array<char, tm> X;
    // Vector A and B are one larger, because it makes it easier to compute
    // the block with the naive algorithm. See compute_block().
    std::array<int, tm + 1> A;
    std::array<int, tn + 1> B;

    // The main loop of the preprocessing step.
    // Compute block for each value of the index and save the results in the LUT.
    // Simple parallelization with omp
    #pragma omp parallel for private(Y, X, A, B) num_threads(thread_cnt)
    for (size_t YXAB = 0; YXAB < LUT_size; ++YXAB)
    {
      LUT[YXAB] = compute_block(Y, X, A, B, YXAB);
    }
    
    // Marking it as complete enables the use of computation_step().
    preprocessing_complete = true;
  }

  // Computes a block according to the values in the index (YXAB) and
  // stores the last column and the last row (A' and B') in the LUT (encoded).
  // All the values that are in the vectors Y, X, A, and B when the function
  // is called are ignored and overwritten.
  template <int tm, int tn>
  typename Four_russians<tm, tn>::LUT_elem Four_russians<tm, tn>::compute_block(
    std::array<char, tn> & Y, std::array<char, tm> & X,
    std::array<int, tm + 1> & A, std::array<int, tn + 1> & B,
    size_t YXAB)
  {
    // Decode the index value. The order of vectors in the index is Y X A B.
    // But the decoding is done in reverse order as it is much faster with modulo encoding.
    // The first elem of A and B is set to 0, as it is required by the decode_ab() function.
    A[0] = 0;
    B[0] = 0;
    decode_ab_vals<tn>(YXAB, B);
    decode_ab_vals<tm>(YXAB, A);
    decode_x(YXAB, X);
    decode_y(YXAB, Y);

    // Now we have to compute the block with the Wagner-Fischer algorithm.
    // The order is row-major and only linear space is consumed, that is vector B.
    // We also need to update vector A so that in the end it is equal to the last column.
    for (int i = 1; i < tm + 1; ++i)
    {
      // Get the new diagonal value
      int diagonal = A[i - 1];
      // Set the first elem in the row
      B[0] = A[i];
      // Update the saved diagonal value to the rightmost elem B
      A[i - 1] = B[tn];
      for (int j = 1; j < tn + 1; ++j)
      {
        // Save the value of the cell right above
        int tmp = B[j];
        // Update the cell in the current row
        B[j] = min3(
          B[j] + 1,
          B[j-1] + 1,
          diagonal + (X[i - 1] != Y[j - 1]));
        // The diagonal value is now the value of the cell, which was previously right above
        diagonal = tmp;
      }
    }
    // Update the last cell in vector A
    A[tm] = B[tn];

    // Encode the vectors A' and B' into two integers.
    LUT_elem result;
    result.a = encode_ab<tm>(A);
    result.b = encode_ab<tn>(B);
    return result;
  }

  // -------------------------------------------------------------------
  // COMPUTATION STEP - COMPUTING MATRIX
  // -------------------------------------------------------------------

  // Public function.
  // Compute the edit distance between U_size characters of U and V_size characters of V.
  // thread_cnt determines the number of threads used by OpenMP (default is 1).
  // If thread_cnt is nonpositive, maximum threads is used.
  // chunk_i and chunk_j are the height and width (in number of blocks) of the chunks.
  // If chunk_i or chunk_j is nonpositive, their value is set to the default (which is 64).
  // Returns -1 if the preprocessing step has not yet been run.
  template <int tm, int tn>
  long Four_russians<tm, tn>::computation_step(const char * U, size_t U_size,
    const char * V, size_t V_size, int thread_cnt, int chunk_i, int chunk_j) const
  {
    // If preprocessing_step() hasn't been run, always return -1.
    if (!preprocessing_complete)
    {
      return -1;
    }
    
    // Then check for easily solvable cases that would break the implementation.
    if (U_size == 0)
    {
      return V_size;
    }
    if (V_size == 0)
    {
      return U_size;
    }
    if (U_size < tm || V_size < tn)
    {
      return wagner_fischer(U, U_size, V, V_size);
    }

    // Prepare the vectors S and T with the required sizes and default values.
    // The default values are all ones, for example tm=3, then default values are 1/1/1 (encoded).
    // As the values are packed in perfectly, we just have to subtract 1 from the total number of values. 
    std::vector<A> S(U_size / tm, ab_vals<tm>() - 1);
    std::vector<B> T(V_size / tn, ab_vals<tn>() - 1);
    
    // Compute the matrix.
    // After this S is the last column of the matrix and T the last row.
    compute_matrix(S, T, U, V, thread_cnt, chunk_i, chunk_j);
    
    // Finish off with computing the edit distance from vectors S and T.
    // If the blocks did not fit into the matrix perfectly,
    // compute also the surplus rows and columns.
    return compute_edit_distance(S, T, U, U_size, V, V_size);
  }

  // The parallel wrapper around the actual computation of the matrix.
  // The matrix is split into chunks of size chunk_i x chunk_j, and chunks
  // on the diagonal are computed concurrently. This is done in practice
  // mainly by keeping the vectors in one piece and just computing the required
  // indices/offsets.
  template <int tm, int tn>
  void Four_russians<tm, tn>::compute_matrix(std::vector<A> & S, std::vector<B> & T,
    const char * U, const char * V, int thread_cnt, int chunk_i, int chunk_j) const
  {
  // ifdef for conditional compilation without OpenMP.
  #ifdef _OPENMP
    // Nonpositive thread_cnt means use max threads.
    if (thread_cnt < 1)
    {
      thread_cnt = omp_get_max_threads();
    }
  #else
    thread_cnt = 1;
  #endif
    
    // If nonpositive chunk sizes were entered, just use the default value.
    if (chunk_i < 1)
    {
      chunk_i = 64;
    }
    if (chunk_j < 1)
    {
      chunk_j = 64;
    }
    
    // If we only use one thread, overwrite the chunk parameters,
    // as they are probably not useful with one thread.
    if (thread_cnt == 1)
    {
      chunk_i = S.size();
      chunk_j = T.size();
    }
    
    // Calculate the number of chunks in one column.
    size_t chunks_column = S.size() / chunk_i + ((S.size() % chunk_i) != 0);
    // Calculate the number of chunks in one row.
    size_t chunks_row = T.size() / chunk_j + ((T.size() % chunk_j) != 0);
    // We start with a diagonal of size 1.
    size_t diag_size = 1;
    
    /*
    For simplicity, the algorithm is split into two for loops.
    The first one calculates the upper-left half of the matrix
    (without the "middle" diagonal), that is for example all the
    chunks here represented with a number:
    
    0,0  1,1  2,2  x,x  x,x
    1,0  2,1  x,x  x,x  x,x
    2,0  x,x  x,x  x,x  x,x
    x,x  x,x  x,x  x,x  x,x
    
    The first number is diag_idx, the second is chunk_idx.
    */
    // Loop over diagonals
    for (size_t diag_idx = 0; diag_idx < chunks_column - 1; ++diag_idx)
    {
      // Loop over chunks in a diagonal, simply parallelized with OpenMP.
      #pragma omp parallel for num_threads(thread_cnt)
      for (size_t chunk_idx = 0; chunk_idx < diag_size; ++chunk_idx)
      {
        // Computation of these loop indices is tricky to explain with words,
        // drawing on a piece of paper works better.
        size_t i_start = diag_idx * chunk_i - chunk_idx * chunk_i;
        size_t i_end = i_start + chunk_i;
        size_t j_start = chunk_idx * chunk_j;
        size_t j_end = j_start + chunk_j > T.size() ? T.size() : j_start + chunk_j;
        // Computes a chunk from (i_start, j_start) to (i_end, j_end).
        compute_chunk(S, T, U, V, i_start, i_end, j_start, j_end);
      }
      // This check makes non-square matrices work (with square matrices we
      // can always just add one).
      if (diag_idx + 1 < chunks_row)
        diag_size += 1;
    }
    /* 
    The second loop computes the lower-right half of the matrix:
    
    x,x  x,x  x,x  0,3  1,3
    x,x  x,x  0,2  1,2  2,2
    x,x  0,1  1,1  2,1  3,1
    0,0  1,0  2,0  3,0  4,0
    
    The diag_idx starts again from 0 and the loop is quite similar to the first one.
    */
    for (size_t diag_idx = 0; diag_idx < chunks_row; ++diag_idx)
    {
      #pragma omp parallel for num_threads(thread_cnt)
      for (size_t chunk_idx = 0; chunk_idx < diag_size; ++chunk_idx)
      {
        size_t i_start = (chunks_column - 1) * chunk_i - chunk_idx * chunk_i;
        size_t i_end = i_start + chunk_i > S.size() ? S.size() : i_start + chunk_i;
        size_t j_start = diag_idx * chunk_j + chunk_idx * chunk_j;
        size_t j_end = j_start + chunk_j > T.size() ? T.size() : j_start + chunk_j;
        compute_chunk(S, T, U, V, i_start, i_end, j_start, j_end);
      }
      if (chunks_row - diag_idx <= chunks_column)
        diag_size -= 1;
    }
  }

  // Computes the matrix, that is vectors S and T, from (i_start, j_start) to (i_end, j_end).
  template <int tm, int tn>
  void Four_russians<tm, tn>::compute_chunk(std::vector<A> & S, std::vector<B> & T,
    const char * U, const char * V, size_t i_start, size_t i_end, size_t j_start, size_t j_end) const
  {
    // First we move the ptrs to the offsets.
    U += i_start * tm;
    V += j_start * tn;
    // Create an alphabet for encode_y() and encode_x()
    // and fill it with zeroes.
    std::array<char, 256> alphabet;
    alphabet.fill(0);
    // Loop over each column.
    for (size_t j = j_start; j < j_end; ++j)
    {
      // We only encode Y once each column.
      size_t Y = encode_y(V, alphabet);
      // Loop over each row.
      for (size_t i = i_start; i < i_end; ++i)
      {
        // We start building the index.
        size_t YABX = Y;
        // encode_x() adds X to the index.
        encode_x(U, alphabet, YABX);
        // The ptr has to be moved manually, to be ready for the next row.
        U += tm;
        // "Shift" the index to the left.
        YABX *= ab_vals<tm>();
        // Add A from position [i, j-1], which is found in S[i].
        YABX += S[i];
        // "Shift" the index to the left again.
        YABX *= ab_vals<tn>();
        // Add B from position [i-1, j], which is found in T[j].
        YABX += T[j];
        // Finally, do the lookup.
        LUT_elem result = LUT[YABX];
        // Assign the results.
        S[i] = result.a;
        T[j] = result.b;
      }
      // At the end of each column, reset the alphabet for characters in Y.
      for (size_t i = 0; i < tn; ++i)
      {
        // static_cast is needed because char might be signed.
        alphabet[static_cast<unsigned char>(V[i])] = 0;
      }
      // Reset U so that it points at the first element again.
      U -= (i_end - i_start) * tm;
      // And move V to the next column.
      V += tn;
    }
  }
  
  // -------------------------------------------------------------------
  // COMPUTATION STEP - COMPUTING EDIT DISTANCE
  // -------------------------------------------------------------------

  // Either computes the resulting edit distance by adding up all steps in T,
  // or, if the blocks did not fit into the matrix perfectly,
  // compute also the missing rows and columns with S, U and V.
  // This is actually slightly more annoying than it seems it should be,
  // and not very interesting either.
  template <int tm, int tn>
  long Four_russians<tm, tn>::compute_edit_distance(std::vector<A> & S, std::vector<B> & T,
    const char * U, size_t U_size, const char * V, size_t V_size)
  {
    // If size of U is divisible by tm and size of V is divisible by tn.
    if (S.size() * tm == U_size && T.size() * tn == V_size)
    {
      // Bottom left corner of the matrix is equal to U_size,
      // so we skip the addition of the first row by setting result to U_size.
      long result = U_size;
      // Prepare an array for decoding.
      std::array<int, tn> block;
      // Loop over the blocks in T.
      for (size_t i = 0; i < T.size(); ++i)
      {
        // Decode the block into individual steps.
        decode_ab_steps<tn>(T[i], block);
        // Loop over the steps.
        for (int j = 0; j < tn; ++j)
        {
          // Add them up to the result.
          result += block[j];
        }
      }
      return result;
    }
    // If only size of U is divisible by tm.
    else if (S.size() * tm == U_size)
    {
      std::vector<long> missing;
      // This functions computes the missing part of the last row (in values, not in steps).
      compute_missing<tm, tn>(S, U, V, V_size, missing);
      // So we just return the last value.
      return missing.back();
    }
    // If only size of V is divisible by tn.
    else if (T.size() * tn == V_size)
    {
      std::vector<long> missing;
      // This functions computes the missing part of the last column (in values, not in steps).
      compute_missing<tn, tm>(T, V, U, U_size, missing);
      // So we just return the last value.
      return missing.back();
    }
    // If neither string is divisible.
    else
    {
      // Compute both the missing part of S and the missing part of T.
      std::vector<long> hori, vert;
      compute_missing<tm, tn>(S, U, V, V_size, hori);
      size_t corner = compute_missing<tn, tm>(T, V, U, U_size, vert);
      // This leaves us a with small rectangle of uncomputed values,
      // which we compute with the Wagner-Fischer algorithm.
      // Start by moving the ptrs to strings to the correct position.
      U += S.size() * tm;
      V += T.size() * tn;
      // We iterate row by row.
      // old_val represents the previous value in the first column
      // and new_val the new one.
      size_t old_val = corner, new_val;
      for (size_t i = 0; i < vert.size(); ++i)
      {
        new_val = vert[i];
        // This functions computes one whole row.
        wagner_fischer_step(old_val, new_val, U[i], V, hori);
        old_val = vert[i];
      }
      return hori.back();
    }
  }

  // Computes a missing part of the matrix.
  // There are two possibilities:
  // If divisible_t = tm, not_divisible_t = tn, T = std::vector<A>,
  // divisible_vec = S, divisible_str = U, not_divisible_str = V
  // the the missing part of the last row is computed.
  // Or if divisible_t = tn, not_divisible_t = tm, T = std::vector<B>,
  // divisible_vec = T, divisible_str = V, not_divisible_str = U
  // the the missing part of the last column is computed.
  // The missing part is then stored in output, and the value just before that
  // returned as the return value (helpful when the last rectangle needs to be computed).
  template <int tm, int tn>
  template <int divisible_t, int not_divisible_t, typename T>
  long Four_russians<tm, tn>::compute_missing(const T & divisible_vec, const char * divisible_str,
    const char * not_divisible_str, size_t not_divisible_str_size, std::vector<long> & output)
  {
    // Find first index which was not computed during Four-Russians.
    size_t missing_size = not_divisible_str_size % not_divisible_t;
    size_t first_index = not_divisible_str_size - missing_size;
    
    // Move the string ptr and prepare the output vector for computation.
    not_divisible_str += first_index;
    output.resize(missing_size);
    for (size_t i = 1; i <= missing_size; ++i)
    {
      output[i - 1] = first_index + i;
    }
    
    // Loop over blocks of divisible_vec and decode them into steps.
    std::array<int, divisible_t> block;
    long old_val = first_index, new_val = 0;
    for (size_t i = 0; i < divisible_vec.size(); ++i)
    {
      decode_ab_steps<divisible_t>(divisible_vec[i], block);
      // Then loop over each step and compute each iteration
      // with the Wagner-Fischer alg.
      for (int j = 0; j < divisible_t; ++j)
      {
        new_val = old_val + block[j];
        wagner_fischer_step(old_val, new_val, divisible_str[i * divisible_t + j], not_divisible_str, output);
        old_val = new_val;
      }
    }
    return new_val;
  }

  // Computes one iteration of the Wagner-Fischer algorithm.
  // Since the size of the values vector is the size of the row/column minus one,
  // the min3 call has to be duplicated. This means this function and the ones that
  // call it could probably be refactored to be more readable.
  /*
  Explanation of the parameters/variables:
  Before:
                              chars[0]  chars[1]  ... chars[size-1]
  
                 prev_corner  values[0] values[1] ... values[size-1]
  corner_char    next_corner
  
  After:
                              chars[0]  chars[1]  ... chars[size-1]
  
                 prev_corner
  corner_char    next_corner  values[0] values[1] ... values[size-1]
  */
  template <int tm, int tn>
  void Four_russians<tm, tn>::wagner_fischer_step(long prev_corner, long next_corner,
    char corner_char, const char * chars, std::vector<long> & values)
  {
    long tmp1, tmp2;
    tmp1 = min3(
      prev_corner + (corner_char != chars[0]),
      next_corner + 1,
      values[0] + 1);
    for (size_t i = 1; i < values.size(); ++i)
    {
      tmp2 = tmp1;
      tmp1 = min3(
        tmp1 + 1,
        values[i] + 1,
        values[i - 1] + (corner_char != chars[i]));
      values[i - 1] = tmp2;
    }
    values.back() = tmp1;
  }

  // Computes the edit distance between U_size characters of U and V_size characters of V.
  // Fallback if m or n is smaller than tm or tn.
  template <int tm, int tn>
  long Four_russians<tm, tn>::wagner_fischer(const char * U, size_t U_size, const char * V, size_t V_size)
  {
    std::vector<long> column(U_size);
    for (size_t i = 0; i < U_size; ++i)
    {
      column[i] = i + 1;
    }
    for (size_t i = 0; i < V_size; ++i)
    {
      wagner_fischer_step(i, i + 1, V[i], U, column);
    }
    return column.back();
  }
} // end of FR namespace
