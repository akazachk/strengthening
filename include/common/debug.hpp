/**
 * @file debug.hpp
 * @brief Useful functions for debugging that are not necessary in the rest of the code
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#include <string>
#include <vector>

namespace StrengtheningParameters {
  struct Parameters;
}
class OsiSolverInterface;
class OsiCuts;
class CoinPackedVectorBase;
class CoinPackedVector;
class CoinPackedMatrix;

void printVector(const CoinPackedVectorBase& vec, const bool use_newline = true);
void printVectors(const std::vector<CoinPackedVector>& vecs, const bool use_newline = true);
template <typename T>
void printVector(const int n, const T* vec, const bool use_newline = true, const bool sparse_mode = true);
template <typename T>
void printVector(const std::vector<T>& vec, const bool use_newline = true, const bool sparse_mode = true);

// For gdb
inline void printVectorInt(const int n, const int* vec, const bool use_newline = true, const bool sparse_mode = true) {
  printVector(n, vec, use_newline, sparse_mode);
}
inline void printVectorDouble(const int n, const double* vec, const bool use_newline = true, const bool sparse_mode = true) {
  printVector(n, vec, use_newline, sparse_mode);
}
inline void printVectorInt(const std::vector<int>& vec, const bool use_newline = true, const bool sparse_mode = true) {
  printVectorInt(vec.size(), vec.data(), use_newline, sparse_mode);
}
inline void printVectorDouble(const std::vector<double>& vec, const bool use_newline = true, const bool sparse_mode = true) {
  printVectorDouble(vec.size(), vec.data(), use_newline, sparse_mode);
}

void printMatrix(const CoinPackedMatrix& mx_in);

/// @brief Test to make sure Gomory methods all produce same outcome
void testGomory(OsiSolverInterface* const solver, const StrengtheningParameters::Parameters params);

/// @brief Currently this is for debugging purposes for bm23
void checkCoefficientForColumn(
    const OsiSolverInterface* const liftingSolver,
    const double* const solution,
    const int term_ind,
    const int col);