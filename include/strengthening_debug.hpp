/**
 * @file strengthening_debug.hpp
 * @brief Useful functions for debugging strengthening code that are not necessary in the rest of the code
 *
 * @author A. M. Kazachkov
 * @date 2024-06-22
 */
#include <string>
#include <vector>

namespace StrengtheningParameters {
  struct Parameters;
}
class OsiSolverInterface;

/// @brief Test to make sure Gomory methods all produce same outcome
void testGomory(OsiSolverInterface* const solver, const StrengtheningParameters::Parameters params);

/// @brief Currently this is for debugging purposes for bm23
void checkCoefficientForColumn(
    const OsiSolverInterface* const liftingSolver,
    const double* const solution,
    const int term_ind,
    const int col);
