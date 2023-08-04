/**
 * @file verify.cpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#include "verify.hpp"

// COIN-OR files
#include <CoinPackedMatrix.hpp>
#include <OsiSolverInterface.hpp>

// Project files
#include "utility.hpp"

/// @details Calculate \p alpha = (A^t) \p v (where A^t is the constraint matrix of \p term_solver and \p v is the Farkas multiplier).
/// 
/// In \p term_solver, A^t is the constraint matrix: first m+m_t rows of v correspond to A;D^t; the next n are bounds on the variables.
void getCutFromCertificate(
    /// [out] calculated cut coefficients
    std::vector<double>& alpha, 
    /// [in] Farkas multipliers
    const std::vector<double>& v,
    /// [in] LP solver corresponding to disjunctive term
    const OsiSolverInterface* const term_solver) {
  alpha.clear();
  alpha.resize(term_solver->getNumCols(), 0.0);

  const CoinPackedMatrix* mat = term_solver->getMatrixByCol();

  // std::vector<double> new_v(v.begin(), v.end());
  // for (int col = 0; col < term_solver->getNumCols(); col++) {
  //   double& val = new_v[term_solver->getNumRows() + col];
  //   val = std::abs(val);
  // }
  for (int col = 0; col < term_solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    alpha[col] += dotProduct(mat->getVectorSize(col), 
        mat->getIndices() + start, mat->getElements() + start, v.data());
    alpha[col] += v[term_solver->getNumRows() + col];
  }

  /*
  for (int col = 0; col < solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    for (int el_ind = 0; el_ind < mat->getVectorSize(col); el_ind++) {
      const int row = mat->getIndices()[start + el_ind];
      const double mult = (solver->getRowSense()[row] == 'L') ? -1. : 1.; //isNonBasicUBSlack(solver, row) ? -1.0 : 1.0;
      alpha[col] += mult * v[row];
    }
    const double mult = isNonBasicUBVar(solver, col) ? -1.0 : 1.0;
    alpha[col] += mult * v[solver->getNumRows() + col];
  }
  */
} /* getCutFromCertificate */

void checkCut(
    /// [out] number of errors
    int& num_errors,
    /// [out] total difference between cut coefficients
    double& total_diff,
    /// [in] cut coefficients
    const std::vector<double>& cut_coeff,
    /// [in] Farkas multipliers
    const std::vector<double>& v,
    /// [in] LP solver corresponding to disjunctive term
    const OsiSolverInterface* const term_solver
) {
  // Obtain the cut that the certificate yields (should be the same as the original cut)
  std::vector<double> new_coeff(term_solver->getNumCols());
  getCutFromCertificate(new_coeff, v, term_solver);

  num_errors = 0;
  total_diff = 0.;
  for (int i = 0; i < term_solver->getNumCols(); i++) {
    const double diff = cut_coeff[i] - new_coeff[i];
    if (greaterThanVal(std::abs(diff), 0.0)) {
      fprintf(stderr, "%d: cut: %g\tcalc: %g\tdiff: %g\n", i, cut_coeff[i], new_coeff[i], diff);
      num_errors++;
      total_diff += std::abs(diff);
    }
  }
} /* checkCut */