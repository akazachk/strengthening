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
#include "Disjunction.hpp"
#include "utility.hpp"

/// @details Return the set of rows and cols that have a nonzero multiplier in the certificate for some disjunctive term.
void getRowsColsFromCertificate(
    /// [out] set of rows
    std::set<int>& rows,
    /// [out] set of cols
    std::set<int>& cols,
    /// [in] Farkas multipliers
    const TermCutCertificate& v,
    /// [in] LP solver corresponding to disjunctive term
    const OsiSolverInterface* const term_solver) {
  // TODO
  rows.clear();
  cols.clear();

  const CoinPackedMatrix* mat = term_solver->getMatrixByCol();

  for (int col = 0; col < term_solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    for (int el_ind = 0; el_ind < mat->getVectorSize(col); el_ind++) {
      const int row = mat->getIndices()[start + el_ind];
      if (greaterThanVal(std::abs(v[row]), 0.0)) {
        rows.insert(row);
        cols.insert(col);
      }
    }
    if (greaterThanVal(std::abs(v[term_solver->getNumRows() + col]), 0.0)) {
      cols.insert(col);
    }
  }
} /* getRowsColsFromCertificate */

/// @details Calculate \p alpha = (A^t) \p v (where A^t is the constraint matrix of \p term_solver and \p v is the Farkas multiplier).
/// 
/// In \p term_solver, A^t is the constraint matrix: first m+m_t rows of v correspond to A;D^t; the next n are bounds on the variables.
void getCutFromCertificate(
    /// [out] calculated cut coefficients
    std::vector<double>& alpha, 
    /// [in] Farkas multipliers
    const TermCutCertificate& v,
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
} /* getCutFromCertificate (term_solver) */

/// @details Calculate \p alpha = (A^t) \p v (where A^t is the constraint matrix of \p solver combined with the globally-valid inequalities from \p disj, and \p v is the Farkas multiplier).
void getCutFromCertificate(
    /// [out] calculated cut coefficients
    std::vector<double>& alpha,
    /// [in] Farkas multipliers
    const TermCutCertificate& v,
    /// [in] Disjunction from which cut was derived
    const Disjunction* const disj,
    /// [in] original MILP instance
    const OsiSolverInterface* const solver) {
  alpha.clear();
  alpha.resize(solver->getNumCols(), 0.0);

  const CoinPackedMatrix* mat = solver->getMatrixByCol();

  for (int col = 0; col < solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    alpha[col] += dotProduct(mat->getVectorSize(col), 
        mat->getIndices() + start, mat->getElements() + start, v.data());
    for (int extra_row_ind = 0; extra_row_ind < disj->common_changed_var.size(); extra_row_ind++) {
      const int var = disj->common_changed_var[extra_row_ind];
      if (var != col) continue;
      const int row = solver->getNumRows() + extra_row_ind;
      const double coeff = disj->common_changed_bound[extra_row_ind] <= 0 ? 1. : -1.;
      alpha[col] += coeff * v[row];
    }
    alpha[col] += v[solver->getNumRows() + col]; // for lb or ub multiplier
  }
} /* getCutFromCertificate (disj, solver) */

/// @details If \p disj is given, then assume that \p solver is original MILP instance;
/// otherwise \p solver is the LP solver corresponding to disjunctive term.
void checkCut(
    /// [out] number of errors
    int& num_errors,
    /// [out] total difference between cut coefficients
    double& total_diff,
    /// [in] cut coefficients
    const std::vector<double>& cut_coeff,
    /// [in] Farkas multipliers
    const TermCutCertificate& v,
    /// [in] LP solver corresponding to disjunctive term OR the original MILP instance
    const OsiSolverInterface* const solver,
    /// [in] Disjunction from which cut was derived (if \p solver is the original MILP instance)
    const Disjunction* const disj) {
  // Obtain the cut that the certificate yields (should be the same as the original cut)
  std::vector<double> new_coeff(solver->getNumCols());
  if (disj == NULL) {
    getCutFromCertificate(new_coeff, v, solver);
  } else {
    getCutFromCertificate(new_coeff, v, disj, solver);
  }

  num_errors = 0;
  total_diff = 0.;
  for (int i = 0; i < solver->getNumCols(); i++) {
    const double diff = cut_coeff[i] - new_coeff[i];
    if (greaterThanVal(std::abs(diff), 0.0)) {
      fprintf(stderr, "%d: cut: %g\tcalc: %g\tdiff: %g\n", i, cut_coeff[i], new_coeff[i], diff);
      num_errors++;
      total_diff += std::abs(diff);
    }
  }
} /* checkCut */

/// @details Run #checkCut and print num errors
///
/// This is a wrapper for #checkCut that creates num_errors and total_diff variables and prints the number of errors.
int checkCutHelper(
    /// [in] cut coefficients
    const std::vector<double>& cut_coeff,
    /// [in] Farkas multipliers
    const TermCutCertificate& v,
    /// [in] LP solver corresponding to disjunctive term OR the original MILP instance
    const OsiSolverInterface* const solver,
    /// [in] Disjunction from which cut was derived (if \p term_solver is the original MILP instance)
    const Disjunction* const disj,
    /// [in] Log file
    FILE* const logfile) {
  int num_errors = 0;
  double total_diff = 0;
  checkCut(num_errors, total_diff, cut_coeff, v, solver, disj);
  
  if (num_errors > 0) {
    const bool should_continue = isZero(total_diff, 1e-3);
    if (should_continue) {
      // Send warning
      warning_msg(warnstring,
          "Number of differences between true and calculated cuts: %d. Total difference: %g. Small enough difference that we will try to continue, but beware of numerical issues.\n",
          num_errors, total_diff);
    } else {
      // Exit
      error_msg(errorstring,
          "Number of differences between true and calculated cuts: %d. Total difference: %g. Exiting.\n",
          num_errors, total_diff);
      writeErrorToLog(errorstring, logfile);
    }
// #ifdef USE_EIGEN
//     fprintf(stderr, "x:\n");
//     for (int i = 0; i < solver->getNumCols(); i++) {
//       fprintf(stderr, "x[%d] = %g\n", i, x(i));
//     }
//     fprintf(stderr, "b:\n");
//     for (int i = 0; i < solver->getNumCols(); i++) {
//       fprintf(stderr, "b[%d] = %g\n", i, b(i));
//     }
// #endif // USE_EIGEN
    fprintf(stderr, "v:\n");
    for (int i = 0; i < (int) v.size(); i++) {
      fprintf(stderr, "v[%d] = %g\n", i, v[i]);
    }
    if (!should_continue) {
      exit(1);
    }
  } // num_errors > 0

  return num_errors;
} /* checkCutHelper */