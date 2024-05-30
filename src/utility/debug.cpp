/**
 * @file debug.cpp
 *
 * @author A. M. Kazachkov
 * @date 2018-12-25
 */
#include "debug.hpp"

#include "gmic.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;
#include "utility.hpp"

#include <CoinPackedMatrix.hpp>
#include <CoinPackedVectorBase.hpp>
#include <CoinPackedVector.hpp>
#include <CoinShallowPackedVector.hpp>
#include <OsiSolverInterface.hpp>
#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>

/**
 * Print sparse vector
 *
 * @param vec           Vector to be printed 
 * @param use_newline   Whether to print the entire vector on one line or not
 */
void printVector(const CoinPackedVectorBase& vec, const bool use_newline) {
  int numElems = vec.getNumElements();
  const int* index = vec.getIndices();
  const double* element = vec.getElements();
  if (use_newline)
    fprintf(stdout, "Num elements is %d.", numElems);
  for (int i = 0; i < numElems; i++) {
    if (use_newline)
      fprintf(stdout, "\n\t");
    else if (i > 0)
      fprintf(stdout, ", ");
    fprintf(stdout, "(%d, %g)", index[i], element[i]);
  }
  fprintf(stdout, "\n");
} /* printVector (CoinPackedVectorBase) */

/** Print several vectors in CoinPackedVectorBase form */
void printVectors(const std::vector<CoinPackedVector>& vecs, const bool use_newline) {
  for (auto& vec : vecs)
    printVector(vec, use_newline);
} /* printVectors */

/**
 * Print dense array
 *
 * @param n     Number of elements in the vector
 * @param vec   Vector to be printed
 */
template <typename T>
void printVector(const int n, const T* vec, const bool use_newline, const bool sparse_mode) {
  for (int i = 0; i < n; ++i) {
    if (sparse_mode && std::abs(vec[i]) < 1e-7)
      continue;
    if (use_newline)
      fprintf(stdout, "\n");
    else if (i > 0)
      fprintf(stdout, ", ");

    if (std::abs(vec[i]-std::floor(vec[i])) != 0.0)
      fprintf(stdout, "\t(%d, %g)", i, static_cast<double>(vec[i]));
    else
      fprintf(stdout, "\t(%d, %d)", i, static_cast<int>(vec[i]));
  }
  fprintf(stdout, "\n");
} /* printVector (int, T*) */

/**
 * Print vector
 */
template <typename T>
void printVector(const std::vector<T>& vec, const bool use_newline, const bool sparse_mode) {
  printVector(vec.size(), vec.data(), use_newline, sparse_mode);
} /* printVector (std::vector<T>) */

/**
 * Print matrix (row-wise)
 *
 * @param mx_in   Matrix to be printed
 */
void printMatrix(const CoinPackedMatrix& mx_in) {
  const CoinPackedMatrix* mx; 
  if (mx_in.isColOrdered()) {
    CoinPackedMatrix* mx_rows = new CoinPackedMatrix;
    mx_rows->reverseOrderedCopyOf(mx_in);
    mx = mx_rows;
  } else {
    mx = &mx_in;
  }
  for (int i = 0; i < mx->getNumRows(); i++) {
    const CoinShallowPackedVector& vec = mx->getVector(i);
    fprintf(stdout, "Row %d: ", i);
    printVector(vec, false);
  }
  if (mx_in.isColOrdered()) {
    delete mx;
  }
} /* printMatrix */

/// @brief Test to make sure Gomory methods all produce same outcome
void testGomory(
    /// [in/out] non-const because we may need to enable factorization inside of \link generateGomoryCuts \endlink
    OsiSolverInterface* const solver,
    /// [in] parameters to use for GMIC generation
    const StrengtheningParameters::Parameters params) {
  OsiCuts currGMICs1, currGMICs3;
  generateGomoryCuts(currGMICs1, solver, 1, params.get(intParam::STRENGTHEN), params.get(doubleConst::AWAY), params.get(doubleConst::DIFFEPS), params.logfile);
  generateGomoryCuts(currGMICs3, solver, 3, params.get(intParam::STRENGTHEN), params.get(doubleConst::AWAY), params.get(doubleConst::DIFFEPS), params.logfile);
  const int num_gmics1 = currGMICs1.sizeCuts();
  const int num_gmics3 = currGMICs3.sizeCuts();
  if (num_gmics1 != num_gmics3) { printf("*** ERROR: %d gmics1 != %d gmics3\n", num_gmics1, num_gmics3); exit(1); }
  for (int i = 0; i < num_gmics1; i++) {
    OsiRowCut* cut = currGMICs1.rowCutPtr(i);
    const double rhs = cut->rhs();
    cut->mutableRow() /= rhs;
    cut->setLb(1.);
    cut->setUb(getInfinity());
  }
  for (int i = 0; i < num_gmics3; i++) {
    OsiRowCut* cut = currGMICs3.rowCutPtr(i);
    const double rhs = cut->rhs();
    cut->mutableRow() /= rhs;
    cut->setLb(1.);
    cut->setUb(getInfinity());
  }
  // Check difference
  for (int i = 0; i < num_gmics1; i++) {
    OsiRowCut cut1 = currGMICs1.rowCut(i);
    OsiRowCut cut2 = currGMICs3.rowCut(i);
    CoinPackedVector vec = cut1.row() - cut2.row();
    printf("Cut %d: sum = %f.\n", i, vec.sum());
  }
  printf("\n## DONE DEBUGGING GMICS ##\n");
} /* testGomory */

/// @brief Currently this is for debugging purposes for bm23
void checkCoefficientForColumn(
    /// [in] RCVMILP instance
    const OsiSolverInterface* const liftingSolver,
    /// [in] Solution to the RCVMILP instance
    const double* const solution,
    /// [in] Term to check
    const int term_ind,
    /// [in] Column to check
    const int col) {
  // Verify column col by taking dot product with solver matrix
  const CoinPackedMatrix* mat = liftingSolver->getMatrixByRow();
  const int num_rows = 20;
  const int num_cols = 27;
  const int mprime = num_rows + num_cols * 2;
  const int num_term_rows = 1; // assuming the same for all
  const int num_terms = 2;
  const int term_uvar_start_ind = 1 + mprime + term_ind * mprime;
  const int term_uvar_lb_start_ind = term_uvar_start_ind + num_rows;
  const int term_uvar_ub_start_ind = term_uvar_lb_start_ind + num_cols;
  const int term_u0var_start_ind = 1 + mprime + num_terms * mprime + term_ind * num_term_rows;
  double val_rows = 0;
  for (int row = 0; row < num_rows; row++) {
    const int rcvmilp_row_ind_for_col = term_ind * (num_cols + 1 + mprime) + col;
    const int rcvmilp_col_ind_for_row = term_uvar_start_ind + row;
    const double coeff = mat->getCoefficient(rcvmilp_row_ind_for_col, rcvmilp_col_ind_for_row);
    val_rows += coeff * solution[rcvmilp_col_ind_for_row];
  }
  val_rows /= -1. * solution[0]; // (recall that liftingSolver has \alpha \theta - v^t A^t = 0 as the constraint, so we negate the coefficients here to get back v^t A^t)
  const double val_lb = solution[term_uvar_lb_start_ind + col] / solution[0];
  const double val_ub = -1. * solution[term_uvar_ub_start_ind + col] / solution[0];
  
  double val_u0 = 0;
  for (int row = 0; row < num_term_rows; row++) {
    const int rcvmilp_row_ind_for_col = term_ind * (num_cols + 1 + mprime) + col;
    const int rcvmilp_col_ind_for_row = term_u0var_start_ind + row;
    const double coeff = mat->getCoefficient(rcvmilp_row_ind_for_col, rcvmilp_col_ind_for_row);
    val_u0 += coeff * solution[rcvmilp_col_ind_for_row];
  }
  
  printf("value for coefficient on cut for column %d = %f.\n", col, val_rows + val_lb + val_ub + val_u0);
} /* checkCoefficientForColumn */