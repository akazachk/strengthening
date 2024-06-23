/**
 * @file strengthening_debug.cpp
 *
 * @author A. M. Kazachkov
 * @date 2024-06-22
 */
#include "strengthening_debug.hpp"

#include "gmic.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;

#include <CoinPackedVector.hpp>
#include <OsiSolverInterface.hpp>
#include <OsiCuts.hpp>

/// @brief Test to make sure Gomory methods all produce same outcome
void testGomory(
    /// [in/out] non-const because we may need to enable factorization inside of \link generateGomoryCuts \endlink
    OsiSolverInterface* const solver,
    /// [in] parameters to use for GMIC generation
    const StrengtheningParameters::Parameters params) {
  OsiCuts currGMICs1, currGMICs3;
  generateGomoryCuts(currGMICs1, solver, 1, params.get(intParam::STRENGTHEN), params.get(intConst::MIN_SUPPORT_THRESHOLD), params.get(doubleParam::MAX_SUPPORT_REL), params.get(doubleParam::AWAY), params.get(doubleConst::DIFFEPS), params.logfile);
  generateGomoryCuts(currGMICs3, solver, 3, params.get(intParam::STRENGTHEN), params.get(intConst::MIN_SUPPORT_THRESHOLD), params.get(doubleParam::MAX_SUPPORT_REL), params.get(doubleParam::AWAY), params.get(doubleConst::DIFFEPS), params.logfile);
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

