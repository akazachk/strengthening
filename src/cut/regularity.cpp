/**
 * @file regularity.cpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#include "regularity.hpp"

// COIN-OR files
#include <CbcModel.hpp>
#include <CoinPackedMatrix.hpp>
#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>
#include <OsiSolverInterface.hpp>

// Project files
// #include "analysis.hpp" // cutCertificate
#include "Disjunction.hpp"
#include "Parameters.hpp"
#include "SolverHelper.hpp" // setLPSolverParameters, setIPSolverParameters
#include "SolverInterface.hpp"
#include "utility.hpp" // isInfinity, error_msg
#include "verify.hpp" // checkCutHelper

#ifdef USE_EIGEN
#include "eigen.hpp"
#endif

#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#include <gurobi_c++.h>
#endif

/// @details Number of rows should be solver->getNumRows() + num_common_rows + number of finite bounds
int calculateNumRowsAtilde(
    /// [in] Disjunction from which to get globally-valid inequalities
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver) {
  // Let mprime be the number of rows of Atilde
  // This is the number of original constraints + number of globally-valid bound changes
  // We also add the number of lower and upper bounds to mprime
  int mprime_bounds = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
      const double lb = solver->getColLower()[col];
      const double ub = solver->getColUpper()[col];

      if (!isInfinity(std::abs(lb))) {
          mprime_bounds++;
      }
      if (!isInfinity(std::abs(ub))) {
          mprime_bounds++;
      }
  }
  const int num_common_rows = disj->common_changed_var.size() + disj->common_ineqs.size();
  return solver->getNumRows() + num_common_rows + mprime_bounds;
} /* calculateNumRowsAtilde */

/// @details Matrix \p Atilde contains the original constraints, globally-valid inequalities, and variable bounds.
/// It will have dimensions as specified in #calculateNumRowsAtilde
/// (# rows + # globally-valid inequalities + # finite bounds)
void prepareAtilde(
    /// [out] Coefficient matrix with original constraints, globally-valid inequalities, and variable bounds
    CoinPackedMatrix& Atilde,
    /// [out] Right-hand side vector with original constraints, globally-valid inequalities, and variable bounds
    std::vector<double>& btilde,
    /// [in] Disjunction from which to get globally-valid inequalities
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver,
    /// [in] Log where to output error messages
    FILE* logfile) {
  Atilde.clear();
  btilde.clear();

  // Put Ax \ge b into Atilde and btilde
  // (It may be the case that solver's rows are not all in >= form!
  // We currently ignore that and handle it by setting the nonpositivity of the corresponding variables.)
  Atilde.copyOf(*(solver->getMatrixByCol())); // copy A into mx, col-ordered, dimesions = [num_rows, num_cols]
  Atilde.reverseOrdering(); // make it row-ordered
  btilde.assign(solver->getRightHandSide(), solver->getRightHandSide() + solver->getNumRows()); // copy b into btilde

  // Add rows corresponding to globally-valid bound changes to Atilde
  // All will be stored in >= form
  const double one = 1.;
  const double negone = -1.;
  for (int i = 0; i < (int) disj->common_changed_var.size(); i++) {
      const int col = disj->common_changed_var[i];
      const double coeff = (disj->common_changed_bound[i] <= 0) ? one : negone;
      const double val = disj->common_changed_value[i];

      Atilde.appendRow(1, &col, &coeff);
      btilde.push_back(coeff * val); // lower bound row
  }

  // Add rows corresponding to globally-valid common_ineqs to Atilde
  // We assume these are all >= constraints
  for (int i = 0; i < (int) disj->common_ineqs.size(); i++) {
      const OsiRowCut* currCut = &(disj->common_ineqs[i]);
      Atilde.appendRow(currCut->row().getNumElements(), currCut->row().getIndices(), currCut->row().getElements());
      btilde.push_back(currCut->rhs());
  }

  // Add rows corresponding to variable lower bounds to Atilde
  // These will be stored as x_j \ge lb_j
  for (int col = 0; col < solver->getNumCols(); col++) {
      const double val = solver->getColLower()[col];
      if (!isInfinity(std::abs(val))) {
          Atilde.appendRow(1, &col, &one);
          btilde.push_back(val);
      }
  }

  // Add rows corresponding to variable upper bounds to Atilde
  // These will be stored as -x_j \ge -ub_j
  for (int col = 0; col < solver->getNumCols(); col++) {
      const double val = solver->getColUpper()[col];
      if (!isInfinity(std::abs(val))) {
          Atilde.appendRow(1, &col, &negone);
          btilde.push_back(negone * val);
      }
  }

  // Verify number of rows of Atilde
  const int mtilde = calculateNumRowsAtilde(disj, solver);
  if (Atilde.getNumRows() != mtilde) {
      error_msg(errorstring,
          "Number of rows of Atilde (%d) does not match expected (%d).\n",
          Atilde.getNumRows(), mtilde);
      writeErrorToLog(errorstring, logfile);
      exit(1);
  }
} /* prepareAtilde */

/**
 * @details Generates the cut-generating linear program in the extended formulation space, given a fixed cut.
 * 
 * The feasible region of the disjunctive program is as follows, over variables x \in \R^n:
 *  \tilde{A} x \ge \tilde{b}
 *  \vee_{t \in T} (D^t x \ge D^t_0).
 * 
 * Let Atilde x \ge btilde encompass all constraints, including bounds 
 * *and* globally-valid bound changes at the root from the disjunction.
 * In practice, for bounds, we will first append the lower bounds, then the upper bounds.
 * 
 * Let m' denote the number of rows of Atilde.
 * Let m_t denote the number of rows of D^t.
 * 
 * The cut \alpha^T x \ge \beta is given in \p cut.
 * 
 * The CGLP has nonnegative variables (u^t,u_0^t) for the certificate for each term t of the disjunction.
 * In addition, we add binary variables \delta_j, j \in [n], to determine which columns of Atilde are used to certify the cut.
 * Finally, we need a nonnegative variable theta to determine if the cut needs to be scaled.
 * 
 * One extra parameter we add is \lambda, which will be used to state that each multiplier u^t_i is at most \lambda * \delta_i (a binary variable).
 * The value \lambda = 1 is default, but if the entries in Atilde are large then it is better to use a larger \lambda to avoid numerical issues.
 * 
 *  min_{(u^t,u_0^t)_{t \in T}, \delta, \theta} 
 *  (0)        -\theta
 *  (0')  \sum_{i \in [m']} \delta_i         
 *  s.t. 
 *  (1)  \alpha \theta               - Atilde^T u^t     - (D^t)^T   u^t_0     = 0,                 \forall t \in T,
 *  (2) -\beta  \theta               + btilde^T u^t     + (D^t_0)^T u^t_0   \ge 0,                 \forall t \in T,
 *  (3)                     \delta_i          - u^t_i                       \ge 0,                 \forall t \in T, \forall i \in [m'],
 *  (4)   \sum_{i \in [m']} \delta_i                                        \le n,
 *  (4')  \sum_{i \in [m']} \delta_i                                        \ge n,
 *  (5)   \sum_{i \in M'}   \delta_i                                        \le rank(Atilde_M'),   \forall M' \subseteq [m'],
 *  (6)                                         u^t,                u^t_0   \ge 0,                 \forall t in T,
 *  (7)                     \delta_i                                        \in {0,1},             \forall i \in [m'],
 *  (8)         \theta                                                      \in [0,1].
 * 
 * Num rows:
 *   (n+1) * |T|  // for the first two constraints
 *   + m' * |T|   // for the third constraint
 *   + 1          // for the fourth constraint
 *   + 2^{m'}.    // for the fifth constraint
 * Num variables:
 *   \theta = 1
 *   \delta = m'
 *      u^t = m'
 *    u^t_0 = m_t
 * = 1 + m' + m' * |T| + \sum_{t \in T} m_t.
 */
void genRCVMILPFromCut(
    /// [out] Cut-generating linear program to be generated
    OsiSolverInterface* const liftingSolver,
    /// [in] Cut \alpha^T x \ge \beta for which we are seeking a certificate
    const OsiRowCut* const cut,
    /// [in] Disjunction for which the cut is being generated
    const Disjunction* const disj,
    /// [in] Original solver for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Parameters (for logfile and deciding where to print the MILP)
    const StrengtheningParameters::Parameters& params,
    /// [in] Whether to use original RCVMILP or one with min sum delta objective and sum delta >= n constraint
    const bool use_min_sum_delta) {
  const int num_terms = disj->terms.size();
  const int num_rows = solver->getNumRows();
  // const int num_common_rows = disj->common_changed_var.size() + disj->common_ineqs.size();
  const int num_cols = solver->getNumCols();
  // const double* rhs = solver->getRightHandSide();
  const double one = 1.;
  const double negone = -1.;
  const int mprime = calculateNumRowsAtilde(disj, solver);

  // Let total_mt be the number of rows of all D^t
  int total_mt = 0;
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      total_mt += disj->terms[term_ind].changed_var.size() + disj->terms[term_ind].ineqs.size();
  }

  // Prepare rows corresponding to the cut
  CoinPackedVector alpha_beta; // will be [alpha; beta; 0]
  alpha_beta.append(cut->row());
  alpha_beta.insert(num_cols, -1. * cut->rhs());
  
  // Create block [ [\alpha; -\beta; 0], [0; 0; I_{m'}] ] repeated for all terms
  CoinPackedMatrix theta_delta; // col-ordered
  theta_delta.appendCol(alpha_beta);
  for (int i = 0; i < mprime; i++) {
      CoinPackedVector delta_i; // will be [0 (n rows); e_i]
      delta_i.insert(num_cols + 1 + i, 1.);
      theta_delta.appendCol(delta_i);
  }
  
  // Prepare block Atilde and btilde
  CoinPackedMatrix Atilde;
  std::vector<double> btilde;
  prepareAtilde(Atilde, btilde, disj, solver, params.logfile);

  // Create common_mx = [-Atilde^T; btilde^T; -I_{m'}]
  // This is repeated for all terms
  CoinPackedMatrix common_mx; // col-ordered
  common_mx.setDimensions(0, mprime);
  common_mx.reverseOrdering(); // make it row-ordered

  // Append -Atilde to common_mx
  // We want Atilde.getVector(i) to return column i of Atilde, which is multiplied by -1. and inputted as row i of common_mx
  Atilde.reverseOrdering(); // make it col-ordered
  if (!Atilde.isColOrdered()) {
      throw std::logic_error("genRCVMILPFromCut:: Atilde is not col-ordered!");
  }
  for (int col = 0; col < Atilde.getNumCols(); col++) {
      common_mx.appendRow(-1. * Atilde.getVector(col));
  }

  // Append btilde to common_mx
  std::vector<int> indices(num_rows);
  for (int i = 0; i < num_rows; i++) {
      indices[i] = i;
  }
  common_mx.appendRow(num_rows, &indices[0], &btilde[0]);

  // Append -I_{m'} to common_mx, to set whether corresponding indicator variables will be used
  // for constraints \delta_i \ge u^t_i
  // NB: if the constraint is <=, then u^t_i <= 0,
  // so we instead need to add \delta_i \ge -u^t_i
  const char* row_sense = solver->getRowSense();
  for (int i = 0; i < mprime; i++) {
    if (i < num_rows && row_sense[i] == 'L') {
      common_mx.appendRow(1, &i, &one);
    } else {
      common_mx.appendRow(1, &i, &negone);
    }
  }

  // Verify size of common_mx, and if incorrect, throw an error giving the sizes
  if (common_mx.getNumRows() != num_cols + 1 + mprime) {
      error_msg(errorstring,
          "genRCVMILPFromCut:: Num rows (mx, actual, predicted) = (Atilde, %d, %d); (common_mx, %d, %d).\n",
          Atilde.getNumRows(), mprime, common_mx.getNumRows(), num_cols + 1 + mprime);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
  }

  // Create the CGLP constraint matrix
  // Constraints: |T| * ( (1) + (2) + (3) ) + 1
  //  = |T| * ( num_cols + 1 + mprime ) + 1
  // Num variables:
  //  \theta = 1
  //  \delta = m'
  //     u^t = m' for each term
  //   u^t_0 = m_t for each term
  //  = 1 + m' + m' * |T| + \sum_{t \in T} m_t.
  const int num_cglp_constraints = num_terms * (num_cols + 1 + mprime) + 1;
  const int num_cglp_vars = 1 + mprime + num_terms * mprime + total_mt;
  CoinPackedMatrix mx; // create a new col-ordered matrix
  mx.reserve(num_cglp_constraints, num_cglp_vars);
  mx.setDimensions(0., num_cglp_vars);

  // For term t, we will prepend t-1 copies of the all-zeroes matrix with dimensions same as common_mx
  CoinPackedMatrix zero_mx; // col-ordered
  zero_mx.setDimensions(common_mx.getNumRows(), common_mx.getNumCols());
  
  // For each term, create the constraints consisting of
  // theta_delta, the appropriate number of copies of zero_mx, and a copy of common_mx to mx,
  // then add the term-specific constraints
  for (int term_ind = 0; term_ind < num_terms; term_ind++) {
    const DisjunctiveTerm& term = disj->terms[term_ind];
    CoinPackedMatrix term_mx; // col-ordered
    term_mx.setDimensions(common_mx.getNumRows(), 0.);

    // Add theta_delta, which is col-ordered and has dimensions [n + 1 + m', 1 + m']
    // Variables theta and delta with indices 0, 1, ..., m'
    term_mx.rightAppendPackedMatrix(theta_delta);

    // Add zero_mx for the previous terms
    // Variables u^{t'} for t' < term_ind
    // Indices are m'+1 to m'+term_ind*m'
    for (int t = 0; t < term_ind; t++) {
      term_mx.rightAppendPackedMatrix(zero_mx);
    }

    // Add common_mx, which is row-ordered and has dimensions [n + 1 + m', m']
    // Variables u^{term_ind}
    // Indices are (m'+term_ind*m')+1, ..., m'+(term_ind+1)*m'
    term_mx.majorAppendOrthoOrdered(common_mx);
    
    // Add zero_mx for the remaining terms
    // Variables u^{t'} for t' > term_ind
    // Indices are m'+(term_ind+1)*m'+1, ..., m'+|T|*m'
    for (int t = term_ind + 1; t < num_terms; t++) {
      term_mx.rightAppendPackedMatrix(zero_mx);
    }

    // Add term-specific zero_mx for the previous terms
    // Variables u^{t'}_0 for t' < term_ind
    for (int t = 0; t < term_ind; t++) {
      const DisjunctiveTerm& curr_term = disj->terms[t];
      const int mt = curr_term.changed_var.size() + curr_term.ineqs.size();
      CoinPackedMatrix term_zero_mx;
      term_zero_mx.setDimensions(common_mx.getNumRows(), mt);
      term_mx.rightAppendPackedMatrix(term_zero_mx);
    }

    // Create matrix for the negative of the disjunctive term constraints
    // We assume D^t x \ge D^t_0 is stored as >= constraints
    // Variables u^{term_ind}_0
    for (int i = 0; i < (int) term.changed_var.size(); i++) {
      const int col = term.changed_var[i];
      const double coeff = (term.changed_bound[i] <= 0) ? one : negone;
      const double val = term.changed_value[i];
      
      CoinPackedVector row;
      row.insert(col, -1. * coeff);
      row.insert(num_cols, -1. * coeff * val);
      term_mx.appendCol(row);
    }
    for (int i = 0; i < (int) term.ineqs.size(); i++) {
      CoinPackedVector row = -1. * term.ineqs[i].row();
      row.insert(num_cols, -1. * term.ineqs[i].rhs());
      term_mx.appendCol(row);
    }

    // Add term-specific zero_mx for the previous terms
    // Variables u^{t'}_0 for t' > term_ind
    for (int t = term_ind + 1; t < num_terms; t++) {
      const DisjunctiveTerm& curr_term = disj->terms[t];
      const int mt = curr_term.changed_var.size() + curr_term.ineqs.size();
      CoinPackedMatrix term_zero_mx;
      term_zero_mx.setDimensions(common_mx.getNumRows(), mt);
      term_mx.rightAppendPackedMatrix(term_zero_mx);
    }

    mx.bottomAppendPackedMatrix(term_mx);
  } // loop over disjunctive terms

  // Append constraint (4)
  // This is the constraint that the sum of the indicator variables adds up to at most num_cols
  // The variables in matrix are 1 through mprime
  CoinPackedVector numcols_rank_constraint;
  for (int i = 0; i < mprime; i++) {
      numcols_rank_constraint.insert(i+1, one);
  }
  mx.appendRow(numcols_rank_constraint);

  // // Prepare bounds (default >= 0)
  // std::vector<double> colLB(num_cglp_vars, 0.);
  // std::vector<double> colUB(num_cglp_vars, liftingSolver->getInfinity());

  // // Set theta as upper-bounded by one
  // colLB[0] = 0.;
  // colUB[0] = 1.;

  // // Set delta as binary variables
  // for (int i = 0; i < mprime; i++) {
  //   colLB[i+1] = 0.;
  //   colUB[i+1] = 1.;
  // }

  // // Set u^t_i as nonpositive if original row is <=
  // // Set u^t_i as unrestricted in sign if original row is =
  // // Set u^t_i as nonnegative if original row is >=
  // const char* row_sense = solver->getRowSense();
  // for (int term_ind = 0; term_ind < num_terms; term_ind++) {
  //   const int offset = 1 + mprime + term_ind * mprime;
  //   for (int i = 0; i < num_rows; i++) {
  //     if (row_sense[i] == 'L') {
  //       colLB[offset + i] = -1. * liftingSolver->getInfinity();
  //       colUB[offset + i] = 0.;
  //     }
  //     else if (row_sense[i] == 'E') {
  //       colLB[offset + i] = -1. * liftingSolver->getInfinity();
  //       colUB[offset + i] = liftingSolver->getInfinity();
  //     }
  //     else if (row_sense[i] == 'G') {
  //       colLB[offset + i] = 0.;
  //       colUB[offset + i] = liftingSolver->getInfinity();
  //     }
  //   }
  // }

  // Create liftingSolver from mx
  // Defaults (in COIN-OR):
  // colLB: 0   **** This should be -inf when row is <= or =
  // colUB: inf **** This should be 1 for binary variables delta (maybe set automatically by marking them binary?), and 0 for <= rows
  // obj coeff: 0
  // rowSense: >= **** This should be -inf (except when changed above)
  // rhs: 0 **** We want to change this to num_cols for (4)
  // rowRange: 0
  // liftingSolver->loadProblem(mx, colLB.data(), colUB.data(), NULL, NULL, NULL, NULL);
  liftingSolver->loadProblem(mx, NULL, NULL, NULL, NULL, NULL, NULL);
  
  // First num_cols rows are = 0 constraints
  for (int col = 0; col < num_cols; col++) {
    liftingSolver->setRowLower(col, 0);
    liftingSolver->setRowUpper(col, 0);
  }

  if (!use_min_sum_delta) {
    // Constraint (4) is <= num_cols
    liftingSolver->setRowLower(num_cglp_constraints - 1, -1. * liftingSolver->getInfinity());
    // liftingSolver->setRowLower(num_cglp_constraints - 1, num_cols); // for = num_cols constraint, if we want to use that instead
    liftingSolver->setRowUpper(num_cglp_constraints - 1, num_cols);

    // Set objective function as min -\theta
    liftingSolver->setObjCoeff(0, -1.); // min -\theta = max \theta
  } else {
    // Constraint (4') is >= num_cols
    liftingSolver->setRowLower(num_cglp_constraints - 1, num_cols);
    liftingSolver->setRowUpper(num_cglp_constraints - 1, liftingSolver->getInfinity());

    // Set objective function as min sum delta
    for (int i = 0; i < mprime; i++) {
      liftingSolver->setObjCoeff(i+1, 1.);
    }
  }

  // Set theta as upper-bounded by one (lb is already 0)
  liftingSolver->setColUpper(0, 1.);

  // Set u^t_i as nonpositive if original row is <=
  // Set u^t_i as unrestricted in sign if original row is =
  // Variable u^t_i is column 1 [for theta] + m' [for delta] + (t-1) * m' [for u^{t'} for t' < t] + i
  for (int i = 0; i < num_rows; i++) {
    if (row_sense[i] == 'L') {
      for (int term_ind = 0; term_ind < num_terms; term_ind++) {
        const int var_ind = i + 1 + mprime + term_ind * mprime;
        liftingSolver->setColLower(var_ind, -1. * liftingSolver->getInfinity());
        liftingSolver->setColUpper(var_ind, 0.);
      }
    } else if (row_sense[i] == 'E') {
      for (int term_ind = 0; term_ind < num_terms; term_ind++) {
        const int var_ind = i + 1 + mprime + term_ind * mprime;
        liftingSolver->setColLower(var_ind, -1. * liftingSolver->getInfinity());
        liftingSolver->setColUpper(var_ind, liftingSolver->getInfinity());
      }
    } else if (row_sense[i] == 'G') {
      // Nothing to do because default is nonnegative
    }
  }

  // Set binary variables
  for (int i = 0; i < mprime; i++) {
      liftingSolver->setInteger(i+1);
      liftingSolver->setColLower(i+1, 0.);
      liftingSolver->setColUpper(i+1, 1.);
  }

  // Set names of variables
  std::vector<std::string> var_names;
  var_names.push_back("theta");
  for (int i = 0; i < mprime; i++) {
    var_names.push_back("delta_" + std::to_string(i));
  }
  // u^t
  for (int t = 0; t < num_terms; t++) {
    for (int i = 0; i < mprime; i++) {
      var_names.push_back("u_" + std::to_string(t) + "_" + std::to_string(i));
    }
  }
  // u^t_0
  for (int t = 0; t < num_terms; t++) {
    for (int i = 0; i < disj->terms[t].changed_var.size(); i++) {
      var_names.push_back("u0_" + std::to_string(t) + "_" + std::to_string(i));
    }
  }
  assert(liftingSolver->getNumCols() == var_names.size());
  liftingSolver->setColNames(var_names, 0, liftingSolver->getNumCols(), 0);
} /* genRCVMILPFromCut */

void updateRCVMILPFromCut(
    /// [in] RCVMILP to be updated
    OsiSolverInterface* const liftingSolver,
    /// [in] Cut to be added to the RCVMILP
    const OsiRowCut* const cut,
    /// [in] Disjunction that the cut is from
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver,
    /// [in] Log file
    FILE* const log_file) {
  const int num_terms = disj->terms.size();

  // The cut goes into column 0 (corresponding to variable \theta)
  // Let m' be the number of rows of Atilde
  // Then the cut coefficients are inputted into rows 0,...,n-1, then t*(n+m'+1),...,t*(n+m'+1)+n-1 for t = 0,...,num_terms-1
  const int theta_col = 0;
  const int mprime = calculateNumRowsAtilde(disj, solver);

  // 2023-08-04 amk: not sure how to replace column 1, since indices may be different (so CoinPackedMatrix::replaceVector is not an option)
  // Create new CoinPackedMatrix with new cut column
  CoinPackedMatrix old_mx(*liftingSolver->getMatrixByCol()); // this is a copy
  old_mx.deleteCols(1, &theta_col);

  // Prepare rows corresponding to the cut
  CoinPackedVector alpha_beta; // will be [alpha; beta; 0]
  alpha_beta.append(cut->row());
  alpha_beta.insert(solver->getNumCols(), -1. * cut->rhs());
  
  // Duplicate this for each term
  CoinPackedVector alpha_beta_col(alpha_beta);
  for (int term_ind = 1; term_ind < num_terms; term_ind++) {
    const int offset = (solver->getNumCols() + mprime + 1);
    // Modify indices of alpha_beta
    const int num_el = alpha_beta.getNumElements();
    int* ind = alpha_beta.getIndices();
    for (int i = 0; i < num_el; i++) {
      ind[i] += offset;
    }
    alpha_beta_col.append(alpha_beta);
  }
  alpha_beta_col.insert(old_mx.getNumRows()-1, 0.);

  // Create new matrix with first column containing the appropriate number of copies of theta_delta
  // this is not right if theta_col != 0
  if (theta_col != 0) throw std::logic_error("Value of theta_col has to be 0 for now.");
  CoinPackedMatrix new_mx;
  new_mx.setDimensions(old_mx.getNumRows(), 0);
  new_mx.appendCol(alpha_beta_col);
  new_mx.rightAppendPackedMatrix(old_mx);

  // Replace matrix in liftingSolver
  liftingSolver->replaceMatrix(new_mx);

#ifdef TRACE
  {
    // Check that all coefficients were correctly set
    const CoinPackedMatrix* mx = liftingSolver->getMatrixByCol();
    const CoinShallowPackedVector alpha_beta_col_check = mx->getVector(theta_col);  
    const CoinPackedVector cut_row = cut->row();
    const double cut_rhs = cut->rhs();
    for (int term_ind = 0; term_ind < num_terms; term_ind++) {
      const int offset = term_ind * (solver->getNumCols() + mprime + 1);

      // Check coefficients of (alpha,-beta) for this term
      for (int el = 0; el < cut_row.getNumElements(); el++) {
        const int row_ind = offset + cut_row.getIndices()[el];
        const double coeff = cut_row.getElements()[el];
        const double coeff_check = alpha_beta_col_check[row_ind];
        if (!isVal(coeff, coeff_check)) {
          error_msg(
              errorstring,
              "*** ERROR: Coefficients do not match for term %d, row %d: %f != %f.\n",
              term_ind, row_ind, coeff, coeff_check);
          throw std::logic_error(errorstring);
        }
      }

      // Check right-hand side row
      const int row_ind = offset + solver->getNumCols();
      const double coeff = -1. * cut_rhs;
      const double coeff_check = alpha_beta_col_check[row_ind];
      if (!isVal(coeff, coeff_check)) {
          error_msg(
              errorstring,
              "*** ERROR: Constant side do not match for term %d, row %d: %f != %f.\n",
              term_ind, row_ind, coeff, coeff_check);
          throw std::logic_error(errorstring);
      }
    } // loop over terms
  } // check that update was performed correctly
#endif
} /* updateRCVMILPFromCut */

/// @details Solves the RCVMILP and populates the certificate
/// @return 0 if problem solved to optimality or terminated by limit, 1 if problem proved infeasible
int solveRCVMILP(
    /// [in/out] RCVMILP instance
    OsiSolverInterface* const liftingSolver,
    /// [out] Solution to the RCVMILP, where order of variables is theta, delta, {v^t}_{t \in T}
    std::vector<double>& solution,
    /// [in] Parameters for choosing solver
    const StrengtheningParameters::Parameters& params,
    /// [in] Index of the cut for which we are checking the certificate
    const int cut_ind) {
  int return_code = 0;
  const bool write_lp = true && !params.get(StrengtheningParameters::LOGFILE).empty();
  std::string lp_filename_stub = "", LP_EXT = ".lp";
  if (write_lp) {
    // Write to file, using logfile as the output directory
    std::string logdir, logname, in_file_ext;
    parseFilename(logdir, logname, in_file_ext, params.get(StrengtheningParameters::stringParam::LOGFILE), params.logfile);
    std::string instdir, instname;
    parseFilename(instdir, instname, in_file_ext, params.get(StrengtheningParameters::stringParam::FILENAME), params.logfile);
    lp_filename_stub = logdir + "/" + instname + "_cglp_" + stringValue(cut_ind, "%d");
  }

  if (write_lp) {
    std::string lp_filename = lp_filename_stub + "_COIN" + LP_EXT;
    printf("\n## Saving CGLP from cut to file: %s\n", lp_filename.c_str());
    liftingSolver->writeLp(lp_filename.c_str());
  }

  if (use_bb_option(params.get(StrengtheningParameters::intParam::BB_STRATEGY),
      StrengtheningParameters::BB_Strategy_Options::gurobi)) {
#ifdef USE_GUROBI
    GRBModel* m = buildGRBModelFromOsi(liftingSolver, params.logfile);
    GRBModel model = *m;
    setStrategyForBBTestGurobi(params, 0, model);

    if (write_lp) {
      std::string lp_filename = lp_filename_stub + "_GUROBI" + LP_EXT;
      printf("\n## Saving CGLP from cut to file: %s\n", lp_filename.c_str());
      m->write(lp_filename);
    }

    return_code = solveGRBModel(model, params.logfile);
    // doBranchAndBoundWithGurobi(liftingSolver,
    //     params.get(StrengtheningParameters::intParam::BB_STRATEGY),
    //     tmp_bb_info, best_bound, &solution);

    if (return_code == GRB_INFEASIBLE || return_code == GRB_UNBOUNDED || return_code == GRB_INF_OR_UNBD) {
      error_msg(errorstring,
          "solveRCVMILP: Branch and bound with Gurobi did not find an optimal solution for cut %d.\n", cut_ind);
      writeErrorToLog(errorstring, params.logfile);
      throw std::logic_error(errorstring);
    } else {
      return_code = 0;
    }

    // Retrieve solution from model
    saveSolution(solution, model);

    { // DEBUG DEBUG
      GRBVar* vars = model.getVars();

      // Verify col 2
      const double val = -1. * (
          -1 * vars[75].get(GRB_DoubleAttr::GRB_DoubleAttr_X)
          + 8 * vars[76].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 5 * vars[77].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 4 * vars[79].get(GRB_DoubleAttr::GRB_DoubleAttr_X)
          - 2 * vars[80].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 3 * vars[81].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 6 * vars[82].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 7 * vars[83].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 3 * vars[85].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 5 * vars[86].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 8 * vars[87].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 7 * vars[88].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 4 * vars[89].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 1 * vars[90].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 6 * vars[91].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 5 * vars[92].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + 8 * vars[93].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - 2 * vars[94].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          - vars[97].get(GRB_DoubleAttr::GRB_DoubleAttr_X) 
          + vars[124].get(GRB_DoubleAttr::GRB_DoubleAttr_X)
     ) / solution[0];

      printf("value for coefficient on cut for column 2 = %f.\n", val);

      if (vars) {
        delete[] vars;
      }
    } // DEBUG DEBUG

    // Free memory
    if (m) { delete m; }
#endif // USE_GUROBI
  } // USE_GUROBI
  else if (use_bb_option(params.get(StrengtheningParameters::intParam::BB_STRATEGY), 
      StrengtheningParameters::BB_Strategy_Options::cbc)) {
#ifdef USE_CBC
    CbcModel cbc_model(*liftingSolver);
    setIPSolverParameters(&cbc_model, params.get(StrengtheningParameters::VERBOSITY));

    if (write_lp) {
      std::string lp_filename = lp_filename_stub + "_COIN" + LP_EXT;
      printf("\n## Saving CGLP from cut to file: %s\n", lp_filename.c_str());
      liftingSolver->writeLp(lp_filename.c_str());
    }
    
    cbc_model.branchAndBound(params.get(StrengtheningParameters::VERBOSITY));
    if (cbc_model.isProvenOptimal()) {
      // Store found solution in solution
      const double* const cbc_sol = cbc_model.getColSolution();
      const int num_cols = cbc_model.getNumCols();
      solution.assign(cbc_sol, cbc_sol + num_cols);
    } else if (cbc_model.status() == 1) { // time limit, max nodes, or max iters reached
      warning_msg(warnstring, "Cbc stopped with status = 1 for cut %d.\n", cut_ind);
    } else {
      return_code = 1;
    }
#endif // USE_CBC
  } // USE_CBC
  else {
    error_msg(errorstring, "solveRCVMILP: Implementation for solving RCVMILP is available only for Cbc and Gurobi.\n");
    writeErrorToLog(errorstring, params.logfile);
    throw std::logic_error(errorstring);
  }

  return return_code;
} /* solveRCVMILP */

void getCertificateFromRCVMILPSolution(
    /// [out] Certificate of cut (vector of length m + m_t + n)
    CutCertificate& v,
    /// [in] Solution to the RCVMILP, where order of variables is theta, delta (length equal to m' := #calculateNumRowsAtilde), {u^t}_{t \in T}, {u^t_0}_{t \in T}
    const std::vector<double>& solution,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Index of the cut for which we are checking the certificate
    const int cut_ind,
    /// [in] Log file
    FILE* const logfile) {
  std::vector<int> rows, cols;
  std::vector<int> delta;
  delta.reserve(solver->getNumCols());

  const int num_nonbound_constr_tilde = solver->getNumRows() + disj->common_changed_var.size() + disj->common_ineqs.size();

  int num_lb = 0;
  int num_ub = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
    const double lb = solver->getColLower()[col];
    const double ub = solver->getColUpper()[col];

    num_lb += !isInfinity(std::abs(lb));
    num_ub += !isInfinity(std::abs(ub));
  }

  // Check the binary delta variables,
  // which are forced to be nonzero if the corresponding row
  // has a nonzero multiplier for any disjunctive term,
  // through the \delta_i \ge u^t_i constraints (or \delta_i \ge -u^t_i for <= constraints)
  const int delta_var_start = 1;
  for (int row_ind = 0; row_ind < num_nonbound_constr_tilde; row_ind++) {
    const int var = delta_var_start + row_ind;
    // assert(cbc_model.isBinary(var));
    if (isZero(solution[var])) {
      continue;
    }
    delta.push_back(row_ind);
    rows.push_back(row_ind);
  }

  // For the original variable bounds,
  // check if the lb or ub \delta variables are nonzero
  // Throw a warning if they are both nonzero...
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const int var_lb = delta_var_start + num_nonbound_constr_tilde + col_ind;
    const int var_ub = delta_var_start + num_nonbound_constr_tilde + num_lb + col_ind;

    // assert(cbc_model.isBinary(var_lb));
    // assert(cbc_model.isBinary(var_ub));
    
    const bool lb_nonzero = !isZero(solution[var_lb]);
    const bool ub_nonzero = !isZero(solution[var_ub]);
    
    if (lb_nonzero && ub_nonzero) {
      warning_msg(warnstring, "Both lower and upper bound delta variables are nonzero for cut %d, col %d.\n", cut_ind, col_ind);
    }
    else if (lb_nonzero) {
      delta.push_back(var_lb);
      cols.push_back(col_ind);
    }
    else if (ub_nonzero) {
      delta.push_back(var_ub);
      cols.push_back(col_ind);
    }
  } // loop over columns

  // Help to keep track of where v^t variables start in solution, for term t
  const int mtilde = calculateNumRowsAtilde(disj, solver);

  // Loop over terms to set the CutCertificate
  // Recall that certificate v is a vector of length m + m_t + n
  // corresponding to original rows (+ globally-valid inequalities) = m
  // then term-specific entries = m_t
  // then variable bounds = n
  // Meanwhile, the solution vector is ordered as theta, delta, {u^t}_{t \in T}, {u^t_0}_{t \in T}
  const int mprime = calculateNumRowsAtilde(disj, solver);
  v.clear();
  v.resize(disj->num_terms);
  
  int m_t_previous = 0;
  const double theta = solution[0];
  if (!greaterThanVal(theta, 0.)) {
    error_msg(errorstring,
      "getCertificateFromRCVMILPSolution: Theta = %e, which is not positive.\n", theta);
    writeErrorToLog(errorstring, logfile);
    throw std::logic_error(errorstring);
  }
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const int num_term_constr = disj->terms[term_ind].changed_var.size();
    const int rcvmilp_term_uvar_start_ind = delta_var_start + mprime + term_ind * mprime;
    const int rcvmilp_term_u0var_start_ind = delta_var_start + mprime + disj->num_terms * mprime + term_ind * m_t_previous;
    m_t_previous += num_term_constr;
    
    const int size_certificate = num_nonbound_constr_tilde + num_term_constr + solver->getNumCols();
    v[term_ind].clear();
    v[term_ind].resize(size_certificate, 0.0);

    // Set multipliers for original (+ globally-valid) constraints
    for (int row_ind = 0; row_ind < num_nonbound_constr_tilde; row_ind++) {
      const int v_ind = row_ind;
      const int var = rcvmilp_term_uvar_start_ind + row_ind;
      v[term_ind][v_ind] = solution[var] / theta;
    }

    // Set multilpliers for term-specific constraints
    for (int row_ind = 0; row_ind < num_term_constr; row_ind++) {
      const int v_ind = num_nonbound_constr_tilde + row_ind;
      const int var = rcvmilp_term_u0var_start_ind + row_ind;
      v[term_ind][v_ind] = solution[var] / theta;
    }

    // Set multipliers for variable bounds
    for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
      const int v_ind = num_nonbound_constr_tilde + num_term_constr + col_ind;
      const int var_lb = rcvmilp_term_uvar_start_ind + num_nonbound_constr_tilde + col_ind;
      const int var_ub = rcvmilp_term_uvar_start_ind + num_nonbound_constr_tilde + num_lb + col_ind;

      const bool lb_nonzero = !isZero(solution[var_lb]);
      const bool ub_nonzero = !isZero(solution[var_ub]);

      if (lb_nonzero && ub_nonzero) {
        warning_msg(warnstring, "Both lower and upper bound delta variables are nonzero for cut %d, col %d.\n", cut_ind, col_ind);
      }
      else if (lb_nonzero) {
        v[term_ind][v_ind] = solution[var_lb] / theta;
      }
      else if (ub_nonzero) {
        v[term_ind][v_ind] = solution[var_ub] / theta;
      }
    }
  } // loop over terms to set CutCertificate
} /* getCertificateFromRCVMILPSolution */

void analyzeCertificateRegularity(
    /// [out] Rank of submatrix associated to the certificate
    int& certificate_rank,
    /// [out] Number of original (+ globally valid) constraints that have nonzero multipliers in the certificate
    int& num_nonzero_multipliers,
    /// [in] Certificate of cut
    const CutCertificate& v,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Matrix containing original + globally-valid constraints
    const CoinPackedMatrix& Atilde,
    /// [in] Parameters
    const StrengtheningParameters::Parameters& params) {
  std::vector<int> rows;
  std::vector<int> cols;
  for (int row = 0; row < solver->getNumRows() + disj->common_changed_var.size(); row++) {
    for (int term = 0; term < disj->num_terms; term++) {
      if (!isZero(v[term][row])) {
        rows.push_back(row);
        break;
      }
    }
  }
  const int first_col_ind = solver->getNumRows() + disj->common_changed_var.size();
  for (int col = 0; col < solver->getNumCols(); col++) {
    for (int term = 0; term < disj->num_terms; term++) {
      if (!isZero(v[term][first_col_ind + col])) {
        cols.push_back(col);
        break;
      }
    }
  }

  // if (Atilde.getNumRows() == 0) {
  //   std::vector<double> btilde;
  //   prepareAtilde(Atilde, btilde, disj, solver, logfile);
  // }

  certificate_rank = computeRank(&Atilde, rows, cols);
  num_nonzero_multipliers = rows.size() + cols.size();
} /* analyzeCertificateRegularity */

void analyzeCutRegularity(
    /// [out] Certificate of cuts that, in the end, per term, will be of dimension rows + disj term ineqs + cols with indices [cut][term][Farkas multiplier]
    std::vector<CutCertificate>& v,
    /// [out] Rank of submatrix associated to the certificate for each cut
    std::vector<int>& certificate_submx_rank,
    /// [in] Set of cuts that are to be analyzed for regularity
    const OsiCuts& cuts,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Parameters for setting verbosity and logfile
    const StrengtheningParameters::Parameters& params) {
  if (cuts.sizeCuts() == 0) return;

  certificate_submx_rank.clear();
  certificate_submx_rank.resize(cuts.sizeCuts());

  // Prepare Atilde and btilde, which are common to all terms
  // This encompasses the original constraints
  // and all globally-valid inequalities identified as part of the disjunction
  CoinPackedMatrix Atilde;
  std::vector<double> btilde;
  prepareAtilde(Atilde, btilde, disj, solver, params.logfile);

  // Compute rank of Atilde
  const int rank_atilde = computeRank(&Atilde, std::vector<int>(), std::vector<int>());

  // Prepare solver for computing certificate
  OsiSolverInterface* liftingSolver = new SolverInterface;
  setLPSolverParameters(liftingSolver, params.get(StrengtheningParameters::VERBOSITY));
  genRCVMILPFromCut(liftingSolver, cuts.rowCutPtr(0), disj, solver, params);
  // cbc_model->setModelOwnsSolver(false);
  // cbc_model->swapSolver(liftingSolver);
  // cbc_model->setModelOwnsSolver(true); // solver will be deleted with cbc object

  v.clear();
  v.resize(cuts.sizeCuts());
  for (int cut_ind = 0; cut_ind < cuts.sizeCuts(); cut_ind++) {
    if (cut_ind > 0) {
      updateRCVMILPFromCut(liftingSolver, cuts.rowCutPtr(cut_ind), disj, solver, params.logfile);
    }

    // Solve the RCVMILP
    std::vector<double> solution;
    const int return_code = solveRCVMILP(liftingSolver, solution, params, cut_ind);
    if (return_code > 0) {
      error_msg(errorstring, "Solver does not terminate with optimal solution for cut %d.\n", cut_ind);
      writeErrorToLog(errorstring, params.logfile);
      
      if (liftingSolver) { delete liftingSolver; }
      
      throw std::logic_error(errorstring);
    }

    // Retrieve certificate from solution
    getCertificateFromRCVMILPSolution(v[cut_ind], solution, disj, solver, cut_ind, params.logfile);

    // Verify the certificate using dense cut coefficient vector
    const OsiRowCut* cut = cuts.rowCutPtr(cut_ind);
    const int num_elem = cut->row().getNumElements();
    const int* ind = cut->row().getIndices();
    const double* coeff = cut->row().getElements();
    std::vector<double> cut_coeff(solver->getNumCols(), 0.0);
    for (int i = 0; i < num_elem; i++) {
      cut_coeff[ind[i]] = coeff[i];
    }
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      checkCutHelper(cut_coeff, v[cut_ind][term_ind], term_ind, disj, solver, params.logfile);
    }
  } // loop over cuts

  // if (cbc_model) { delete cbc_model; }
  // if (liftingSolver && !cbc_model->modelOwnsSolver()) { delete liftingSolver; }
  if (liftingSolver) { delete liftingSolver; }
} /* analyzeCutRegularity */
