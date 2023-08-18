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
#include "CutHelper.hpp" // getParallelism
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

#ifdef DEBUG
#include "debug.hpp"
#endif

/// @brief Status of RCVMIP after termination of #solveRCVMIP
enum class RCVMIPStatus {
  OPTIMAL_REG = 0,      ///< Optimal regular solution found
  OPTIMAL_IRREG,        ///< Optimal irregular solution found
  OPTIMAL_UNCONVERGED,  ///< Optimal solution found, but not converged
  INFEASIBLE,           ///< Infeasible
  UNBOUNDED,            ///< Unbounded
  INF_OR_UNBD,          ///< Infeasible or unbounded
  SOLVER_LIMIT,         ///< Cutoff, iteration, node, time, solution, user_obj limit
  RCVMIP_ITER_LIMIT,    ///< RCVMIP iteration limit
  ERROR                 ///< Error
};

/// @brief Calculate number of finite lower and upper bounds
int calculateNumFiniteBounds(
  const OsiSolverInterface* const solver,
  int* const num_lb = NULL, int* const num_ub = NULL) {
  if (num_lb) { *num_lb = 0; }
  if (num_ub) { *num_ub = 0; }
  int num_finite_bounds = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
    const bool lb_finite = !isInfinity(std::abs(solver->getColLower()[col]));
    const bool ub_finite = !isInfinity(std::abs(solver->getColUpper()[col]));
    num_finite_bounds += lb_finite + ub_finite;
    
    if (num_lb && lb_finite) {
        (*num_lb)++;
    }

    if (num_ub && ub_finite) {
        (*num_ub)++;
    }
  }

  return num_finite_bounds;
} /* calculateNumFiniteBounds */

/// @brief Get index of variable theta in RCVMIP
int getThetaIndex() {
  return 0;
} /* getThetaIndex */

/// @brief Get index of row \p row_ind in Atilde
/// @return -1 if row is not one of the nonbound constraints in Atilde
int getRowAtildeIndex(
    /// [in] Solver
    const OsiSolverInterface* const solver,
    /// [in] Disjunction
    const Disjunction* const disj,
    /// [in] Row index
    const int row_ind) {
  const int num_orig_rows = (solver) ? solver->getNumRows() : 0;
  const int num_common_rows = (disj) ? disj->common_changed_var.size() + disj->common_ineqs.size() : 0;
  if (row_ind < num_orig_rows + num_common_rows) {
    return row_ind;
  } else {
    return -1;
  }
} /* getRowAtildeIndex */

/// @brief Get index of row in Atilde for lb on column \p col_ind (if negative, just return where these start)
/// @details This is useful in case we later change Atilde, e.g., if we do not have to calculate the number of finite bounds
/// @return -1 if no lb var exists
int getLBAtildeIndex(
    /// [in] Solver
    const OsiSolverInterface* const solver,
    /// [in] Disjunction
    const Disjunction* const disj,
    /// [in] Column index (if negative return where these start)
    const int col_ind = -1,
    /// [in] Number of previous lb delta variables; if < 0, that means we do not have that information available and need to recalc
    const int num_prev_bound = -1) {
  const int num_nonbound_constrs = solver->getNumRows() + disj->common_changed_bound.size() + disj->common_ineqs.size();
  if (col_ind < 0) {
    return num_nonbound_constrs;
  }

  if (isInfinity(std::abs(solver->getColLower()[col_ind]))) {
    return -1;
  }

  if (num_prev_bound >= 0) {
    return num_nonbound_constrs + num_prev_bound;
  } else {
    int num_prev_lb = 0;
    for (int i = 0; i < col_ind; i++) {
      const double val = solver->getColLower()[i];
      num_prev_lb += !isInfinity(std::abs(val));
    }
    return num_nonbound_constrs + num_prev_lb;
  }
} /* getLBAtildeIndex */

/// @brief Get index of row in Atilde for ub on column \p col_ind (if negative, just return where these start)
/// @details This is useful in case we later change Atilde, e.g., if we do not have to calculate the number of finite bounds
/// @return -1 if no lb var exists
int getUBAtildeIndex(
    /// [in] Solver
    const OsiSolverInterface* const solver,
    /// [in] Disjunction
    const Disjunction* const disj,
    /// [in] Column index (if negative return where these start)
    const int col_ind = -1,
    /// [in] Number of previous lb+ub delta variables; if < 0, that means we do not have that information available and need to recalc
    const int num_prev_bound = -1) {
  const int num_nonbound_constrs = solver->getNumRows() + disj->common_changed_bound.size() + disj->common_ineqs.size();
  if (col_ind < 0) {
    if (num_prev_bound >= 0) {
      return num_nonbound_constrs + num_prev_bound;
    } else {
      int num_lb = 0;
      calculateNumFiniteBounds(solver, &num_lb);
      return num_nonbound_constrs + num_lb;
    }
  }

  if (isInfinity(solver->getColUpper()[col_ind])) {
    return -1;
  }

  if (num_prev_bound >= 0) {
    return num_nonbound_constrs + num_prev_bound;
  } else {
    int num_lb = 0, num_prev_ub = 0;
    for (int i = 0; i < solver->getNumCols(); i++) {
      const double lb = solver->getColLower()[i];
      num_lb += !isInfinity(std::abs(lb));

      if (i < col_ind) {
        const double ub = solver->getColUpper()[i];
        num_prev_ub += !isInfinity(std::abs(ub));
      }
    }
    return num_nonbound_constrs + num_lb + num_prev_ub;
  }
} /* getUBAtildeIndex */

/// @brief Get start index of variable u^t for term \p term_ind in RCVMIP
int getUtStartIndex(
    /// [in] Solver
    const OsiSolverInterface* const solver,
    /// [in] Disjunction
    const Disjunction* const disj,
    /// [in] Term index
    const int term_ind,
    /// [in] Number of rows of Aprime
    const int mtilde_in = -1) {
  const int mtilde = (mtilde_in >= 0) ? mtilde_in : calculateNumRowsAtilde(disj, solver);
  return 1 + mtilde + term_ind * mtilde;
} /* getUtStartIndex */

/// @brief Get start index of variable u_0^t for term \p term_ind in RCVMIP
int getUt0StartIndex(
    /// [in] Solver
    const OsiSolverInterface* const solver,
    /// [in] Disjunction
    const Disjunction* const disj,
    /// [in] Term index
    const int term_ind,
    /// [in] Number of rows of Aprime
    const int mtilde_in = -1,
    /// [in] Previous terms u^t_0 variables counted already
    const int prev_mt_in = -1) {
  const int mtilde = (mtilde_in >= 0) ? mtilde_in : calculateNumRowsAtilde(disj, solver);
  const int first_u0_var = 1 + mtilde + disj->num_terms * mtilde;
  int prev_mt = (prev_mt_in < 0) ? 0 : prev_mt_in;
  if (prev_mt_in < 0) {
    for (int t = 0; t < term_ind; t++) {
      prev_mt += disj->terms[t].changed_var.size() + disj->terms[t].ineqs.size();
    }
  }
  return first_u0_var + prev_mt;
} /* getUt0StartIndex */

/// @details Number of rows should be solver->getNumRows() + num_common_rows + number of finite bounds
int calculateNumRowsAtilde(
    /// [in] Disjunction from which to get globally-valid inequalities
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver) {
  // Let mprime be the number of rows of Atilde
  // This is the number of original constraints + number of globally-valid bound changes
  // We also add the number of lower and upper bounds to mprime
  const int num_common_rows = disj->common_changed_var.size() + disj->common_ineqs.size();
  const int mprime_bounds = calculateNumFiniteBounds(solver);
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
void genRCVMIPFromCut(
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
    /// [in] Whether to use original RCVMIP or one with min sum delta objective and sum delta >= n constraint
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
      throw std::logic_error("genRCVMIPFromCut:: Atilde is not col-ordered!");
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
          "genRCVMIPFromCut:: Num rows (mx, actual, predicted) = (Atilde, %d, %d); (common_mx, %d, %d).\n",
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
      row.insert(num_cols, coeff * val);
      term_mx.appendCol(row);
    }
    for (int i = 0; i < (int) term.ineqs.size(); i++) {
      CoinPackedVector row = -1. * term.ineqs[i].row();
      row.insert(num_cols, term.ineqs[i].rhs());
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

  // Now, for each row, see if there is another row that is a multiple of it,
  // and add the constraint that the sum of the indicator variables for the two rows is at most 1
  // Check only the rows corresponding to the original constraints and globally-valid inequalities
  // not the rows corresponding to the variable bounds
  const double last_row_to_check = num_rows + disj->common_changed_bound.size() + disj->common_ineqs.size();
  Atilde.reverseOrdering(); // make it row-ordered
  assert( !Atilde.isColOrdered() );
  for (int r1 = 0; r1 < last_row_to_check; r1++) {
    // Although we do not check the variable bounds for the first index,
    // in case the original or globally-valid constraints are parallel to a variable bound,
    // we check the variable bounds for the second index
    // (This is indeed going to be the case for the globally-valid cosntraints when those are bound improvements derived from strong branching)
    for (int r2 = r1 + 1; r2 < Atilde.getNumRows(); r2++) {
      const CoinPackedVector row1 = Atilde.getVector(r1);
      const CoinPackedVector row2 = Atilde.getVector(r2);
      
      // Two vectors are parallel if u . v / ||u|| ||v|| = 1
      const double how_parallel = getParallelism(row1, row2);
      const bool is_parallel = isVal(how_parallel, 1.) || isVal(how_parallel, -1.);
      if (!is_parallel) {
        continue;
      }
      
      // Add to mx the constraint that the sum of the indicator variables for the two rows is at most 1
      CoinPackedVector constraint;
      constraint.insert(r1 + 1, 1.);
      constraint.insert(r2 + 1, 1.);
      mx.appendRow(constraint);
    } // inner loop over rows
  } // outer loop over rows

  // Add constraints, when appropriate, that variable lb and ub delta vars are at most one 
  int num_lb = 0;
  calculateNumFiniteBounds(solver, &num_lb);
  int lb_ind = 1 + last_row_to_check;
  int ub_ind = lb_ind + num_lb;
  for (int col_ind = 0; col_ind < num_cols; col_ind++) {
    const bool lb_exists = !isInfinity(std::abs(solver->getColLower()[col_ind]));
    const bool ub_exists = !isInfinity(std::abs(solver->getColUpper()[col_ind]));

    if (lb_exists && ub_exists) {
      CoinPackedVector constraint;
      constraint.insert(lb_ind, 1.);
      constraint.insert(ub_ind, 1.);
      mx.appendRow(constraint);
    }
    
    lb_ind += lb_exists;
    ub_ind += ub_exists;
  } // loop over columns to make sure lb + ub delta vars for same column are at most one

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
  
  // Set constant sides for the rows
  // First num_cols rows for each term are = 0 constraints
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const int term_rows_start = term_ind * (num_cols + 1 + mprime);
    for (int col = 0; col < num_cols; col++) {
      liftingSolver->setRowLower(term_rows_start + col, 0);
      liftingSolver->setRowUpper(term_rows_start + col, 0);
    }
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

  // Set rows as <= 1 for extra constraints from parallel rows of Atilde
  for (int extra_row_ind = num_cglp_constraints; extra_row_ind < liftingSolver->getNumRows(); extra_row_ind++) {
    liftingSolver->setRowLower(extra_row_ind, -1. * liftingSolver->getInfinity());
    liftingSolver->setRowUpper(extra_row_ind, 1.);
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

  // Set names of rows
  std::vector<std::string> row_names;
  for (int t = 0; t < num_terms; t++) {
    for (int i = 0; i < num_cols; i++) {
      row_names.push_back("term_" + std::to_string(t) + "_alpha_" + std::to_string(i));
    }
    row_names.push_back("term_" + std::to_string(t) + "_beta");
    for (int i = 0; i < mprime; i++) {
      row_names.push_back("term_" + std::to_string(t) + "_delta_" + std::to_string(i));
    }
  }
  row_names.push_back("sum_delta");
  for (int extra_row_ind = num_cglp_constraints; extra_row_ind < liftingSolver->getNumRows(); extra_row_ind++) {
    row_names.push_back("par_rows_" + std::to_string(extra_row_ind - num_cglp_constraints));
  }
  assert(liftingSolver->getNumRows() == row_names.size());
  liftingSolver->setRowNames(row_names, 0, liftingSolver->getNumRows(), 0);
} /* genRCVMIPFromCut */

int computeRankOfRCVMIPSolution(
    /// [out] Indices of nonzero multipliers
    std::vector<int>& delta_var_inds,
    /// [in] Solution to RCVMIP, where variables are theta then delta then {u^t}_{t \in T} then {u^t_0}_{t \in T}
    const double* const solution,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Matrix containing original + globally-valid constraints
    const CoinPackedMatrix& Atilde,
    /// [in] Parameters
    const StrengtheningParameters::Parameters& params,
    /// [in] Cut index
    const int cut_ind) {
  std::vector<int> rows, cols;

  delta_var_inds.clear();
  delta_var_inds.reserve(solver->getNumCols());

  const int num_nonbound_constr_tilde = solver->getNumRows() + disj->common_changed_var.size() + disj->common_ineqs.size();
  const int mtilde = calculateNumRowsAtilde(disj, solver);

  int num_lb = 0;
  calculateNumFiniteBounds(solver, &num_lb);

  // Check the binary delta variables,
  // which are forced to be nonzero if the corresponding row
  // has a nonzero multiplier for any disjunctive term,
  // through the \delta_i \ge u^t_i constraints (or \delta_i \ge -u^t_i for <= constraints)
  const int delta_var_start = 1;
  for (int row_ind = 0; row_ind < num_nonbound_constr_tilde; row_ind++) {
    const int atilde_row_ind = getRowAtildeIndex(solver, disj, row_ind);
    assert( atilde_row_ind >= 0 );
    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;
    
    if (isZero(solution[rcvmip_delta_ind])) {
      continue;
    }

    // It is possible that the delta var is 1 but the corresponding row is not used
    // Check the multipliers for each term
    bool is_used = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const int rcvmip_var_ind = getUtStartIndex(solver, disj, term_ind, mtilde) + atilde_row_ind;
      if (!isZero(solution[rcvmip_var_ind])) {
        is_used = true;
        break;
      }
    }

    if (!is_used) {
      continue;
    }

    delta_var_inds.push_back(rcvmip_delta_ind);
    rows.push_back(row_ind);
  }

  // For the original variable bounds,
  // check if the lb or ub \delta variables are nonzero
  // Throw an error if they are both nonzero...
  // (cannot be independent rows, so we have already said this is not possible with a constraint)
  int num_prev_bound = 0;
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const int atilde_row_ind = getLBAtildeIndex(solver, disj, col_ind, num_prev_bound);
    if (atilde_row_ind < 0) {
      // There are no RCVMIP variables associated to this column
      continue;
    } else {
      num_prev_bound += 1;
    }
    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;

    if (isZero(solution[rcvmip_delta_ind])) {
      continue;
    }

    // It is possible that the delta var is 1 but the corresponding row is not used
    // Check the multipliers for each term
    bool is_used = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const int rcvmip_var_ind = getUtStartIndex(solver, disj, term_ind, mtilde) + atilde_row_ind;
      if (!isZero(solution[rcvmip_var_ind])) {
        is_used = true;
        break;
      }
    }

    if (!is_used) {
      continue;
    }

    delta_var_inds.push_back(rcvmip_delta_ind);
    cols.push_back(col_ind);
  } // loop over column for lb

  // Repeat for ub
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const int atilde_row_ind = getUBAtildeIndex(solver, disj, col_ind, num_prev_bound);
    if (atilde_row_ind < 0) {
      // There are no RCVMIP variables associated to this column
      continue;
    } else {
      num_prev_bound += 1;
    }
    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;

    if (isZero(solution[rcvmip_delta_ind])) {
      continue;
    }

    // It is possible that the delta var is 1 but the corresponding row is not used
    // Check the multipliers for each term
    bool is_used = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const int rcvmip_var_ind = getUtStartIndex(solver, disj, term_ind, mtilde) + atilde_row_ind;
      if (!isZero(solution[rcvmip_var_ind])) {
        is_used = true;
        break;
      }
    }

    if (!is_used) {
      continue;
    }

    delta_var_inds.push_back(rcvmip_delta_ind);
    cols.push_back(col_ind);
  } // loop over column for lb

  const int certificate_rank = computeRank(&Atilde, rows, cols);
  return certificate_rank;
} /* computeRankOfRCVMIPSolution */

#ifdef USE_GUROBI
/// @brief Add warm start if existing certificate provided
void setRCVMIPStart(
    /// [in/out] Gurobi RCVMIP instance 
    GRBModel* const model,
    /// [out] Solution derived from certificate, having indices for variables theta then delta (m') then {u^t}_t then {u^t_0}_t
    std::vector<double>& solution,
    /// [in] Certificate of cut (in-version used to set MIP start), [term][Farkas multiplier]; per term, m + m_t + n indices correspond to rows (including globally-valid constraints) + disj term ineqs + cols
    const CutCertificate& v,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver) {    
  if (v.size() == 0) return;

  // Recall that the certificate v is a vector of length m + m_t + n
  // The delta variable for constraint row_ind will be set to one if v[t][row_ind] > 0 for some t \in num_terms
  const int num_nonbound_constrs = solver->getNumRows() + disj->common_changed_var.size() + disj->common_ineqs.size();
  assert( v[0].size() == num_nonbound_constrs + disj->terms[0].changed_var.size() + solver->getNumCols() );
  const int mtilde = calculateNumRowsAtilde(disj, solver);
  // assert( Atilde.getNumRows() == mtilde );

  // If need to scale, find largest value in v
  const bool SHOULD_SCALE = true;
  double scale = 1.;
  if (SHOULD_SCALE) {
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      for (int i = 0; i < v[term_ind].size(); i++) {
        if (std::abs(v[term_ind][i]) > scale) {
          scale = std::abs(v[term_ind][i]);
        }
      }
    }
  }

  // Get variables (delete at end)
  GRBVar* const vars = model->getVars(); // vector of variables, theta then delta (length m') then {u^t}_t then {u^t_0}_t
  int num_vars_set = 0;

  // Set theta
  vars[getThetaIndex()].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 1. / scale);
  num_vars_set++;

  // Set delta and u^t variables for original + globally-valid constraints
  const int delta_var_start = 1;
  for (int row_ind = 0; row_ind < num_nonbound_constrs; row_ind++) {
    const int atilde_row_ind = getRowAtildeIndex(solver, disj, row_ind);

    bool set_delta = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      // Get RCVMIP variable index for this term and constraint
      const int rcvmip_term_var_start_ind = getUtStartIndex(solver, disj, term_ind, mtilde);
      const int curr_rcvmip_ind = rcvmip_term_var_start_ind + atilde_row_ind;
      
      // Get the index inside of the certificate v
      const int term_cert_var_ind = row_ind;
      const double v_val = v[term_ind][term_cert_var_ind];

      vars[curr_rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, v_val / scale);
      num_vars_set++;

      if (!set_delta && !isZero(v_val)) {
        set_delta = true;
      }
    }

    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;
    if (set_delta) {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 1.);
    } else {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 0.);
    }
    num_vars_set++;
  }

  // Set delta and u^t variables for variable lower bounds
  int num_prev_bound = 0;
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const int atilde_row_ind = getLBAtildeIndex(solver, disj, col_ind, num_prev_bound);
    if (atilde_row_ind < 0) {
      // There are no RCVMIP variables associated to this column
      continue;
    } else {
      num_prev_bound += 1;
    }

    bool set_delta = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      // Get RCVMIP variable index for this term and constraint
      const int rcvmip_term_var_start_ind = getUtStartIndex(solver, disj, term_ind, mtilde);
      const int curr_rcvmip_ind = rcvmip_term_var_start_ind + atilde_row_ind;
      
      // Get the index inside of the certificate v
      const int term_cert_var_ind = num_nonbound_constrs + disj->terms[term_ind].changed_var.size() + col_ind;
      const double v_val = v[term_ind][term_cert_var_ind];

      // If positive then this is corresponds to a multiplier on the lower bound
      if (greaterThanVal(v_val, 0.)) {
        vars[curr_rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, v_val / scale);
        set_delta = true;
      } else {
        vars[curr_rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 0.);
      }
      num_vars_set++;
    } // loop over terms for lb
    
    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;
    if (set_delta) {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 1.);
    } else {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 0.);
    }
    num_vars_set++;
  } // loop over columns for lb
  
  // Set delta and u^t variables for variable upper bounds
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const int atilde_row_ind = getUBAtildeIndex(solver, disj, col_ind, num_prev_bound);
    if (atilde_row_ind < 0) {
      // There are no RCVMIP variables associated to this column
      continue;
    } else {
      num_prev_bound += 1;
    }

    bool set_delta = false;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      // Get RCVMIP variable index for this term and constraint
      const int rcvmip_term_var_start_ind = getUtStartIndex(solver, disj, term_ind, mtilde);
      const int curr_rcvmip_ind = rcvmip_term_var_start_ind + atilde_row_ind;
      
      // Get the index inside of the certificate v
      const int term_cert_var_ind = num_nonbound_constrs + disj->terms[term_ind].changed_var.size() + col_ind;
      const double v_val = v[term_ind][term_cert_var_ind];

      // If positive then this is corresponds to a multiplier on the lower bound
      if (lessThanVal(v_val, 0.)) {
        vars[curr_rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, v_val / scale);
        set_delta = true;
      } else {
        vars[curr_rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 0.);
      }
      num_vars_set++;
    } // loop over terms for lb
    
    const int rcvmip_delta_ind = delta_var_start + atilde_row_ind;
    if (set_delta) {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 1.);
    } else {
      vars[rcvmip_delta_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, 0.);
    }
    num_vars_set++;
  } // loop over columns for ub

  // Set u^t_0 variables
  int prev_mt = 0;
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const int rcvmip_term_var_start_ind = getUt0StartIndex(solver, disj, term_ind, mtilde, prev_mt);
    const int curr_mt = disj->terms[term_ind].changed_var.size() + disj->terms[term_ind].ineqs.size();
    prev_mt += curr_mt;

    for (int i = 0; i < curr_mt; i++) {
      const int rcvmip_ind = rcvmip_term_var_start_ind + i;
      const int cert_ind = num_nonbound_constrs + i;
      const double v_val = v[term_ind][cert_ind];
      vars[rcvmip_ind].set(GRB_DoubleAttr::GRB_DoubleAttr_Start, v_val / scale);
      num_vars_set++;
    }
  }

  model->update();

  const int num_rcvmip_vars = model->get(GRB_IntAttr::GRB_IntAttr_NumVars);
  assert( num_rcvmip_vars == num_vars_set );

  solution.resize(num_rcvmip_vars);
  for (int i = 0; i < num_rcvmip_vars; i++) {
    solution[i] = vars[i].get(GRB_DoubleAttr::GRB_DoubleAttr_Start);
  }

  if (vars) { delete[] vars; }
} /* setRCVMIPStart */

/// @brief Append row to Gurobi \p model restricting sum of delta variables to be <= \p rank
void addRankConstraint(
    /// [in/out] Gurobi RCVMIP instance 
    GRBModel* const model,
    /// [in] Indices of delta variables
    const std::vector<int>& delta,
    /// [in] Rank of solution
    const int rank,
    /// [in] Iteration number
    const int iter_num = -1) {
  const int num_nonzero_multipliers = delta.size();
  if (num_nonzero_multipliers <= rank) {
    return;
  }

  // Create new GRBLinExpr
  GRBLinExpr expr;

  std::vector<double> coeffs(num_nonzero_multipliers, 1.);
  std::vector<GRBVar> vars;

  GRBVar* all_vars = model->getVars(); // An array of all variables in the model. Note that this array is heap-allocated, and must be returned to the heap by the user.

  for (int i = 0; i < num_nonzero_multipliers; i++) {
    const int var_ind = delta[i];
    vars.push_back(all_vars[var_ind]);
  }

  expr.addTerms(coeffs.data(), vars.data(), num_nonzero_multipliers);

  // Add constraint
  std::string constr_name = iter_num >= 0 ? "rank_" + std::to_string(iter_num) : "";
  model->addConstr(expr, GRB_LESS_EQUAL, rank, constr_name);

  // Update model
  model->update();

  if (all_vars) {
    delete[] all_vars;
  }
} /* addRankConstraint (Gurobi) */

void updateRCVMIPFromCut(
    /// [in] RCVMIP to be updated
    GRBModel* const model,
    /// [in] Cut to be added to the RCVMIP
    const OsiRowCut* const cut,
    /// [in] Disjunction that the cut is from
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver,
    /// [in] Log file
    FILE* const log_file) {
  // Create dense vector from cut
  std::vector<double> alpha_beta(solver->getNumCols() + 1, 0.);

  const CoinPackedVector cut_row = cut->row();
  for (int i = 0; i < cut_row.getNumElements(); i++) {
    const int col_ind = cut_row.getIndices()[i];
    const double coeff = cut_row.getElements()[i];
    alpha_beta[col_ind] = coeff;
  }
  alpha_beta[solver->getNumCols()] = -1. * cut->rhs();

  // The cut goes into column 0 (corresponding to variable \theta)
  // Let m' be the number of rows of Atilde
  // Then the cut coefficients are inputted into rows 0,...,n-1, then t*(n+m'+1),...,t*(n+m'+1)+n-1 for t = 0,...,num_terms-1
  const int num_terms = disj->terms.size();
  const int theta_col = 0;
  const int mprime = calculateNumRowsAtilde(disj, solver);

  GRBConstr* constrs = model->getConstrs(); // An array of all linear constraints in the model. Note that this array is heap-allocated, and must be returned to the heap by the user.
  GRBVar* vars = model->getVars();
  GRBVar theta_var = vars[theta_col];

  // Set each term's coefficients with new cut
  for (int term_ind = 0; term_ind < num_terms; term_ind++) {
    const int offset = term_ind * (solver->getNumCols() + 1 + mprime);

    for (int col_ind = 0; col_ind <= solver->getNumCols(); col_ind++) {
      const int constr_ind = col_ind + offset;
      const double newvalue = alpha_beta[col_ind];
      model->chgCoeff(constrs[constr_ind], theta_var, newvalue);
    }
  }
  
  model->update();

#ifdef TRACE
  {
    // Check coefficients of (alpha,-beta) for each term
    for (int term_ind = 0; term_ind < num_terms; term_ind++) {
      const int offset = term_ind * (solver->getNumCols() + 1 + mprime);

      for (int col_ind = 0; col_ind <= solver->getNumCols(); col_ind++) {
        const int constr_ind = col_ind + offset;
        const double coeff = alpha_beta[col_ind];
        const double coeff_check = model->getCoeff(constrs[constr_ind], theta_var);
        if (!isVal(coeff, coeff_check)) {
          error_msg(
              errorstring,
              "*** ERROR: Coefficients do not match for term %d, col %d/%d (last = rhs): %f != %f.\n",
              term_ind, col_ind, solver->getNumCols(), coeff, coeff_check);
          throw std::logic_error(errorstring);
        }
      }
    } // loop over terms
  } // check that update was performed correctly
#endif

  if (constrs) { delete[] constrs; }
  if (vars) { delete[] vars; }
} /* updateRCVMIPFromCut (Gurobi) */

/// @details Solves the RCVMIP and populates the certificate
/// @return 0 if problem solved to optimality, 1 if problem proved infeasible, 2 if terminated by limit
RCVMIPStatus solveRCVMIP(
    /// [in/out] RCVMIP instance
    GRBModel* const model,
    /// [out] Solution to the RCVMIP, where order of variables is theta, delta, {v^t}_{t \in T}
    std::vector<double>& solution,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Matrix containing original + globally-valid constraints
    const CoinPackedMatrix& Atilde,
    /// [in] Parameters for choosing solver
    const StrengtheningParameters::Parameters& params,
    /// [in] Index of the cut for which we are checking the certificate
    const int cut_ind,
    /// [out] Number of iterations taken
    int& num_iters) {
  RCVMIPStatus return_code = RCVMIPStatus::ERROR;

#ifdef TRACE
  const bool write_lp = true && !params.get(StrengtheningParameters::LOGFILE).empty();
#else
  const bool write_lp = false && !params.get(StrengtheningParameters::LOGFILE).empty();
#endif
  std::string lp_filename_stub = "", LP_EXT = ".lp";
  if (write_lp) {
    // Write to file, using logfile as the output directory
    std::string logdir, logname, in_file_ext;
    parseFilename(logdir, logname, in_file_ext, params.get(StrengtheningParameters::stringParam::LOGFILE), params.logfile);
    std::string instdir, instname;
    parseFilename(instdir, instname, in_file_ext, params.get(StrengtheningParameters::stringParam::FILENAME), params.logfile);
    lp_filename_stub = logdir + "/" + instname + "_cglp_" + stringValue(cut_ind, "%d") + "_GUROBI";
  }

  bool reached_feasibility = false;
  num_iters = 0;
  const int MAX_ITERS = params.get(StrengtheningParameters::intParam::RCVMIP_MAX_ITERS);
  printf("\n## solveRCVMIP (Gurobi): Solving CGLP from cut %d. ##\n", cut_ind);
  while (!reached_feasibility && num_iters < MAX_ITERS) {
    const int grb_return_code = solveGRBModel(*model, params.logfile);
    num_iters++;

    if (grb_return_code == GRB_OPTIMAL) {
      return_code = RCVMIPStatus::OPTIMAL_UNCONVERGED;
    }
    else if (grb_return_code == GRB_CUTOFF 
        || grb_return_code == GRB_ITERATION_LIMIT
        || grb_return_code == GRB_NODE_LIMIT
        || grb_return_code == GRB_TIME_LIMIT
        || grb_return_code == GRB_SOLUTION_LIMIT
        || grb_return_code == GRB_USER_OBJ_LIMIT) {
      return_code = RCVMIPStatus::SOLVER_LIMIT;
      break;
    }
    else if (grb_return_code == GRB_INFEASIBLE) {
      return_code = RCVMIPStatus::INFEASIBLE;
      break;
    }
    else if (grb_return_code == GRB_UNBOUNDED) {
      return_code = RCVMIPStatus::UNBOUNDED;
      break;
    }
    else if (grb_return_code == GRB_INF_OR_UNBD) {
      return_code = RCVMIPStatus::INF_OR_UNBD;
      break;
    }
    else {
      // Other options include GRB_LOADED, GRB_INTERRUPTED, GRB_NUMERIC, GRB_SUBOPTIMAL, GRB_INPROGRESS, GRB_WORK_LIMIT
      return_code = RCVMIPStatus::ERROR;
      break;
    }

    // Retrieve solution from model
    saveSolution(solution, *model);

    // Check rank of solution vs number of nonzero multipliers
    std::vector<int> delta_var_inds;
    const int certificate_rank = computeRankOfRCVMIPSolution(delta_var_inds, solution.data(), disj, solver, Atilde, params, cut_ind);
    if (certificate_rank < (int) delta_var_inds.size()) {
      addRankConstraint(model, delta_var_inds, certificate_rank, num_iters);
    } else {
      reached_feasibility = true;
    }

#ifdef TRACE
    {
      // Print value of theta variable (should this decrease across rounds?)
      const double theta_val = solution.size() > 0 ? solution[0] : -1.;
      printf("Iter %d/%d: theta = %f\tcert_rank = %d\tcert_size = %d\n", num_iters, MAX_ITERS, theta_val, certificate_rank, (int) delta_var_inds.size());
    }
#endif
  } // iterate while !reached_feasibility
  printf("solveRCVMIP (Gurobi): Terminated in %d / %d iterations. Reached feasibility: %d.\n", num_iters, MAX_ITERS, reached_feasibility);

  if (reached_feasibility) {
    const double theta_val = solution.size() > 0 ? solution[0] : -1.;
    if (greaterThanVal(theta_val, 0.)) {
      return_code = RCVMIPStatus::OPTIMAL_REG; // May not be feasible, due to missing rank constraints
    }
    else if (isVal(theta_val, 0.)) {
      return_code = RCVMIPStatus::OPTIMAL_IRREG;
    } else {
      error_msg(errorstring, "*** ERROR: Theta value is negative: %f.\n", theta_val);
      writeErrorToLog(errorstring, params.logfile);
      throw std::logic_error(errorstring);
    }
  }

  if (num_iters >= MAX_ITERS && !reached_feasibility) {
    return_code = RCVMIPStatus::RCVMIP_ITER_LIMIT;
  }
  
  if (write_lp) {
    const std::string lp_filename = lp_filename_stub + LP_EXT;
    printf("Saving RCVMIP (Gurobi) to file: %s\n", lp_filename.c_str());
    model->write(lp_filename);
  }

  return return_code;
} /* solveRCVMIP (Gurobi) */
#endif // USE_GUROBI

#ifdef USE_CBC
/// @brief Append row to \p model restricting sum of delta variables to be <= \p rank
void addRankConstraint(
    /// [in/out] CbcModel RCVMIP instance
    CbcModel* const model,
    /// [in] Indices of delta variables
    const std::vector<int>& delta,
    /// [in] Rank of solution
    const int rank,
    /// [in] Iteration number
    const int iter_num = -1) {
  const int num_nonzero_multipliers = delta.size();
  if (num_nonzero_multipliers <= rank) {
    return;
  }

  // Add new constraint to Cbc model to restrict sum of delta variables to be <= rank
  // Create new CoinPackedVector
  CoinPackedVector row;
  row.reserve(num_nonzero_multipliers);
  for (int i = 0; i < num_nonzero_multipliers; i++) {
    const int var_ind = delta[i];
    row.insert(var_ind, 1.);
  }

  // Add constraint
  const double lb = -1. * model->getInfinity();
  const double ub = rank;
  // const char sense = 'L';
  const std::string name = iter_num > 0 ? "rank_" + std::to_string(iter_num) : "";

  model->solver()->addRow(row, lb, ub, name);
} /* addRankConstraint (Cbc) */

void updateRCVMIPFromCut(
    /// [in] RCVMIP to be updated
    CbcModel* const cbc_model,
    /// [in] Cut to be added to the RCVMIP
    const OsiRowCut* const cut,
    /// [in] Disjunction that the cut is from
    const Disjunction* const disj,
    /// [in] Original solver
    const OsiSolverInterface* const solver,
    /// [in] Log file
    FILE* const log_file) {
  OsiSolverInterface* liftingSolver = cbc_model->solver();

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

  OsiSolverInterface* oldSolver = cbc_model->swapSolver(liftingSolver);
  if (oldSolver && cbc_model->modelOwnsSolver()) {
    delete oldSolver;
  }
} /* updateRCVMIPFromCut (Cbc) */

/// @details Solves the RCVMIP and populates the certificate
/// @return 0 if problem solved to optimality or terminated by limit, 1 if problem proved infeasible
RCVMIPStatus solveRCVMIP(
    /// [in/out] RCVMIP instance
    CbcModel* const cbc_model,
    /// [out] Solution to the RCVMIP, where order of variables is theta, delta, {v^t}_{t \in T}
    std::vector<double>& solution,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Matrix containing original + globally-valid constraints
    const CoinPackedMatrix& Atilde,
    /// [in] Parameters for choosing solver
    const StrengtheningParameters::Parameters& params,
    /// [in] Index of the cut for which we are checking the certificate
    const int cut_ind,
    /// [out] Number of iterations taken
    int& num_iters) {
  RCVMIPStatus return_code = RCVMIPStatus::ERROR;
  return return_code; // not yet implemented

#ifdef TRACE
  const bool write_lp = true && !params.get(StrengtheningParameters::LOGFILE).empty();
#else
  const bool write_lp = false && !params.get(StrengtheningParameters::LOGFILE).empty();
#endif
  std::string lp_filename_stub = "", LP_EXT = ".lp";
  if (write_lp) {
    // Write to file, using logfile as the output directory
    std::string logdir, logname, in_file_ext;
    parseFilename(logdir, logname, in_file_ext, params.get(StrengtheningParameters::stringParam::LOGFILE), params.logfile);
    std::string instdir, instname;
    parseFilename(instdir, instname, in_file_ext, params.get(StrengtheningParameters::stringParam::FILENAME), params.logfile);
    lp_filename_stub = logdir + "/" + instname + "_cglp_" + stringValue(cut_ind, "%d") + "_COIN";
  }

  bool reached_feasibility = false;
  num_iters = 0;
  const int MAX_ITERS = params.get(StrengtheningParameters::intParam::RCVMIP_MAX_ITERS);
  while (!reached_feasibility && num_iters < MAX_ITERS) {
    printf("\n## solveRCVMIP (Cbc): Solving CGLP (iter %d) from cut %d. ##\n", num_iters, cut_ind);

    if (write_lp) {
      const std::string lp_filename = lp_filename_stub + LP_EXT;
      printf("File: %s\n", lp_filename.c_str());
      cbc_model->referenceSolver()->writeLp(lp_filename.c_str());
    }
    
    cbc_model->branchAndBound(params.get(StrengtheningParameters::VERBOSITY));
    num_iters++;
    
    if (cbc_model->status() == 1) { // time limit, max nodes, or max iters reached
      // warning_msg(warnstring, "solveRCVMIP (Cbc): Cbc stopped with status = 1 for cut %d.\n", cut_ind);
      return_code = RCVMIPStatus::SOLVER_LIMIT;
      break;
    }
    else if (cbc_model->status() == 0 && cbc_model->isProvenInfeasible()) {
      return_code = RCVMIPStatus::INFEASIBLE;
      break;
    }
    else if (cbc_model->status() == 0 && cbc_model->isProvenDualInfeasible()) {
      return_code = RCVMIPStatus::UNBOUNDED;
      break;
    }
    else {
      return_code = RCVMIPStatus::ERROR;
      break;
    }

    // Retrieve solution from model
    const double* const cbc_sol = cbc_model->getColSolution();
    const int num_cols = cbc_model->getNumCols();
    solution.assign(cbc_sol, cbc_sol + num_cols);

    // Check rank of solution vs number of nonzero multipliers
    std::vector<int> delta;
    const int certificate_rank = computeRankOfRCVMIPSolution(delta, solution.data(), disj, solver, Atilde, params, cut_ind);
    if (certificate_rank < (int) delta.size()) {
      addRankConstraint(cbc_model, delta, certificate_rank, num_iters);
    } else {
      reached_feasibility = true;
    }
  } // iterate while !reached_feasibility
  printf("solveRCVMIP (Cbc): Terminated in %d / %d iterations. Reached feasibility: %d.\n", num_iters, MAX_ITERS, reached_feasibility);

  if (reached_feasibility) {
    const double theta_val = solution.size() > 0 ? solution[0] : -1.;
    if (greaterThanVal(theta_val, 0.)) {
      return_code = RCVMIPStatus::OPTIMAL_REG; // May not be feasible, due to missing rank constraints
    }
    else if (isVal(theta_val, 0.)) {
      return_code = RCVMIPStatus::OPTIMAL_IRREG;
    } else {
      error_msg(errorstring, "*** ERROR: Theta value is negative: %f.\n", theta_val);
      writeErrorToLog(errorstring, params.logfile);
      throw std::logic_error(errorstring);
    }
  }

  if (num_iters >= MAX_ITERS && !reached_feasibility) {
    return_code = RCVMIPStatus::RCVMIP_ITER_LIMIT;
  }

  return return_code;
} /* solveRCVMIP (Cbc) */
#endif // USE_CBC

void getCertificateFromRCVMIPSolution(
    /// [out] Certificate of cut, [term][Farkas multiplier]; per term, m + m_t + n indices correspond to rows + disj term ineqs + cols
    CutCertificate& v,
    /// [in] Solution to the RCVMIP, where order of variables is theta, delta (length equal to m' := #calculateNumRowsAtilde), {u^t}_{t \in T}, {u^t_0}_{t \in T}
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
  calculateNumFiniteBounds(solver, &num_lb);

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
  int lb_ind = delta_var_start + num_nonbound_constr_tilde;
  int ub_ind = delta_var_start + num_nonbound_constr_tilde + num_lb;
  for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
    const bool lb_exists = !isInfinity(std::abs(solver->getColLower()[col_ind]));
    const bool ub_exists = !isInfinity(std::abs(solver->getColUpper()[col_ind]));

    const bool lb_nonzero = lb_exists && !isZero(solution[lb_ind]);
    const bool ub_nonzero = ub_exists && !isZero(solution[ub_ind]);
    
    if (lb_nonzero && ub_nonzero) {
      error_msg(errorstring, "Delta variables %d (lb) and %d (ub) are nonzero for cut %d, col %d.\n", lb_ind, ub_ind, cut_ind, col_ind);
      writeErrorToLog(errorstring, logfile);
      throw std::logic_error(errorstring);
    }
    else if (lb_nonzero) {
      delta.push_back(lb_ind);
      cols.push_back(col_ind);
    }
    else if (ub_nonzero) {
      delta.push_back(ub_ind);
      cols.push_back(col_ind);
    }

    lb_ind += lb_exists;
    ub_ind += ub_exists;
  } // loop over columns

  // Help to keep track of where v^t variables start in solution, for term t
  // const int mtilde = calculateNumRowsAtilde(disj, solver);

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
      "getCertificateFromRCVMIPSolution: Theta = %e, which is not positive.\n", theta);
    writeErrorToLog(errorstring, logfile);
    throw std::logic_error(errorstring);
  }
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const int num_term_constr = disj->terms[term_ind].changed_var.size();
    const int RCVMIP_term_uvar_start_ind = delta_var_start + mprime + term_ind * mprime;
    const int RCVMIP_term_u0var_start_ind = delta_var_start + mprime + disj->num_terms * mprime + m_t_previous;
    m_t_previous += num_term_constr;
    
    const int size_certificate = num_nonbound_constr_tilde + num_term_constr + solver->getNumCols();
    v[term_ind].clear();
    v[term_ind].resize(size_certificate, 0.0);

    // Set multipliers for original (+ globally-valid) constraints
    for (int row_ind = 0; row_ind < num_nonbound_constr_tilde; row_ind++) {
      const int v_ind = row_ind;
      const int var = RCVMIP_term_uvar_start_ind + row_ind;
      v[term_ind][v_ind] = solution[var] / theta;
    }

    // Set multilpliers for term-specific constraints
    for (int row_ind = 0; row_ind < num_term_constr; row_ind++) {
      const int v_ind = num_nonbound_constr_tilde + row_ind;
      const int var = RCVMIP_term_u0var_start_ind + row_ind;
      v[term_ind][v_ind] = solution[var] / theta;
    }

    // Set multipliers for variable bounds
    int term_lb_ind = RCVMIP_term_uvar_start_ind + num_nonbound_constr_tilde;
    int term_ub_ind = RCVMIP_term_uvar_start_ind + num_nonbound_constr_tilde + num_lb;
    for (int col_ind = 0; col_ind < solver->getNumCols(); col_ind++) {
      const bool lb_exists = !isInfinity(std::abs(solver->getColLower()[col_ind]));
      const bool ub_exists = !isInfinity(std::abs(solver->getColUpper()[col_ind]));

      const int v_ind = num_nonbound_constr_tilde + num_term_constr + col_ind;

      const bool lb_nonzero = lb_exists && !isZero(solution[term_lb_ind]);
      const bool ub_nonzero = ub_exists && !isZero(solution[term_ub_ind]);

      if (lb_nonzero && ub_nonzero) {
        // warning_msg(warnstring, "Both lower and upper bound delta variables are nonzero for cut %d, col %d.\n", cut_ind, col_ind);
        error_msg(errorstring, "Delta variables %d (lb) and %d (ub) are nonzero for cut %d, col %d.\n", lb_ind, ub_ind, cut_ind, col_ind);
        writeErrorToLog(errorstring, logfile);
        throw std::logic_error(errorstring);
      }
      else if (lb_nonzero) {
        v[term_ind][v_ind] = solution[term_lb_ind] / theta;
      }
      else if (ub_nonzero) {
        v[term_ind][v_ind] = -1. * solution[term_ub_ind] / theta;
      }

      term_lb_ind += lb_exists;
      term_ub_ind += ub_exists;
    } // loop over columns
  } // loop over terms to set CutCertificate
} /* getCertificateFromRCVMIPSolution */

/// @return Regularity status (-1 if certificate rank is less than number of nonzero multipliers, 0 if equal, 1 if greater)
RegularityStatus analyzeCertificateRegularity(
    /// [out] Rank of submatrix associated to the certificate
    int& certificate_rank,
    /// [out] Number of original (+ globally valid) constraints that have nonzero multipliers in the certificate
    int& num_nonzero_multipliers,
    /// [in] Certificate of cut, [term][Farkas multiplier]; per term, m + m_t + n indices correspond to rows (including globally-valid constraints) + disj term ineqs + cols
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

  for (int col = 0; col < solver->getNumCols(); col++) {
    for (int term = 0; term < disj->num_terms; term++) {
      const int first_col_ind = solver->getNumRows() + disj->common_changed_var.size() + disj->terms[term].changed_var.size();
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

  if (certificate_rank == num_nonzero_multipliers) {
    return RegularityStatus::REG;
  } else if (certificate_rank < num_nonzero_multipliers) {
    return RegularityStatus::IRREG_LESS;
  } else {
    return RegularityStatus::IRREG_MORE;
  }
} /* analyzeCertificateRegularity */

void analyzeCutRegularity(
    /// [in/out] Certificate of cut (in-version used to set MIP start), [term][Farkas multiplier]; per term, m + m_t + n indices correspond to rows (including globally-valid constraints) + disj term ineqs + cols
    std::vector<CutCertificate>& v,
    /// [out] Rank of submatrix associated to the certificate for each cut
    std::vector<int>& certificate_submx_rank,
    /// [out] Number of original (+ globally valid) constraints that have nonzero multipliers in the certificate
    std::vector<int>& num_nonzero_multipliers,
    /// [out] Regularity status (-1 if rank < num_nonzero_multipliers, 0 if equal, 1 if greater)
    std::vector<RegularityStatus>& regularity_status,
    /// [out] Number of iterations to converge to certificate
    std::vector<int>& num_iters,
    /// [in] Set of cuts that are to be analyzed for regularity
    const OsiCuts& cuts,
    /// [in] Disjunction from which cuts were generated
    const Disjunction* const disj,
    /// [in] Solver corresponding to instance for which cuts are valid
    const OsiSolverInterface* const solver,
    /// [in] Parameters for setting verbosity and logfile
    const StrengtheningParameters::Parameters& params) {
  if (cuts.sizeCuts() == 0) return;
  
   // Check that disjunction has not been lost
  if (disj == NULL) return;
  if (disj->terms.size() == 0) return;

  // Prepare vector sizes
  certificate_submx_rank.clear();
  certificate_submx_rank.resize(cuts.sizeCuts());

  num_nonzero_multipliers.clear();
  num_nonzero_multipliers.resize(cuts.sizeCuts());

  regularity_status.clear();
  regularity_status.resize(cuts.sizeCuts(), RegularityStatus::UNKNOWN);

  num_iters.clear();
  num_iters.resize(cuts.sizeCuts());

  // Prepare Atilde and btilde, which are common to all terms
  // This encompasses the original constraints
  // and all globally-valid inequalities identified as part of the disjunction
  CoinPackedMatrix Atilde;
  std::vector<double> btilde;
  prepareAtilde(Atilde, btilde, disj, solver, params.logfile);

  // Compute rank of Atilde
  // const int rank_atilde = computeRank(&Atilde, std::vector<int>(), std::vector<int>());

  // Prepare solver for computing certificate
  // TODO: Probably should just input directly to other solver if we are not using Cbc...
  OsiSolverInterface* liftingSolver = new SolverInterface;
  setLPSolverParameters(liftingSolver, params.get(StrengtheningParameters::VERBOSITY));
  genRCVMIPFromCut(liftingSolver, cuts.rowCutPtr(0), disj, solver, params);

  // Check that if Gurobi is selected, then USE_GUROBI is defined
  const bool use_gurobi = use_bb_option(params.get(StrengtheningParameters::intParam::BB_STRATEGY),
      StrengtheningParameters::BB_Strategy_Options::gurobi);
  const bool use_cbc = use_bb_option(params.get(StrengtheningParameters::intParam::BB_STRATEGY),
      StrengtheningParameters::BB_Strategy_Options::cbc);
      
#ifdef USE_GUROBI
  GRBModel* grbSolver;
#endif
  if (use_gurobi) {
#ifndef USE_GUROBI
    error_msg(errorstring, "analyzeCutRegularity: Gurobi is selected as the branch-and-bound solver, but USE_GUROBI is not defined.\n");
    writeErrorToLog(errorstring, params.logfile);
    throw std::logic_error(errorstring);
#else
    grbSolver = buildGRBModelFromOsi(liftingSolver, params.logfile); 
    setStrategyForBBTestGurobi(params, 0, *grbSolver);
#endif
  }

#ifdef USE_CBC
  CbcModel* cbcSolver;
#endif
  if (use_cbc) {
#ifndef USE_CBC
    error_msg(errorstring, "analyzeCutRegularity: Cbc is selected as the branch-and-bound solver, but USE_CBC is not defined.\n");
    writeErrorToLog(errorstring, params.logfile);
    throw std::logic_error(errorstring);
#else
    cbcSolver = new CbcModel(*liftingSolver);
    setIPSolverParameters(cbcSolver, params.get(StrengtheningParameters::VERBOSITY));
#endif
  }

  if (!use_gurobi && !use_cbc) {
    error_msg(errorstring, "analyzeCutRegularity: Implementation for solving RCVMIP is available only for Cbc and Gurobi.\n");
    writeErrorToLog(errorstring, params.logfile);
    throw std::logic_error(errorstring);
  }

  v.resize(cuts.sizeCuts());
  for (int cut_ind = 0; cut_ind < cuts.sizeCuts(); cut_ind++) {
    std::vector<double> solution; // RCVMIP solution
    bool reached_feasibility = false;

    // Set theta column for the current cut
    if (cut_ind > 0) {
      if (use_gurobi) {
#ifdef USE_GUROBI
        updateRCVMIPFromCut(grbSolver, cuts.rowCutPtr(cut_ind), disj, solver, params.logfile);
#endif
      } else {
        updateRCVMIPFromCut(cbcSolver, cuts.rowCutPtr(cut_ind), disj, solver, params.logfile);
      }
    }

    // Set MIP start for current cut
#ifdef USE_GUROBI
    if (use_gurobi) {
      setRCVMIPStart(grbSolver, solution, v[cut_ind], disj, solver);

      std::vector<int> delta_var_inds;
      const int certificate_rank = computeRankOfRCVMIPSolution(delta_var_inds, solution.data(), disj, solver, Atilde, params, cut_ind);

      if (certificate_rank < (int) delta_var_inds.size()) {
        addRankConstraint(grbSolver, delta_var_inds, certificate_rank, 0);
      } else {
        reached_feasibility = true;
      }
    }
#endif

    if (!reached_feasibility) {
      // Clear the certificate for the current cut, if one was provided
      if (v.size() > cut_ind) {
        v[cut_ind].clear();
      }

      // Solve the RCVMIP
      RCVMIPStatus return_code;
      if (use_gurobi) {
  #ifdef USE_GUROBI
        return_code = solveRCVMIP(grbSolver, solution, disj, solver, Atilde, params, cut_ind, num_iters[cut_ind]);
  #endif
      } else {
        return_code = solveRCVMIP(cbcSolver, solution, disj, solver, Atilde, params, cut_ind, num_iters[cut_ind]);
      }
      
      if (return_code == RCVMIPStatus::INFEASIBLE
          || return_code == RCVMIPStatus::UNBOUNDED
          || return_code == RCVMIPStatus::INF_OR_UNBD
          || return_code == RCVMIPStatus::ERROR) {
        error_msg(errorstring,
            "Solver terminates after %d iterations with status %d (%d: infeasible; %d: unbounded; %d: error) cut %d.\n",
            num_iters[cut_ind],
            static_cast<int>(return_code),
            static_cast<int>(RCVMIPStatus::INFEASIBLE),
            static_cast<int>(RCVMIPStatus::UNBOUNDED),
            static_cast<int>(RCVMIPStatus::ERROR),
            cut_ind);
        writeErrorToLog(errorstring, params.logfile);
        
        if (liftingSolver) { if (!use_cbc || !cbcSolver || !cbcSolver->modelOwnsSolver()) { delete liftingSolver; } }
  #ifdef USE_GUROBI
        if (use_gurobi && grbSolver) { delete grbSolver; }
  #endif
  #ifdef USE_CBC
        if (use_cbc && cbcSolver) { delete cbcSolver; }
  #endif
        
        throw std::logic_error(errorstring);
      } // exit out if infeasible or unbounded

      if (return_code == RCVMIPStatus::OPTIMAL_UNCONVERGED || return_code == RCVMIPStatus::RCVMIP_ITER_LIMIT) {
        regularity_status[cut_ind] = RegularityStatus::UNCONVERGED;
      }

      // Retrieve certificate from solution
      getCertificateFromRCVMIPSolution(v[cut_ind], solution, disj, solver, cut_ind, params.logfile);
    } // make sure initial solution is not already feasible

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

    // Compute rank of submatrix associated to the certificate
    const RegularityStatus curr_status = analyzeCertificateRegularity(certificate_submx_rank[cut_ind],
            num_nonzero_multipliers[cut_ind], v[cut_ind],
            disj, solver, Atilde, params);
    
    if (regularity_status[cut_ind] != RegularityStatus::UNCONVERGED) {
      regularity_status[cut_ind] = curr_status;
    }
  } // loop over cuts

  if (liftingSolver) { if (!use_cbc || !cbcSolver || !cbcSolver->modelOwnsSolver()) { delete liftingSolver; } }
#ifdef USE_GUROBI
  if (use_gurobi && grbSolver) { delete grbSolver; }
#endif
#ifdef USE_CBC
  if (use_cbc && cbcSolver) { delete cbcSolver; }
#endif
} /* analyzeCutRegularity */
