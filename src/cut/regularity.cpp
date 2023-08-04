/**
 * @file cglp.cpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#include "cglp.hpp"

#include <vector>

// COIN-OR
#include <CoinPackedMatrix.hpp>
#include <OsiRowCut.hpp>
#include <OsiSolverInterface.hpp>

// Project files
#include "Disjunction.hpp"
#include "utility.hpp" // isInfinity, error_msg

// #include "Parameters.hpp" // params.logfile
// using namespace StrengtheningParameters;

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
    /// [out] cut-generating linear program to be generated
    OsiSolverInterface* liftingSolver,
    /// [in] cut \alpha^T x \ge \beta for which we are seeking a certificate
    const OsiRowCut* const cut,
    /// [in] disjunction for which the cut is being generated
    const Disjunction* const disj,
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] logfile
    FILE* const logfile,
    /// [in] whether to use original RCVMILP or one with min sum delta objective and sum delta >= n constraint
    const bool use_min_sum_delta) {
  const int num_terms = disj->num_terms;
  const int num_rows = solver->getNumRows();
  const int num_common_rows = disj->common_changed_var.size() + disj->common_ineqs.size();
  const int num_cols = solver->getNumCols();
  const double* rhs = solver->getRightHandSide();

  // Let mprime be the number of rows of Atilde
  // This is the number of original constraints + number of globally-valid bound changes
  // We also add the number of lower and upper bounds to mprime
  int mprime_bounds = 0;
  for (int col = 0; col < num_cols; col++) {
      const double lb = solver->getColLower()[col];
      const double ub = solver->getColUpper()[col];

      if (!isInfinity(std::abs(lb))) {
          mprime_bounds++;
      }
      if (!isInfinity(std::abs(ub))) {
          mprime_bounds++;
      }
  }
  const int mprime = num_rows + num_common_rows + mprime_bounds;

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

  // Prepare block [-Atilde^T; btilde^T, -I_{m'}] repeated for all terms
  CoinPackedMatrix Atilde;
  std::vector<double> btilde;

  // Put Ax \ge b into Atilde and btilde
  Atilde.copyOf(*(solver->getMatrixByCol())); // copy A into mx, col-ordered, dimesions = [num_rows, num_cols]
  Atilde.reverseOrdering(); // make it row-ordered
  btilde.assign(rhs, rhs + num_rows); // copy b into btilde

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
  for (int col = 0; col < num_cols; col++) {
      const double val = solver->getColLower()[col];
      if (!isInfinity(std::abs(val))) {
          Atilde.appendRow(1, &col, &one);
          btilde.push_back(val);
      }
  }

  // Add rows corresponding to variable upper bounds to Atilde
  // These will be stored as -x_j \ge -ub_j
  for (int col = 0; col < num_cols; col++) {
      const double val = solver->getColUpper()[col];
      if (!isInfinity(std::abs(val))) {
          Atilde.appendRow(1, &col, &negone);
          btilde.push_back(negone * val);
      }
  }

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
  for (int i = 0; i < mprime; i++) {
      common_mx.appendRow(1, &i, &negone);
  }

  // Verify size of common_mx, and if incorrect, throw an error giving the sizes
  if (common_mx.getNumRows() != num_cols + 1 + mprime) {
      error_msg(errorstring,
          "genRCVMILPFromCut:: Num rows (mx, actual, predicted) = (Atilde, %d, %d); (common_mx, %d, %d).\n",
          Atilde.getNumRows(), mprime, common_mx.getNumRows(), num_cols + 1 + mprime);
      writeErrorToLog(errorstring, logfile);
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
    // These are variables 0, 1, ..., m'
    term_mx.rightAppendPackedMatrix(theta_delta);

    // Add zero_mx for the previous terms
    for (int t = 0; t < term_ind; t++) {
      term_mx.rightAppendPackedMatrix(zero_mx);
    }

    // Add common_mx, which is row-ordered and has dimensions [n + 1 + m', m']
    term_mx.majorAppendOrthoOrdered(common_mx);
    
    // Add zero_mx for the remaining terms
    for (int t = term_ind + 1; t < num_terms; t++) {
      term_mx.rightAppendPackedMatrix(zero_mx);
    }

    // Add term-specific zero_mx for the previous terms
    for (int t = 0; t < term_ind; t++) {
      const DisjunctiveTerm& curr_term = disj->terms[t];
      const int mt = curr_term.changed_var.size() + curr_term.ineqs.size();
      CoinPackedMatrix term_zero_mx;
      term_zero_mx.setDimensions(common_mx.getNumRows(), mt);
      term_mx.rightAppendPackedMatrix(term_zero_mx);
    }

    // Create matrix for the negative of the disjunctive term constraints
    // We assume D^t x \ge D^t_0 is stored as >= constraints
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
  const char* row_sense = solver->getRowSense();
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

  // Write to file
  const bool write_lp = true;
  if (write_lp) {
    printf("\n## Saving CGLP from cut to file: %s\n", "cglp_from_cut.lp");
    printf("The cut is:\n");
    cut->print();

    std::string lp_filename = "cglp_from_cut";
    liftingSolver->writeLp(lp_filename.c_str());
  }
} /* genRCVMILPFromCut */