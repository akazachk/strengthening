/**
 * @file strengthen.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#include "strengthen.hpp"

#include <OsiRowCut.hpp> // for setting a disjunctive term

#include "CutHelper.hpp" // for setObjective and addToObjective
#include "Disjunction.hpp"
//#include "Parameters.hpp"
#include "SolverHelper.hpp" // isBasicVar, checkSolverOptimality, isNonBasicUBVar
#include "utility.hpp"

#ifdef USE_EIGEN
// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;

void insertRow(
    /// [out] sparse matrix, in row-major form
    Eigen::SparseMatrix<double,Eigen::RowMajor>& M,
    /// [in] sparse matrix from COIN-OR, in row-major form
    const CoinPackedMatrix* const mat, 
    /// [in] row we are inserting
    const int row,
    /// [in] row in which to insert this entry
    const int tmp_row) {
  const int* cols = mat->getIndices();
  const double* elem = mat->getElements();

  const int start = mat->getVectorFirst(row);
  const int end = mat->getVectorLast(row);
  for (int ind = start; ind < end; ind++) {
    const double j = cols[ind];
    const double el = elem[ind];
    M.insert(tmp_row,j) = el;
  }
} /* insertRow */

void prepareRow(
    /// [out] Eigen::Triplet<double> that we are updating
    std::vector<Eigen::Triplet<double> >& tripletList,
    /// [in] sparse matrix from COIN-OR, in row-major form
    const CoinPackedMatrix* const mat, 
    /// [in] row we are inserting
    const int row,
    /// [in] row in which to insert this value
    const int tmp_row) {
  const int* cols = mat->getIndices();
  const double* elem = mat->getElements();

  const int start = mat->getVectorFirst(row);
  const int end = mat->getVectorLast(row);
  for (int ind = start; ind < end; ind++) {
    const double j = cols[ind];
    const double el = elem[ind];
    tripletList.push_back(Eigen::Triplet<double>(tmp_row,j,el));
  }
} /* prepareRow */

void createEigenMatrix(
    /// [out] sparse matrix, in row-major form
    Eigen::SparseMatrix<double,Eigen::RowMajor>& M, 
    /// [in] COIN-OR solver
    const OsiSolverInterface* const solver, 
    /// [in] set of rows that we want to consider
    const std::vector<int>& rows,
    /// [in] set of variable bounds we want to explicitly use
    const std::vector<int>& cols) {
  // Get sparse matrix from COIN-OR, in row-major form
  const CoinPackedMatrix* const mat = solver->getMatrixByRow();
#ifdef TRACE
  assert(!mat->isColOrdered());
#endif

  const bool doAll = (rows.size() + cols.size()) == 0;
  const int num_cols = mat->getNumCols();
  const int num_selected_rows = doAll ? mat->getNumRows() + mat->getNumCols() : rows.size() + cols.size();
  int num_elem_removed = 0, num_rows_removed = 0;

  // Prepare sparse matrix
  M.resize(num_selected_rows, num_cols);
  M.reserve(mat->getNumElements() + cols.size());

  const bool batchInsert = true;
  std::vector<Eigen::Triplet<double> > tripletList;
  if (batchInsert) {
    tripletList.reserve(mat->getNumElements() + cols.size());
  }
  unsigned tmp_ind = 0;
  for (int row = 0; row < mat->getNumRows(); row++) {
    if (doAll || (tmp_ind < rows.size() && rows[tmp_ind] == row)) {
      if (batchInsert) prepareRow(tripletList, mat, row, tmp_ind);
      else insertRow(M, mat, row, tmp_ind);
      tmp_ind++;
    } else {
      num_elem_removed += mat->getVectorSize(row);
      num_rows_removed++;
    }
  }
  // Add explicit column lower bound rows
  // We *do not* ``complement'' the upper bounded variables because we would then have to
  // also complement ``alpha'' when we use it for the right-hand side;
  // these two negations cancel each other out
  for (const int& col : cols) {
    const double val = 1.0; // isNonBasicUBVar(solver, col) ? -1.0 : 1.0;
    if (batchInsert) tripletList.push_back(Eigen::Triplet<double>(tmp_ind, col, val));
    else M.insert(tmp_ind, col) = val;
    tmp_ind++;
  }
  if (doAll) {
    for (int col = 0; col < solver->getNumCols(); col++) {
      const double val = 1.0; // isNonBasicUBVar(solver, col) ? -1.0 : 1.0;
      if (batchInsert) tripletList.push_back(Eigen::Triplet<double>(tmp_ind, col, val));
      else M.insert(tmp_ind, col) = val;
      tmp_ind++;
    }
  }
  if (batchInsert) {
    M.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  M.makeCompressed();

#ifdef TRACE
  // Check matrix
  assert(M.rows() == num_selected_rows);
  assert(M.cols() == num_cols);
  assert(M.nonZeros() == mat->getNumElements() - num_elem_removed + (doAll ? num_cols : cols.size()));
#endif
} /* createEigenMatrix (sparse) */

void createEigenMatrix(
    /// [out] dense matrix
    MatrixXd& M,
    /// [in] sparse matrix from COIN-OR, in column-major form
    const CoinPackedMatrix* const mat, 
    /// [in] set of rows that we want to consider
    const std::vector<int>& rows) {
  error_msg(errorstring, "Creating dense eigen matrix not currently implemented.\n");
  exit(1);
} /* createEigenMatrix (dense) */

void solveLinearSystem(
    /// [out] solution to the linear system (if one was successfully found)
    VectorXd& x,
    /// [in] A matrix
    const SparseMatrix<double>& A,
    /// [in] right-hand side to the linear system
    const VectorXd& b) {
  SparseLU<SparseMatrix<double> > solver;
  solver.compute(A);
  if (solver.info()!=Success) {
    // decomposition failed
    error_msg(errorstring, "solveLinearSystem: failed to create decomposition.\n");
    exit(1);
  }
  x = solver.solve(b);
  if (solver.info()!=Success) {
    // solving failed
    error_msg(errorstring, "solveLinearSystem: solving failed.\n");
    exit(1);
  }
} /* solveLinearSystem */
#endif // USE_EIGEN

/// Given problem is (with A an [m x n] matrix)
///   Ax - s = b
///   x \ge 0
///   s \ge 0
/// Suppose the basis is B (m columns); note that some of these might be slack variables
/// Denote by N the remaining n columns (the cobasis)
///
/// If solver is proven optimal, and the cut is valid for the basis cone,
/// then it suffices to take the inverse of the optimal basis
/// In particular, we want \alpha^T = v^T A
/// 
/// Note that \alpha is n-dimensional, v is m-dimensional
/// The components of v corresponding to rows in which the slack is basic will be 0
/// We need to get access to the inverse of the nxn submatrix of A corresponding to nonbasic slacks
///
/// * Invertible case
/// If A is invertible,
/// then we can take v^T = \alpha^T A^{-1}
/// or, equivalently, v = (A^{-1})^T \alpha
/// This means that v_j can be calcuated as the dot product of the j-th column of A^{-1} with \alpha
///
/// * General case
/// If A is not invertible, then we need only consider the relaxed corner polyhedron
/// That is, we only need the basis inverse applied to the matrix
void getCertificate(
    /// [out] Farkas multipliers (vector of length m + n)
    std::vector<double>& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] LP solver corresponding to disjunctive term
    OsiSolverInterface* const solver) {
  v.clear();
  v.resize(solver->getNumCols() + solver->getNumRows(), 0.0);

  // Get dense cut coefficients so that we can set the objective vector
  std::vector<double> cut_coeff(solver->getNumCols(), 0.0);
  for (int i = 0; i < num_elem; i++) {
    cut_coeff[ind[i]] = coeff[i];
  }
  solver->setObjective(cut_coeff.data());
  solver->resolve();

  if (!solver->isProvenOptimal()) {
    return;
  }

  solver->enableFactorization();

  // Collect nonbasic variables
  // TODO Can / should we do this with getBInv calls instead?
  std::vector<int> rows, cols;
  //std::vector<int> NBVarIndex;
  //NBVarIndex.reserve(solver->getNumCols());
  for (int var = 0; var < solver->getNumCols() + solver->getNumRows(); var++) {
    if (!isBasicVar(solver, var)) {
      //NBVarIndex.push_back(var);
      if (var < solver->getNumCols()) {
        cols.push_back(var);
      } else {
        rows.push_back(var - solver->getNumCols());
      }
    }
  }

#ifdef USE_EIGEN
  Eigen::SparseMatrix<double,Eigen::RowMajor> A;
  createEigenMatrix(A, solver, rows, cols);
#ifdef TRACE
  assert(A.rows() == solver->getNumCols());
#endif

  // Now set up right-hand side for solver (should be equal to alpha)
  Eigen::VectorXd b(solver->getNumCols());
  for (int tmp_ind = 0; tmp_ind < num_elem; tmp_ind++) {
    b(ind[tmp_ind]) = coeff[tmp_ind];
  }

  /*int tmp_ind = 0;
  for (const int& row : rows) {
    b(tmp_ind) = solver->getRightHandSide()[row];
    tmp_ind++;
  }
  for (const int& col : cols) {
    double val = 0.0;
    if (isVal(solver->getColSolution()[col], solver->getColLower()[col])) {
      val = solver->getColLower()[col];
    }
    else if (isVal(solver->getColSolution()[col], solver->getColUpper()[col])) {
      val = solver->getColUpper()[col];
    }
    else {
      error_msg(errorstring, 
          "Unable to identify which bound is met by nonbasic variable %d with value %.6g, and bounds [%.6g,%.6g]. May be a free nonbasic variable; check and then possibly add right-hand side equal to its current value instead of exiting.\n",
          col, solver->getColSolution()[col], solver->getColLower()[col], solver->getColUpper()[col]);
      exit(1);
    }
    b(tmp_ind) = val;
    tmp_ind++;
  }*/

  Eigen::VectorXd x(solver->getNumCols());
  solveLinearSystem(x, A.transpose(), b);

  /*{ // DEBUG
    Eigen::VectorXd tmp;
    tmp = A.transpose() * x;

    printf("Orig\tCalc\n");
    for (int col = 0; col < solver->getNumCols(); col++) {
      //printf("%g\t%g\n", solver->getColSolution()[col], x(col));
      printf("%g\t%g\n", b(col), tmp(col));
    }
  }*/ // DEBUG

  int tmp_ind = 0;
  for (const int& row : rows) {
    const double mult = isNonBasicUBSlack(solver, row) ? -1. : 1.;
    v[row] = mult * x(tmp_ind);
    tmp_ind++;
  }
  for (const int& col : cols) {
    const double mult = isNonBasicUBVar(solver, col) ? -1. : 1.;
    v[solver->getNumRows() + col] = mult * x(tmp_ind);
    tmp_ind++;
  }

  /*{ // DEBUG
    Eigen::SparseMatrix<double,Eigen::RowMajor> fullA;
    createEigenMatrix(fullA, solver, std::vector<int>(), std::vector<int>());

    Eigen::VectorXd fullv(solver->getNumCols() + solver->getNumRows());
    for (int i = 0; i < (int) v.size(); i++) {
      fullv(i) = v[i];
    }

    Eigen::VectorXd tmp;
    tmp = fullA.transpose() * fullv;

    printf("Orig\tCalc\n");
    for (int col = 0; col < solver->getNumCols(); col++) {
      printf("%g\t%g\n", b(col), tmp(col));
    }
    exit(1);
  }*/ // DEBUG

  /*{ // DEBUG (print v)
    for (int i = 0; i < (int) v.size(); i++) {
      printf("(%d, %g)\n", i, v[i]);
    }
  }*/ // DEBUG
#endif // USE_EIGEN

#ifdef TRACE
  // Obtain the cut that the certificate yields (should be the same as the original cut)
  std::vector<double> new_coeff(solver->getNumCols());
  verifyCertificate(new_coeff, v, solver);

  int num_errors = 0;
  double total_diff = 0.;
  for (int i = 0; i < solver->getNumCols(); i++) {
    const double diff = cut_coeff[i] - new_coeff[i];
    if (greaterThanVal(std::abs(diff), 0.0)) {
      fprintf(stderr, "%d: cut: %.6f\tcalc: %.6f\tdiff: %g\n", i, cut_coeff[i], new_coeff[i], diff);
      num_errors++;
      total_diff += std::abs(diff);
    }
  }
  if (num_errors > 0) printf("Number of differences between true and calculated cuts: %d. Total difference: %g.\n", num_errors, total_diff);
  //else printf("No errors found.\n");
#endif

  solver->disableFactorization();
} /* getCertificate */

void getCertificateForTerm(
    /// [out] Farkas multipliers
    std::vector<double>& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] disjunctive term specification
    const DisjunctiveTerm* const term,
    /// [in] tolerance for reconstructing disjunctive term solution
    const double DIFFEPS,
    /// [in] logfile for error printing
    FILE* logfile) {
  OsiSolverInterface* termSolver = solver->clone();
  const int curr_num_changed_bounds = term->changed_var.size();
  std::vector < std::vector<int> > commonTermIndices(curr_num_changed_bounds);
  std::vector < std::vector<double> > commonTermCoeff(curr_num_changed_bounds);
  std::vector<double> commonTermRHS(curr_num_changed_bounds);
  for (int i = 0; i < curr_num_changed_bounds; i++) {
    const int col = term->changed_var[i];
    const double coeff = (term->changed_bound[i] <= 0) ? 1. : -1.;
    const double val = term->changed_value[i];
    commonTermIndices[i].resize(1, col);
    commonTermCoeff[i].resize(1, coeff);
    commonTermRHS[i] = coeff * val;
    if (term->changed_bound[i] <= 0) {
      termSolver->setColLower(col, val);
    } else {
      termSolver->setColUpper(col, val);
    }
  }

  const int curr_num_added_ineqs = term->ineqs.size();
  for (int i = 0; i < curr_num_added_ineqs; i++) {
    const OsiRowCut* currCut = &term->ineqs[i];
    termSolver->applyRowCuts(1, currCut); // hopefully this works
  }

  // Set the warm start
  if (term->basis && !(termSolver->setWarmStart(term->basis))) {
    error_msg(errorstring,
        "Warm start information not accepted for term.\n");
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  // Resolve and check the objective matches
#ifdef TRACE
  printf("\n## Solving for term ##\n");
#endif
  termSolver->resolve();
  const bool calcAndFeasTerm = checkSolverOptimality(termSolver, true);

  // If something went wrong
  if (!calcAndFeasTerm) {
    printf("\n## Term is not proven optimal. Exiting from this term. ##\n");
    delete termSolver;
    return;
  }

  // Sometimes we run into a few issues getting the ``right'' value
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    termSolver->resolve();
  }
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    double ratio = termSolver->getObjValue() / term->obj;
    if (ratio < 1.) {
      ratio = 1. / ratio;
    }
    // Allow it to be up to 3% off without causing an error
    if (greaterThanVal(ratio, 1.03)) {
      error_msg(errorstring,
          "Objective at disjunctive term is incorrect. Before, it was %s, now it is %s.\n",
          stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
      writeErrorToLog(errorstring, logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Objective at disjunctive term is incorrect. Before, it was %s, now it is %s.\n",
          stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
    }
#ifdef TRACE
    std::string commonName;
    Disjunction::setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
        commonTermCoeff, commonTermRHS, false);
    printf("Bounds changed: %s.\n", commonName.c_str());
#endif
  } // check that objective value matches

  getCertificate(v, num_elem, ind, coeff, termSolver);
  if (termSolver) { delete termSolver; }
} /* getCertificateForTerm */

/// Use Theorem B.5 of Kazachkov 2018 dissertation to get certificate with trivial normalization
void getCertificateTrivial(
    /// [out] Farkas multipliers
    std::vector<double>& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] disjunction
    const Disjunction* const disj,
    // [in] initial LP solver
    const OsiSolverInterface* const solver) {
  // TODO
} /* getCertificateTrivial */

/// First m+m_t rows of v correspond to A;D^t; the next n are bounds on the variables
void verifyCertificate(
    /// [out] calculated cut coefficients
    std::vector<double>& alpha, 
    /// [in] Farkas multipliers
    const std::vector<double>& v, 
    /// [in] LP solver corresponding to disjunctive term
    const OsiSolverInterface* const solver) {
  alpha.clear();
  alpha.resize(solver->getNumCols(), 0.0);

  const CoinPackedMatrix* mat = solver->getMatrixByCol();

  std::vector<double> new_v(v.begin(), v.end());
  for (double& val : new_v) {
    val = std::abs(val);
  }
  for (int col = 0; col < solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    alpha[col] += dotProduct(mat->getVectorSize(col), 
        mat->getIndices() + start, mat->getElements() + start, v.data());
    alpha[col] += v[solver->getNumRows() + col];
  }

  /*
  for (int col = 0; col < solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    for (int el_ind = 0; el_ind < mat->getVectorSize(col); el_ind++) {
      const int row = mat->getIndices()[start + el_ind];
      const double mult = isNonBasicUBSlack(solver, row) ? -1.0 : 1.0;
      alpha[col] += mult * v[row];
    }
    const double mult = isNonBasicUBVar(solver, col) ? -1.0 : 1.0;
    alpha[col] += mult * v[solver->getNumRows() + col];
  }
  */
} /* verifyCertificate */

/**
 * @brief Attempt to strengthen coefficients of given cut
 *
 * Required:
 * (1) original cut
 * (2) disjunction for which the cut is valid
 * (3) Farkas certificate for the cut for this disjunction, 
 *      where the first m rows correspond to the original constraint matrix, 
 *      the next m_t rows correspond to the multipliers on the disjunctive term inequalities
 *      the last n rows are bounds on the variables, 
 * (4) globally valid lower bounds on the disjunctive term inequalties
 *
 * With those, the strengthening requires finding values m_1,...,m_T
 *  str_coeff[k] = coeff[k] + max_t { -u^t_k + u^t_0 (D^t_0 - \ell^t) m_t }
 *
 * TODO If k is nonbasic at a nonzero bound or is at its upper bound, what do we do?
 */
void strengthenCut(
    /// [out] strengthened cut coefficients
    std::vector<double>& str_coeff,
    /// [out] strengthened cut rhs
    double& str_rhs,
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    /// [in] original cut rhs
    const double rhs,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const std::vector<std::vector<double> >& v, 
    /// [in] original solver (used to get globally-valid lower bounds for the disjunctive terms)
    const OsiSolverInterface* const solver) {
  str_coeff.clear();
  str_coeff.resize(solver->getNumCols(), 0.0);
  if (!disj) return;

  // Set up original cut coeff and rhs
  str_rhs = rhs;
  for (int i = 0; i < num_elem; i++) {
    str_coeff[ind[i]] = coeff[i];
  }

  // Get globally lower bound for each of the disjunctive terms
  // and use this to set up term-by-term coefficients for strengthening
  // TODO get this to work with generic inequality description of disjunction, not just the bound changes version
  std::vector<std::vector<double> > disj_lb_diff(disj->num_terms);
  std::vector<double> lb_term(disj->num_terms, 0.0); // u^t_0 (D^t_0 - \ell^t)
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const DisjunctiveTerm& term = disj->terms[term_ind];
    const int num_disj_ineqs = (int) term.changed_var.size();
    // Resize current term vector
    disj_lb_diff[term_ind].resize(num_disj_ineqs, 0.0);

    for (int bound_ind = 0; bound_ind < num_disj_ineqs; bound_ind++) {
      const int var = term.changed_var[bound_ind];
      const double mult = (term.changed_bound[bound_ind] <= 0) ? 1.0 : -1.0;
      const double bd = (term.changed_bound[bound_ind] <= 0) ? solver->getColLower()[var] : solver->getColUpper()[var];
      double lb = mult * bd;
      if (isInfinity(std::abs(lb))) {
        {
          // No bound readily available; try to optimize to find one
          OsiSolverInterface* tmpSolver = solver->clone();
          addToObjectiveFromPackedVector(tmpSolver, NULL, true);
          setConstantObjectiveFromPackedVector(tmpSolver, mult, 1, &var);
          tmpSolver->resolve();
          checkSolverOptimality(tmpSolver, false);
          if (tmpSolver->isProvenOptimal()) {
            lb = tmpSolver->getObjValue();
          }

          if (tmpSolver) { delete tmpSolver; }
        }
        if (isInfinity(std::abs(lb))) {
          fprintf(stderr,
              "*** ERROR: Cannot strengthen cut using this disjunction because variable %d is missing its %s bound.\n",
              var, term.changed_bound[bound_ind] <= 0 ? "lower" : "upper");
          exit(1);
        }
      } // check if lb is infinite

      disj_lb_diff[term_ind][bound_ind] = mult * term.changed_value[bound_ind] - lb;
#ifdef TRACE
      assert(disj_lb_diff[term_ind][bound_ind] > -1e-3);
#endif
      lb_term[term_ind] += v[term_ind][solver->getNumRows() + bound_ind] * disj_lb_diff[term_ind][bound_ind];
    } // loop over bounds
  } // loop over terms

  // Now try to strengthen the coefficients on the integer-restricted variables
  for (int col = 0; col < solver->getNumCols(); col++) {
    if (!solver->isInteger(col)) continue;
    strengthenCutCoefficient(str_coeff[col], str_rhs, col, coeff[col], disj, lb_term, v, solver);
  }
} /* strengthenCut */

/// Returns new cut coefficient (we need to minimize over monoids m)
void strengthenCutCoefficient(
    /// [out] strengthened cut coeff
    double& str_coeff,
    /// [out] strengthened cut rhs
    double& str_rhs,
    /// [in] variable for which the coefficient is being changed
    const int var,
    /// [in] original coeff
    const double coeff,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] u^t_0 (D^t_0 - \ell^t) for all t
    const std::vector<double>& lb_term,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const std::vector<std::vector<double> >& v, 
    /// [in] original solver
    const OsiSolverInterface* const solver) {
  const int num_terms = disj->num_terms;
  if (num_terms > 2) {
    fprintf(stderr, "Currently strengthening is limited to 2-term disjunctions.\n");
    exit(1);
  }

  // Check if complementing is necessary
  double mult = 1.;
  if (v[0][solver->getNumRows() + disj->terms[0].changed_var.size() + var] < 0) {
    mult = -1;
  }
  if (v[1][solver->getNumRows() + disj->terms[1].changed_var.size() + var] < 0) {
    mult = -1;
  }
  str_coeff = mult * coeff;

  std::vector<int> m1(num_terms, 0.0);
  std::vector<int> m2(num_terms, 0.0);
  std::vector<int> m3(num_terms, 0.0);
  if (num_terms == 2) {
    m2[0] = 1.0;
    m2[1] = -1.0;
    m3[0] = -1.0;
    m3[1] = 1.0;
  }
  std::vector<std::vector<int> > m_options = { m1, m2, m3 };
  double min_max_term_val = std::numeric_limits<double>::max();
  for (const auto& m : m_options) {
    double max_term_val = std::numeric_limits<double>::lowest(); 
    for (int term_ind = 0; term_ind < num_terms; term_ind++) {
      const DisjunctiveTerm& term = disj->terms[term_ind];
      const int num_disj_ineqs = (int) term.changed_var.size();
      const double utk = v[term_ind][solver->getNumRows() + num_disj_ineqs + var];
      const double curr_val = -utk + lb_term[term_ind] * m[term_ind];
      if (curr_val > max_term_val) {
        max_term_val = curr_val;
      }
    }
    if (max_term_val < min_max_term_val) {
      min_max_term_val = max_term_val;
    }
  } // loop over monoid options
  str_coeff += min_max_term_val;

  if (mult < 0) {
    str_coeff *= -1;
    str_rhs = str_rhs + (str_coeff - coeff) * solver->getColUpper()[var];
  }
} /* strengthenCutCoefficient */
