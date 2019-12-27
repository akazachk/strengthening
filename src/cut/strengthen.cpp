/**
 * @file strengthen.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#include "strengthen.hpp"

#include "Disjunction.hpp"
#include "Parameters.hpp" // includes SolverInterface
#include "SolverHelper.hpp" // isBasicVar
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

  const bool doAll = rows.size() + cols.size() == 0;
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
    /// [out] Farkas multipliers
    std::vector<double>& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] LP solver corresponding to disjunctive term
    const OsiSolverInterface* const solver) {
  v.clear();
  v.resize(solver->getNumCols() + solver->getNumRows(), 0.0);

  if (!solver->isProvenOptimal()) {
    return;
  }

  solver->enableFactorization();

  // Collect nonbasic variables
  // TODO Can we do this with getBInv calls instead?
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
    v[row] = x(tmp_ind);
    tmp_ind++;
  }
  for (const int& col : cols) {
    v[solver->getNumRows() + col] = x(tmp_ind);
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

  solver->disableFactorization();
} /* getCertificate */

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
    // [in] initial LP solver
    const OsiSolverInterface* const solver,
    // [in] disjunction
    const Disjunction* const disj) {
} /* getCertificateTrivial */

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

  for (int col = 0; col < solver->getNumCols(); col++) {
    const int start = mat->getVectorFirst(col);
    alpha[col] += dotProduct(mat->getVectorSize(col), 
        mat->getIndices() + start, mat->getElements() + start, v.data());
    //const double mult = isNonBasicUBVar(solver, col) ? -1.0 : 1.0;
    alpha[col] += v[solver->getNumRows() + col];
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

    printf("Calc\tAlpha\n");
    for (int col = 0; col < solver->getNumCols(); col++) {
      printf("%g\t%g\n", tmp(col), alpha[col]);
    }
  }*/ // DEBUG
} /* verifyCertificate */
