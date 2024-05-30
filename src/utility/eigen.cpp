/**
 * @file eigen.cpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#include "eigen.hpp"

// COIN-OR
#include <CoinPackedMatrix.hpp>
#include <OsiSolverInterface.hpp>

// Project
#include "utility.hpp" // isVal

#ifdef USE_EIGEN
// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

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
    /// [in] sparse matrix from COIN-OR, in row-major form
    const CoinPackedMatrix* const mat,
    /// [in] set of rows that we want to consider
    const std::vector<int>& rows,
    /// [in] set of variable bounds we want to explicitly use
    const std::vector<int>& cols) {
  if (mat->isColOrdered()) {
    throw std::logic_error("createEigenMatrix (sparse): Matrix must be row-ordered.");
  }

  const bool doAll = (rows.size() + cols.size()) == 0;
  const int num_cols = mat->getNumCols();
  const int num_selected_rows = doAll ? mat->getNumRows() + mat->getNumCols() : rows.size() + cols.size();
  int num_elem_removed = 0;
  //int num_rows_removed = 0;

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
      //num_rows_removed++;
    }
  }

  // Add explicit column lower bound rows
  // We *do not* ``complement'' the upper bounded variables because we would then have to
  // also complement ``alpha'' when we use it for the right-hand side;
  // these two negations cancel each other out
  for (const int& col : cols) {
    const double val = 1.0;
    if (batchInsert) tripletList.push_back(Eigen::Triplet<double>(tmp_ind, col, val));
    else M.insert(tmp_ind, col) = val;
    tmp_ind++;
  }

  if (doAll) {
    for (int col = 0; col < num_cols; col++) { // NB: num_cols is from mat, not the original solver, so there could be more columns if the last ones are empty, but this is unlikely and should not matter anyway
      const double val = 1.0;
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
  assert(M.nonZeros() == (long int) (mat->getNumElements() - num_elem_removed + (doAll ? num_cols : cols.size())));
#endif
} /* createEigenMatrix (sparse, CoinPackedMatrix) */

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
  if (mat->isColOrdered()) {
    throw std::logic_error("createEigenMatrix (sparse, solver): Matrix must be row-ordered.");
  }
#endif
  assert(mat->getNumCols() == solver->getNumCols());
  createEigenMatrix(M, mat, rows, cols);
} /* createEigenMatrix (sparse, solver) */

void createEigenMatrix(
    /// [out] dense matrix
    Eigen::MatrixXd& M,
    /// [in] sparse matrix from COIN-OR, in column-major form
    const CoinPackedMatrix* const mat, 
    /// [in] set of rows that we want to consider
    const std::vector<int>& rows) {
  error_msg(errorstring, "Creating dense eigen matrix not currently implemented.\n");
  exit(1);
} /* createEigenMatrix (dense) */

void solveLinearSystem(
    /// [out] solution to the linear system (if one was successfully found)
    Eigen::VectorXd& x,
    /// [in] A matrix
    const Eigen::SparseMatrix<double>& A,
    /// [in] right-hand side to the linear system
    const Eigen::VectorXd& b) {
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    // decomposition failed
    error_msg(errorstring, "solveLinearSystem: failed to create decomposition.\n");
    exit(1);
  }
  x = solver.solve(b);
  if (solver.info() != Eigen::Success) {
    // solving failed
    error_msg(errorstring, "solveLinearSystem: solving failed.\n");
    exit(1);
  }

#ifdef TRACE
  Eigen::VectorXd tmp = A * x;
  for (int i = 0; i < b.size(); i++) {
    if (!isVal(tmp(i), b(i))) {
      fprintf(stderr, "Calculated A_i . x = %f instead of correct value of b_i = %f.\n", tmp(i), b(i));
    }
  }
#endif
} /* solveLinearSystem */

void calculateCertificateEigen(
    /// [out] Farkas multipliers (vector of length m + m_t + n)
    TermCutCertificate& v,
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] LP solver corresponding to disjunctive term
    OsiSolverInterface* const term_solver,
    /// [in] set of rows that we want to consider
    const std::vector<int>& rows,
    /// [in] set of variable bounds we want to explicitly use
    const std::vector<int>& cols,
    /// [in] logfile for error printing
    FILE* const logfile) {
  Eigen::SparseMatrix<double,Eigen::RowMajor> A;
  createEigenMatrix(A, term_solver, rows, cols);
#ifdef TRACE
  assert(A.rows() == term_solver->getNumCols());
#endif

  // Set up right-hand side for term_solver (should be equal to alpha)
  Eigen::VectorXd b(term_solver->getNumCols());
  b.setZero();
  for (int tmp_ind = 0; tmp_ind < num_elem; tmp_ind++) {
    b(ind[tmp_ind]) = coeff[tmp_ind];
  }

  /*int tmp_ind = 0;
  for (const int& row : rows) {
    b(tmp_ind) = term_solver->getRightHandSide()[row];
    tmp_ind++;
  }
  for (const int& col : cols) {
    double val = 0.0;
    if (isVal(term_solver->getColSolution()[col], term_solver->getColLower()[col])) {
      val = term_solver->getColLower()[col];
    }
    else if (isVal(term_solver->getColSolution()[col], term_solver->getColUpper()[col])) {
      val = term_solver->getColUpper()[col];
    }
    else {
      error_msg(errorstring, 
          "Unable to identify which bound is met by nonbasic variable %d with value %.6g, and bounds [%.6g,%.6g]. May be a free nonbasic variable; check and then possibly add right-hand side equal to its current value instead of exiting.\n",
          col, term_solver->getColSolution()[col], term_solver->getColLower()[col], term_solver->getColUpper()[col]);
      exit(1);
    }
    b(tmp_ind) = val;
    tmp_ind++;
  }*/

  Eigen::VectorXd x(term_solver->getNumCols());
  solveLinearSystem(x, A.transpose(), b);

  /*{ // DEBUG
    Eigen::VectorXd tmp;
    tmp = A.transpose() * x;

    printf("Orig\tCalc\n");
    for (int col = 0; col < term_solver->getNumCols(); col++) {
      //printf("%g\t%g\n", term_solver->getColSolution()[col], x(col));
      printf("%g\t%g\n", b(col), tmp(col));
    }
  }*/ // DEBUG

  int tmp_ind = 0;
  for (const int& row : rows) {
    const double mult = 1.; //(term_solver->getRowSense()[row] == 'L') ? -1. : 1.;
    v[row] = mult * x(tmp_ind);
    tmp_ind++;
  }
  for (const int& col : cols) {
    const double mult = 1.; //isNonBasicUBVar(term_solver, col) ? -1. : 1.;
    v[term_solver->getNumRows() + col] = mult * x(tmp_ind);
    tmp_ind++;
  }

  /*{ // DEBUG
    Eigen::SparseMatrix<double,Eigen::RowMajor> fullA;
    createEigenMatrix(fullA, term_solver, std::vector<int>(), std::vector<int>());

    Eigen::VectorXd fullv(term_solver->getNumCols() + term_solver->getNumRows());
    for (int i = 0; i < (int) v.size(); i++) {
      fullv(i) = v[i];
    }

    Eigen::VectorXd tmp;
    tmp = fullA.transpose() * fullv;

    printf("Orig\tCalc\n");
    for (int col = 0; col < term_solver->getNumCols(); col++) {
      printf("%g\t%g\n", b(col), tmp(col));
    }
    exit(1);
  }*/ // DEBUG

  /*{ // DEBUG (print v)
    for (int i = 0; i < (int) v.size(); i++) {
      printf("(%d, %g)\n", i, v[i]);
    }
  }*/ // DEBUG
} /* calculateCertificateEigen */

int computeRank(
    const Eigen::MatrixXd& M,
    /// [in] Precision when taking submatrix rank
    const double MATRIX_EPS) {
  const bool USE_FULLPIVLU = false;
  if (USE_FULLPIVLU) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
    lu.setThreshold(MATRIX_EPS);
    return (int) lu.rank();
  } else {
    return M.colPivHouseholderQr().rank();
  }
} /* computeRank (Eigen::MatrixXd) */

int computeRank(
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& A,
    /// [in] Precision when taking submatrix rank
    const double MATRIX_EPS) {
  return computeRank(Eigen::MatrixXd(A), MATRIX_EPS);
} /* computeRank (Eigen::SparseMatrix) */

int computeRank(
    /// [in] CoinPackedMatrix from which rank will be computed
    const CoinPackedMatrix* const mat,
    /// [in] Rows of the matrix to use
    const std::vector<int>& rows,
    /// [in] Cols of the matrix to use
    const std::vector<int>& cols,
    /// [in] Precision when taking submatrix rank
    const double MATRIX_EPS) {
  if (mat->isColOrdered()) {
    throw std::logic_error("computeRank: Matrix must be row-ordered.");
  }

  // Create Eigen matrix from CoinPackedMatrix
  Eigen::SparseMatrix<double,Eigen::RowMajor> M;
  createEigenMatrix(M, mat, rows, cols);
  return computeRank(M, MATRIX_EPS);

  // const int num_rows = rows.size();
  // const int num_cols = mat->getNumCols();
  // std::vector<Eigen::Triplet<double> > tripletList;
  // tripletList.reserve(mat->getNumElements());

  // int tmp_ind = 0;
  // for (int i = 0; i < num_rows; i++) {
  //   const CoinShallowPackedVector row = mat->getVector(rows[i]);
  //   for (int j = 0; j < row.getNumElements(); j++) {
  //     tripletList.push_back(Eigen::Triplet<double>(tmp_ind, row.getIndices()[j], row.getElements()[j]));
  //     tmp_ind++;
  //   }
  // }
  // M.resize(num_rows, num_cols);
  // M.setFromTriplets(tripletList.begin(), tripletList.end());

  // Compute rank
  // Eigen::FullPivLU<Eigen::SparseMatrix<double> > lu(M);
  // return (int) lu.rank();
} /* computeRank */

#endif // USE_EIGEN
