/**
 * @file eigen.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#pragma once

#include <vector>

// Project files
#include "CutCertificate.hpp" // TermCutCertificate and CutCertifcate

// COIN-OR files
class CoinPackedMatrix;
class OsiSolverInterface;

#ifdef USE_EIGEN
#include <Eigen/Sparse>
#include <Eigen/Dense>

void insertRow(
    Eigen::SparseMatrix<double,Eigen::RowMajor>& M,
    const CoinPackedMatrix* const mat,
    const int row,
    const int tmp_row);

void prepareRow(
    std::vector<Eigen::Triplet<double> >& tripletList,
    const CoinPackedMatrix* const mat,
    const int row,
    const int tmp_row);

void createEigenMatrix(
    Eigen::SparseMatrix<double,Eigen::RowMajor>& M,
    const CoinPackedMatrix* const mat,
    const std::vector<int>& rows,
    const std::vector<int>& cols);

void createEigenMatrix(
    Eigen::SparseMatrix<double,Eigen::RowMajor>& M,
    const OsiSolverInterface* const solver,
    const std::vector<int>& rows,
    const std::vector<int>& cols);

void createEigenMatrix(
    Eigen::MatrixXd& M,
    const CoinPackedMatrix* const mat,
    const std::vector<int>& rows);

void solveLinearSystem(
    Eigen::VectorXd& x,
    const Eigen::SparseMatrix<double>& A,
    const Eigen::VectorXd& b);

void calculateCertificateEigen(
    TermCutCertificate& v,
    const int num_elem,
    const int* const indices,
    const double* const elements,
    OsiSolverInterface* const solver,
    const std::vector<int>& rows,
    const std::vector<int>& cols,
    FILE* const logfile);

int computeRank(
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& M,
  const double MATRIX_EPS = 1e-10);

int computeRank(
  const Eigen::MatrixXd& M,
  const double MATRIX_EPS = 1e-10);

int computeRank(
  const CoinPackedMatrix* const mat,
  const std::vector<int>& rows,
  const std::vector<int>& cols,
  const double MATRIX_EPS = 1e-10);

#endif // USE_EIGEN