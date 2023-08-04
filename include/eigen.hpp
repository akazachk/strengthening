/**
 * @file eigen.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#pragma once

#include <vector>

class CoinPackedMatrix;
class OsiSolverInterface;

#ifdef USE_EIGEN
#include <Eigen/Sparse>

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
    std::vector<double>&v,
    const int num_elem,
    const int* const indices,
    const double* const elements,
    OsiSolverInterface* const solver,
    const std::vector<int>& rows,
    const std::vector<int>& cols,
    FILE* const logfile);

int computeRank(
  const CoinPackedMatrix* const mat,
  const std::vector<int>& rows);

#endif // USE_EIGEN