/**
 * @file strengthen.hpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#pragma once

#include <vector>
#include <cstdio>

class OsiSolverInterface;
class CoinPackedMatrix;

class Disjunction;
class DisjunctiveTerm;

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
#endif

/**
 * @brief Return Farkas certificate given a valid cut (in sparse form) 
 * and a linear program describing the feasible region of a polyhedron
 */
void getCertificate(std::vector<double>& v, const int num_elem, const int* const ind, 
    const double* const coeff, OsiSolverInterface* const solver);

void getCertificateForTerm(
    std::vector<double>& v, 
    const int num_elem, 
    const int* const ind, 
    const double* const coeff,
    const OsiSolverInterface* const si,
    const DisjunctiveTerm* const term,
    const double DIFFEPS,
    FILE* logfile);

void getCertificateTrivial(std::vector<double>& v, const int num_elem, const int* const ind, 
    const double* const coeff, const OsiSolverInterface* const solver, const Disjunction* const disj);

/**
 * @brief Return Farkas certificate for a given valid cut (in sparse form) and disjunction
 */
void getCertificateForDisjunction();

/**
 * @brief Take dot product with optimal basis and find resulting cut coefficients
 */
void verifyCertificate(std::vector<double>& alpha, const std::vector<double>& v, 
    const OsiSolverInterface* const solver);

/**
 * @brief Attempt to strengthen coefficients of given cut
 */
void strengthenCut(
    std::vector<double>& str_coeff,
    double& str_rhs,
    const int num_elem, 
    const int* const ind, 
    const double* const coeff,
    const double rhs,
    const Disjunction* const disj,
    const std::vector<std::vector<double> >& v, 
    const OsiSolverInterface* const solver);

/**
 * @brief Returns new cut coefficient
 */
double strengthenCutCoefficient(
    const int var,
    const double coeff,
    const Disjunction* const disj,
    const std::vector<double>& lb_term,
    const std::vector<std::vector<double> >& v, 
    const OsiSolverInterface* const solver);
