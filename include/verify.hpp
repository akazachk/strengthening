/**
 * @file verify.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#pragma once

#include <vector>

// COIN-OR files
class OsiSolverInterface;

// Project files
#include "CutCertificate.hpp" // TermCutCertificate and CutCertifcate
class Disjunction;

/// @brief Take dot product of certificate \p v with \p term_solver rows to find resulting cut coefficients
void getCutFromCertificate(
    std::vector<double>& alpha, 
    const TermCutCertificate& v, 
    const OsiSolverInterface* const term_solver);

/// @brief Take dot product of certificate \p v with \p solver rows and \p disj globally-valid constraints to find resulting cut coefficients
void getCutFromCertificate(
    std::vector<double>& alpha,
    const TermCutCertificate& v,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver);

/// @brief Count number of times (and by how much) the original cut coefficients are different from those calculated using the certificate \p v
void checkCut(
    int& num_errors, double& total_diff,
    const std::vector<double>& cut_coeff,
    const TermCutCertificate& v,
    const OsiSolverInterface* const solver,
    const Disjunction* const disj);

/// @brief Verify that the certificate yields the same cut coefficients as the original cut (wrapper for #checkCut)
int checkCutHelper(
    const std::vector<double>& cut_coeff,
    const TermCutCertificate& v,
    const OsiSolverInterface* const solver,
    const Disjunction* const disj,    
    FILE* const logfile);