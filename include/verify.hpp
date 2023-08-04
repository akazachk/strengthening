/**
 * @file verify.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#pragma once

#include <vector>

class OsiSolverInterface;

/// @brief Take dot product with optimal basis and find resulting cut coefficients
void getCutFromCertificate(
    std::vector<double>& alpha, 
    const std::vector<double>& v, 
    const OsiSolverInterface* const solver);

/// @brief Count number of times (and by how much) the original cut coefficients are different from those calculated using the certificate \p v
void checkCut(
    int& num_errors, double& total_diff,
    const std::vector<double>& cut_coeff,
    const std::vector<double>& v,
    const OsiSolverInterface* const term_solver);