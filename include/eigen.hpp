/**
 * @file eigen.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-04
 */
#pragma once

#include <cstdio> // FILE
#include <vector>

// Project files
#include "CutCertificate.hpp" // TermCutCertificate and CutCertifcate

// COIN-OR files
class CoinPackedMatrix;
class OsiSolverInterface;

#ifdef USE_EIGEN
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
  const CoinPackedMatrix* const mat,
  const std::vector<int>& rows,
  const std::vector<int>& cols,
  const double MATRIX_EPS = 1e-10);
#endif // USE_EIGEN
