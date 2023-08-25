/**
 * @file regularity.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#pragma once

#include <cstdio> // FILE
#include <vector>
#include <string>

// Project files
#include "CutCertificate.hpp"
class Disjunction;
namespace StrengtheningParameters {
  struct Parameters;
}

// COIN-OR files
class CoinPackedMatrix;
class OsiSolverInterface;
class OsiRowCut;
class OsiCuts;

/// @brief Regularity status enumeration to track whether we have identified that a cut is regular or not
enum class RegularityStatus {
  IRREG_LESS = -1,          ///< All certificates have n > nnz > submx_rank
  REG = 0,                  ///< There exists a certificate with submx_rank == nnz
  IRREG_MORE = 1,           ///< All certificates have nnz > n
  // TENTATIVE_IRREG_LESS = 2, ///< We have not yet verified that the cut is irregular, but we have a certificate that it is irregular< (i.e., cert submatrix rank < cert num nonzero multipliers < num)
  TENTATIVE_IRREG_MORE = 3, ///< This could be either a certificate that has more than n nnz, or RCVMIP with optimal value 0 but in the presence of prior/extra rank constraints only
  UNCONVERGED = 4,          ///< RCVMIP did not converge (or we ran out of time before trying for that cut)
  NUMERICALLY_UNSTABLE = 5, ///< RCVMIP leads to nonzero, but very small, theta value
  UNKNOWN = 6
};

/// @brief Return string with regularity status name
const std::string getRegularityStatusName(const RegularityStatus& status);

/// @brief Calculate number of finite lower and upper bounds
int calculateNumFiniteBounds(
  const OsiSolverInterface* const solver,
  int* const num_lb = NULL, int* const num_ub = NULL);

/// @brief Calculate the number of rows of Atilde
int calculateNumRowsAtilde(
    const Disjunction* const disj,
    const OsiSolverInterface* const solver);

/// @brief Create CoinPackedMatrix containing original constraints, globally-valid inequalities, and variable bounds
void prepareAtilde(
    CoinPackedMatrix& Atilde,
    std::vector<double>& btilde,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    FILE* logfile);

/// @brief Generate a CGLP from a cut
void genRCVMIPFromCut(
    OsiSolverInterface* const liftingSolver,
    const OsiRowCut* const cut,
    const Disjunction* const disj, 
    const OsiSolverInterface* const solver,
    const StrengtheningParameters::Parameters& params,
    const bool use_min_sum_delta = false);

/// @brief Given an existing RCVMIP, update the coefficients related to the cut
void updateRCVMIPFromCut(
    OsiSolverInterface* const liftingSolver,
    const OsiRowCut* const cut,
    const Disjunction* const disj, 
    const OsiSolverInterface* const solver);

/// @brief Obtain solution to RCVMIP
int solveRCVMIP(
    OsiSolverInterface* const liftingSolver,
    std::vector<double>& solution,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const CoinPackedMatrix& Atilde,
    const StrengtheningParameters::Parameters& params,
    const int cut_ind);

/// @brief Given a solution to the RCVMIP, extract the Farkas multipliers
void getCertificateFromRCVMIPSolution(
    CutCertificate& v,
    const std::vector<double>& solution,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const int cut_ind,
    FILE* const logfile);

/// @brief Use existing \p Atilde matrix (or recalculate it) to compute rank of submatrix given by the #CutCertificate \p v
RegularityStatus analyzeCertificateRegularity(
    int& certificate_rank,
    int& num_nonzero_multipliers,
    const CutCertificate& v,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const CoinPackedMatrix& Atilde,
    const int Atilderank,
    const StrengtheningParameters::Parameters& params);

/// @brief Given a set of \p cuts, identify regular ones and find their certificates to store in \p v
void analyzeCutRegularity(
    std::vector<CutCertificate>& v,
    std::vector<int>& certificate_submx_rank,
    std::vector<int>& num_nonzero_multipliers,
    std::vector<RegularityStatus>& regularity_status,
    std::vector<int>& num_iters,
    std::vector<double>& rcvmip_time,
    const OsiCuts& cuts,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const CoinPackedMatrix& Atilde,
    const int Atilderank,
    const StrengtheningParameters::Parameters& params,
    const bool USE_INPUT_CERTIFICATE = false);
