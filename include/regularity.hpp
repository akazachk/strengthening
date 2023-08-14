/**
 * @file regularity.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#pragma once

#include <cstdio> // FILE
#include <vector>

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
void genRCVMILPFromCut(
    OsiSolverInterface* const liftingSolver,
    const OsiRowCut* const cut,
    const Disjunction* const disj, 
    const OsiSolverInterface* const solver,
    const StrengtheningParameters::Parameters& params,
    const bool use_min_sum_delta = false);

/// @brief Given an existing RCVMILP, update the coefficients related to the cut
void updateRCVMILPFromCut(
    OsiSolverInterface* const liftingSolver,
    const OsiRowCut* const cut,
    const Disjunction* const disj, 
    const OsiSolverInterface* const solver);

/// @brief Obtain solution to RCVMILP
int solveRCVMILP(
  OsiSolverInterface* const liftingSolver,
  std::vector<double>& solution, 
  const StrengtheningParameters::Parameters& params,
  const int cut_ind);

/// @brief Given a solution to the RCVMILP, extract the Farkas multipliers
void getCertificateFromRCVMILPSolution(
  CutCertificate& v,
  const std::vector<double>& solution,
  const Disjunction* const disj,
  const OsiSolverInterface* const solver,
  const int cut_ind,
  FILE* const logfile);

/// @brief Use existing \p Atilde matrix (or recalculate it) to compute rank of submatrix given by the #CutCertificate \p v
void analyzeCertificateRegularity(
  int& certificate_rank,
  int& num_nonzero_multipliers,
  const CutCertificate& v,
  const Disjunction* const disj,
  const OsiSolverInterface* const solver,
  const CoinPackedMatrix& Atilde,
  const StrengtheningParameters::Parameters& params);

/// @brief Given a set of \p cuts, identify regular ones and find their certificates to store in \p v
void analyzeCutRegularity(
  std::vector<CutCertificate>& v,
  std::vector<int>& certificate_submx_rank,
  std::vector<int>& num_nonzero_multipliers,
  const OsiCuts& cuts,
  const Disjunction* const disj,
  const OsiSolverInterface* const solver,
  const StrengtheningParameters::Parameters &params);
