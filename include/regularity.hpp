/**
 * @file regularity.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#pragma once

#include <cstdio> // FILE
#include <vector>

namespace StrengtheningParameters {
  struct Parameters;
}

// COIN-OR files
class CoinPackedMatrix;
class OsiSolverInterface;
class OsiRowCut;
class OsiCuts;

class Disjunction;

/// [term][Farkas multiplier]
using CutCertificate = std::vector<std::vector<double> >; // this is also defined in analysis.hpp

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
    FILE* const logfile,
    const bool use_min_sum_delta = true);

/// @brief Given an existing RCVMILP, update the coefficients related to the cut
void updateRCVMILPFromCut(
    OsiSolverInterface* const liftingSolver,
    const OsiRowCut* const cut,
    const Disjunction* const disj, 
    const OsiSolverInterface* const solver);

/// @brief Given a set of \p cuts, identify regular ones and find their certificates to store in \p v
void analyzeRegularity(
  std::vector<CutCertificate>& v,
  std::vector<int>& certificate_submx_rank,
  const OsiCuts& cuts,
  const Disjunction* const disj,
  const OsiSolverInterface* const solver,
  const StrengtheningParameters::Parameters &params);
