/**
 * @file analysis.hpp
 * @author A. M. Kazachkov
 * @date 2019-11-24
 */
#pragma once

#include <string>
#include <vector>
#include <limits>

class OsiSolverInterface;
class OsiCuts;

#include "CglAdvCut.hpp" // CutType, ObjectiveType
namespace StrengtheningParameters {
  struct Parameters;
}

struct SummaryBBInfo; // BBHelper.hpp

/// @brief Information about objective value at various points in the solution process
/// @details Gives objective for the LP and IP, and after adding GMICs, L&PCs, VPCs, and combinations of these cuts
/// and also keeps number of GMICs, L&PCs, and VPCs applied, and number of VPCs strengthened
struct SummaryBoundInfo {
  double lp_obj = std::numeric_limits<double>::max();
  double best_disj_obj = std::numeric_limits<double>::lowest();
  double worst_disj_obj = std::numeric_limits<double>::lowest();
  double ip_obj = std::numeric_limits<double>::max();
  double gmic_obj = std::numeric_limits<double>::max();
  double lpc_obj = std::numeric_limits<double>::max();
  double unstr_mycut_obj = std::numeric_limits<double>::max();
  double unstr_gmic_mycut_obj = std::numeric_limits<double>::max();
  double unstr_all_cuts_obj = std::numeric_limits<double>::max();
  double mycut_obj = std::numeric_limits<double>::max();
  double gmic_mycut_obj = std::numeric_limits<double>::max();
  double all_cuts_obj = std::numeric_limits<double>::max();
  int num_gmic = 0, num_lpc = 0, num_mycut = 0, num_str_cuts = 0;
}; /* SummaryBoundInfo */

/// @brief Summary statistics for the disjunction generated
struct SummaryDisjunctionInfo {
  int num_disj = 0;
  int num_integer_sol = 0;
  double avg_num_terms = 0;
  double avg_density_prlp = 0.;
  double avg_num_rows_prlp = 0.;
  double avg_num_cols_prlp = 0.;
  double avg_num_points_prlp = 0.;
  double avg_num_rays_prlp = 0.;
  double avg_explored_nodes = 0.;
  double avg_pruned_nodes = 0.;
  double avg_min_depth = 0.;
  double avg_max_depth = 0.;
}; /* SummaryDisjunctionInfo */

struct SummaryCutInfo {
  int num_cuts = 0;
  int num_active_gmic = 0, num_active_lpc = 0, num_active_mycut, num_active_all = 0;
  int num_obj_tried = 0, num_failures = 0;
  int num_rounds = 0;
  int min_support = std::numeric_limits<int>::max();
  int max_support = 0;
  double avg_support = 0.;
  std::vector<CglAdvCut::CutType> cutType; // one entry per cut
  std::vector<CglAdvCut::ObjectiveType> objType; // one entry per cut

  std::vector<int> numCutsOfType;
  std::vector<int> numCutsFromHeur, numObjFromHeur, numFailsFromHeur, numActiveFromHeur;
  std::vector<int> numFails;
}; /* SummaryCutInfo */

enum class Stat { total = 0, avg, stddev, min, max, num_stats };
struct SummaryStrengtheningInfo {
  /// number of cuts affected by strengthening
  int num_str_cuts = 0;
  /// use multipliers to find number of distinct facets of the cut-generating set S
  double avg_num_cgs_facets = 0.;
  /// number of cuts with less than n nonzero mutipliers
  int num_irreg_less = 0;
  /// number of cuts with more than n nonzero mutipliers
  int num_irreg_more = 0;
  /// number of coefficients changed (one entry for each index in #Stat -- total, avg, stddev, min, max)
  std::vector<double> num_coeffs_strengthened = std::vector<double>(static_cast<int>(Stat::num_stats), 0.);
}; /* SummaryStrengtheningInfo */

void printHeader(const StrengtheningParameters::Parameters& params,
    const std::vector<std::string>& time_name,
    const char SEP = ',');
void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile,
    const char SEP = ',');
void printSummaryBBInfo(const std::vector<SummaryBBInfo>& info, FILE* myfile,
    const bool print_blanks = false, const char SEP = ',');
void printFullBBInfo(const std::vector<SummaryBBInfo>& info, FILE* myfile,
    const bool print_blanks = false, const char SEP = ',');
void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP = ',');
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const SummaryCutInfo& cutInfoGMICs, const SummaryCutInfo& cutInfo,
    FILE* logfile, const char SEP = ',');
/// @brief Write to log statistics about the disjunction stored in #SummaryDisjunctionInfo
void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile,
    const char SEP = ',');
/// @brief Write to log statistics about the strengthening that was done as stored in #SummaryStrengtheningInfo
void printStrInfo(const SummaryStrengtheningInfo& info, FILE* myfile,
    const char SEP = ',');
void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, const SummaryCutInfo& cutInfoUnstr,
    FILE* logfile, const char SEP = ',');

void analyzeStrength(const StrengtheningParameters::Parameters& params, 
    const OsiSolverInterface* const solver_gmic,
    const OsiSolverInterface* const solver_mycut,
    const OsiSolverInterface* const solver_all,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfo, 
    const OsiCuts* const gmics, const OsiCuts* const mycuts,
    const SummaryBoundInfo& boundInfo, std::string& output);
void analyzeBB(const StrengtheningParameters::Parameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
double getNumGomoryRounds(const StrengtheningParameters::Parameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver);

void updateCutInfo(SummaryCutInfo& cutInfo, const CglAdvCut* const gen);
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds, const SummaryCutInfo* const oldCutInfos);

