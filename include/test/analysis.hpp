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
struct SummaryBoundInfo {
  double lp_obj = std::numeric_limits<double>::max();
  double best_disj_obj = std::numeric_limits<double>::lowest();
  double worst_disj_obj = std::numeric_limits<double>::lowest();
  double ip_obj = std::numeric_limits<double>::max();
  double gmic_obj = std::numeric_limits<double>::max();
  double lpc_obj = std::numeric_limits<double>::max();
  double mycut_obj = std::numeric_limits<double>::max();
  double gmic_mycut_obj = std::numeric_limits<double>::max();
  double all_cuts_obj = std::numeric_limits<double>::max();
  int num_gmic = 0, num_lpc = 0, num_mycut = 0;
}; /* SummaryBoundInfo */

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
void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, FILE* logfile, const char SEP = ',');

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

