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
#include "CutCertificate.hpp"
namespace StrengtheningParameters {
  struct Parameters;
}

class CglVPC; // CglVPC.hpp
class Disjunction; // Disjunction.hpp
struct SummaryBBInfo; // BBHelper.hpp

// Defined here
struct SummaryBoundInfo; // analysis.hpp
struct SummaryDisjunctionInfo; // analysis.hpp
struct SummaryCutInfo; // analysis.hpp
struct SummaryStrengtheningInfo; // analysis.hpp
struct SummaryCertificateInfo; // analysis.hpp

/// @brief Container for types of statistics we want to keep
enum class Stat { total = 0, avg, stddev, min, max, num_stats };

/// @brief Compute statistics in #Stat about given templated vector
template <typename T>
std::vector<double> computeStats(const std::vector<T>& v);

/// @brief Information about objective value at various points in the solution process
/// @details Gives objective for the LP and IP, and after adding GMICs, L&PCs, VPCs, and combinations of these cuts
/// and also keeps number of GMICs, L&PCs, and VPCs applied, and number of VPCs strengthened
struct SummaryBoundInfo {
  double lp_obj = std::numeric_limits<double>::max();
  double best_disj_obj = std::numeric_limits<double>::lowest();
  double worst_disj_obj = std::numeric_limits<double>::lowest();
  double root_obj = std::numeric_limits<double>::lowest();
  double ip_obj = std::numeric_limits<double>::max();
  double lpc_obj = std::numeric_limits<double>::max();

  double unstr_gmic_obj = std::numeric_limits<double>::max();
  double unstr_mycut_obj = std::numeric_limits<double>::max();
  double unstr_gmic_mycut_obj = std::numeric_limits<double>::max();
  double unstr_all_cuts_obj = std::numeric_limits<double>::max();
  
  double gmic_obj = std::numeric_limits<double>::max();
  double mycut_obj = std::numeric_limits<double>::max();
  double gmic_mycut_obj = std::numeric_limits<double>::max();
  double all_cuts_obj = std::numeric_limits<double>::max();
  
  double rcvmip_mycut_obj = std::numeric_limits<double>::max(); ///< objective after adding strengthened cuts using RCVMILP certificate
  double rcvmip_gmic_mycut_obj = std::numeric_limits<double>::max(); ///< objective after adding GMICs and strengthened cuts using RCVMILP certificate
  double rcvmip_all_cuts_obj = std::numeric_limits<double>::max(); ///< objective after adding all cuts, including strengthened cuts using RCVMILP certificate and strengthened with original certificate
  
  int num_root_bounds_changed = 0, num_gmic = 0, num_lpc = 0, num_mycut = 0;
  int num_str_affected_cuts = 0; ///< number of cuts for which strengthening using original certificate changes at least one coefficient
  int num_rcvmip_str_affected_cuts = 0; ///< number of cuts for which strengthening using RCVMIP certificate changes at least one coefficient
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

 /// @brief Track statistics for a particular family of cuts
struct SummaryCutInfo {
  int num_cuts = 0;
  int num_active_gmic = 0, num_active_lpc = 0, num_active_mycut, num_active_all = 0;
  int num_obj_tried = 0, num_failures = 0;
  int num_rounds = 0;
  int min_support = std::numeric_limits<int>::max();
  int max_support = 0;
  double avg_support = 0.;
  std::vector<CglAdvCut::CutType> cutType; ///< one entry per cut
  std::vector<CglAdvCut::ObjectiveType> objType; ///< one entry per cut

  std::vector<int> numCutsOfType;
  std::vector<int> numCutsFromHeur, numObjFromHeur, numFailsFromHeur, numActiveFromHeur;
  std::vector<int> numFails;
}; /* SummaryCutInfo */

/// @brief Summary statistics for strengthening attempts
struct SummaryStrengtheningInfo {
  /// number of cuts for which strengthening using the corresponding certificate changes at least one coefficient
  int num_str_affected_cuts = 0;
  /// number of coefficients changed (one entry for each index in #Stat -- total, avg, stddev, min, max)
  std::vector<double> num_coeffs_strengthened = std::vector<double>(static_cast<int>(Stat::num_stats), 0.);
}; /* SummaryStrengtheningInfo */

/// @brief Summary statistics for counting regularity / irregularity of cuts and certificates
struct SummaryCertificateInfo {
  std::vector<int> submx_rank; ///< rank of the submatrix of the certificate
  std::vector<int> num_nnz_mult; ///< number of original problem constraints with nonzero multipliers in certificate

  /// number of times certificate uses nonzero multipliers on both the upper and lower bounds on a variable
  int num_unmatched_bounds = 0;
  /// use multipliers to find number of distinct facets of the cut-generating set S using certificates
  double avg_num_cgs_facet = 0.;
  
  // /// number of certificates seemingly (but not proved) leading to submx rank less than number of nonzero multipliers used
  // int num_tentative_irreg_less = 0;
  /// number of certificates for which the first iteration has RCVMIP optimal value = 0, but there were rank constraints, and resolving was not possible (e.g., for limit reasons)
  int num_tentative_irreg_more = 0;
  /// number of certificates leading to submx rank less than number of nonzero multipliers used
  int num_irreg_less = 0;
  /// number of regular certificates (submx rank equals number of nonzero multipliers used)
  int num_reg = 0;
  /// number of certificates leading to submx rank more than number of nonzero multipliers used
  int num_irreg_more = 0;
  /// number of cuts for which certificate could not be ascertained
  int num_unconverged = 0;
  /// number of cuts for which certificate could not be computed due to numerical instability
  int num_numerically_unstable = 0;
  /// number of iterations needed to compute certificate
  std::vector<int> num_iterations;
  /// time spent per cut on RCVMIP
  std::vector<double> rcvmip_time;
}; /* SummaryCertificateInfo */

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
void printStrInfo(const SummaryStrengtheningInfo& orig_info, const SummaryStrengtheningInfo& rcvmip_info, FILE* const logfile,
    const char SEP = ',');
/// @brief Write to log statistics about the certificate investigation summarized in #SummaryRegularityInfo
void printCertificateInfo(const SummaryCertificateInfo& orig_info, const SummaryCertificateInfo& rcvmip_info, const int RCVMIP_ITER_LIMIT,
    FILE* const logfile, const char SEP = ',');
void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, const SummaryCutInfo& cutInfoUnstr,
    FILE* logfile, const char SEP = ',');

/// @brief Check cut density and update min/max support in \p cutInfo
int checkCutDensity(
    SummaryCutInfo& cutInfo,
    const OsiRowCut* const cut,
    const double EPS = 1e-14);

/// @brief Check cut activity in solver and return true if activity of cut equals rhs
bool checkCutActivity(
    const OsiSolverInterface* const solver,
    const OsiRowCut* const cut);

/// @brief Checks how many cuts remove a given feasible solution
int checkCutsAgainstFeasibleSolution(
    const OsiCuts& currCuts,
    const std::vector<double> ip_solution);

/// @brief Compute gap closed and active cuts
void analyzeStrength(const StrengtheningParameters::Parameters& params, 
    const OsiSolverInterface* const solver_gmic,
    const OsiSolverInterface* const solver_mycut,
    const OsiSolverInterface* const solver_all,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfo, 
    const OsiCuts* const gmics, const OsiCuts* const mycuts,
    const SummaryBoundInfo& boundInfo, std::string& output);
/// @brief Prepare summary of #SummaryBBInfo
void analyzeBB(const StrengtheningParameters::Parameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output);
/// @brief Get number rounds of SICs needed to meet bound from VPCs+SICs
double getNumGomoryRounds(const StrengtheningParameters::Parameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver);

/// @brief Update average number of terms, density, rows, cols, points, rays, and partial tree information if applicable
void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen);
/// @brief Add to cut information after a round, such as number of cuts, objectives, failures, etc.
void updateCutInfo(SummaryCutInfo& cutInfo, const CglAdvCut& gen, const OsiCuts* cuts = NULL, const double EPS = 1e-14);
/// @brief Use this to merge cut info from multiple rounds
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds, const SummaryCutInfo* const oldCutInfos);

/// @brief Find |K| (number of nonzero multipliers in \p v) and info about the implied convex cgs for this cut
void setCertificateInfo(
    SummaryCertificateInfo& info,
    const Disjunction* const disj,
    const std::vector<CutCertificate>& v_vec,
    const int num_rows, const int num_cols,
    const std::vector<int>& str_cut_ind,
    const int num_str_cuts,
    const double EPS);