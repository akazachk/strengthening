/**
 * @file main.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-14
 */

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <chrono> // for timing
#include <limits> // numeric_limits
#include <memory> // for smart pointers

// For option handling
#include <getopt.h> // getopt, getopt_long
#define no_argument 0
#define required_argument 1
#define optional_argument 2

// COIN-OR
#include <OsiCuts.hpp>
#include <CglGMI.hpp>
#include <OsiSolverInterface.hpp>

// Project files
#include "analysis.hpp" // SummaryBoundInfo, SummaryCutInfo, Stat, SummaryStrengtheningInfo, SummaryCertificateInfo, printing to logfile
#include "BBHelper.hpp"
#include "CglAdvCut.hpp"
#include "CutCertificate.hpp" // TermCutCertificate and CutCertificate
#include "CutHelper.hpp"
#include "Disjunction.hpp" // DisjunctiveTerm, Disjunction, getSolverForTerm
//#include "disjcuts.hpp"
#include "gmic.hpp"
#include "regularity.hpp"
#include "verify.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "strengthen.hpp"
#include "utility.hpp"

#ifdef USE_EIGEN
#include "eigen.hpp" // computeRank
#endif

// For disjInfo
#include "CglVPC.hpp"
#include "PartialBBDisjunction.hpp"
#include "PRLP.hpp"

#ifdef USE_GUROBI
#include "GurobiHelper.hpp" // for obtaining ip opt
#endif

#ifdef TRACE
#include "debug.hpp"
#endif

enum OverallTimeStats {
  INIT_SOLVE_TIME,
  GOMORY_GEN_TIME,
  GOMORY_APPLY_TIME,
  CUT_TOTAL_TIME,
  CUT_GEN_TIME,
  CUT_APPLY_TIME,
  VPC_INIT_SOLVE_TIME,
  VPC_DISJ_SETUP_TIME,
  VPC_DISJ_GEN_TIME,
  VPC_PRLP_SETUP_TIME,
  VPC_PRLP_SOLVE_TIME,
  VPC_GEN_CUTS_TIME,
  STR_TOTAL_TIME,
  STR_CALC_CERT_TIME,
  STR_APPLY_CERT_TIME,
  REG_TOTAL_TIME,
  REG_GEN_ATILDE_TIME,
  REG_RANK_ATILDE_TIME,
  REG_ANALYZE_ORIG_CERT_TIME,
  REG_CALC_CERT_TIME,
  REG_APPLY_CERT_TIME,
  BB_TIME,
  TOTAL_APPLY_TIME,
  TOTAL_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "INIT_SOLVE_TIME",
  "GOMORY_GEN_TIME",
  "GOMORY_APPLY_TIME",
  "CUT_TOTAL_TIME",
  "CUT_GEN_TIME",
  "CUT_APPLY_TIME",
  "VPC_INIT_SOLVE_TIME",
  "VPC_DISJ_SETUP_TIME",
  "VPC_DISJ_GEN_TIME",
  "VPC_PRLP_SETUP_TIME",
  "VPC_PRLP_SOLVE_TIME",
  "VPC_GEN_CUTS_TIME",
  "STR_TOTAL_TIME",
  "STR_CALC_CERT_TIME",
  "STR_APPLY_CERT_TIME",
  "REG_TOTAL_TIME",
  "REG_GEN_ATILDE_TIME",
  "REG_RANK_ATILDE_TIME",
  "REG_ANALYZE_ORIG_CERT_TIME",
  "REG_CALC_CERT_TIME",
  "REG_APPLY_CERT_TIME",
  "BB_TIME",
  "TOTAL_APPLY_TIME",
  "TOTAL_TIME"
}; /* OverallTimeStatsName */

// Main file variables
Parameters params;

OsiSolverInterface *solver;               ///< stores the cuts we actually want to "count" (i.e., the base LP from which we generate cuts in each round)
OsiSolverInterface *origSolver;           ///< original solver in case we wish to come back to it later
OsiSolverInterface* GMICSolver = NULL;    ///< only GMICs
OsiSolverInterface* strCutSolver = NULL;  ///< if GMICs count, this is only mycuts; otherwise, it is both GMICs and mycuts
OsiSolverInterface* allCutSolver = NULL;  ///< all generated cuts (same as solver if no GMICs, in which case not generated)

OsiCuts gmics, mycuts, rcvmip_cuts;

std::string dir = "", filename = "", instname = "", in_file_ext = "";

CglAdvCut::ExitReason exitReason;
TimeStats timer;
std::time_t start_time_t, end_time_t;
char start_time_string[25];

SummaryBoundInfo boundInfo;
std::vector<SummaryBoundInfo> boundInfoVec;
SummaryBBInfo info_nocuts, info_mycuts, info_allcuts;
SummaryDisjunctionInfo disjInfo;
std::vector<SummaryCutInfo> cutInfoVec;
SummaryCutInfo cutInfo, cutInfoGMICs, cutInfoUnstr;
SummaryStrengtheningInfo strInfo, rcvmipStrInfo;
std::vector<SummaryCertificateInfo> origCertInfoVec;
std::vector<SummaryCertificateInfo> rcvmipCertInfoVec;

// For output
std::string cut_output = "", bb_output = "";

#ifdef CODE_VERSION
const std::string CODE_VERSION_STRING = x_macro_to_string(CODE_VERSION);
#endif
#ifdef VPC_VERSION
const std::string VPC_VERSION_STRING = x_macro_to_string(VPC_VERSION);
#endif
#ifdef VPC_CBC_VERSION
const std::string CBC_VERSION_STRING = x_macro_to_string(VPC_CBC_VERSION);
#endif
#ifdef VPC_CLP_VERSION
const std::string CLP_VERSION_STRING = x_macro_to_string(VPC_CLP_VERSION);
#endif
#ifdef USE_GUROBI
#include <gurobi_c++.h>
const std::string GUROBI_VERSION_STRING = std::to_string(GRB_VERSION_MAJOR) + "." + std::to_string(GRB_VERSION_MINOR) + std::to_string(GRB_VERSION_TECHNICAL);
#endif // USE_GUROBI
#ifdef USE_CPLEX
#include "ilcplex/cpxconst.h"
const std::string CPLEX_VERSION_STRING = std::to_string(CPX_VERSION_VERSION) + "." + std::to_string(CPX_VERSION_RELEASE) + "." + std::to_string(CPX_VERSION_MODIFICATION);
#endif // USE_CPLEX

// Catch abort signal if it ever gets sent
/**
 * Catch a signal. You can output debugging info.
 * If you return from this function, and it was called
 * because abort() was called, your program will exit or crash anyway
 * (with a dialog box on Windows).
 */
#include <csignal>
void signal_handler_with_error_msg(int signal_number) {
  error_msg(errorstring, "Abort or seg fault message received. Signal number: %d.\n", signal_number);
  writeErrorToLog(errorstring, params.logfile);
  exit(1);
} /* signal_handler_with_error_msg */

int startUp(int argc, char** argv);
int processArgs(int argc, char** argv);
void initializeSolver(OsiSolverInterface* &solver);
int wrapUp(int retCode, int argc, char** argv);

/// @brief Update average number of terms, density, rows, cols, points, rays, and partial tree information if applicable
void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen);

/// @brief Strengthen \p currCuts using the Farkas multipliers that will be computed and stored in \p v
void strengtheningHelper(
    OsiCuts& currCuts,
    std::vector<CutCertificate>& v,
    std::vector<int>& str_cut_ind,
    SummaryStrengtheningInfo& strInfo,
    SummaryBoundInfo& boundInfo,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const std::vector<double>& ip_solution,
    const bool is_rcvmip = false
);

/// @brief Calculate the strengthening certificate (Farkas multipliers) \p v for \p currCuts
void calcStrengtheningCertificateHelper(
    const OsiCuts& currCuts,
    std::vector<CutCertificate>& v,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver
);

/// @brief Apply the strengthening certificate (Farkas multipliers) \p v to \p currCuts
void applyStrengtheningCertificateHelper(
    OsiCuts& currCuts,
    const std::vector<CutCertificate>& v,
    std::vector<int>& str_cut_ind,
    SummaryStrengtheningInfo& strInfo,
    SummaryBoundInfo& boundInfo,
    const Disjunction* const disj,
    const OsiSolverInterface* const solver,
    const std::vector<double>& ip_solution,
    const bool is_rcvmip = false
);

/****************** MAIN FUNCTION **********************/
int main(int argc, char** argv) {
  // Do this early in your program's initialization
  std::signal(SIGABRT, signal_handler_with_error_msg);
  std::signal(SIGSEGV, signal_handler_with_error_msg);
  assert( OverallTimeStatsName.size() == OverallTimeStats::NUM_TIME_STATS );

  //====================================================================================================//
  // Set up timing
  for (int t = 0; t < OverallTimeStats::NUM_TIME_STATS; t++) {
    timer.register_name(OverallTimeStatsName[t]);
  }

  //====================================================================================================//
  // Print welcome message, set up logfile
  timer.start_timer(OverallTimeStats::TOTAL_TIME);
  int status = startUp(argc, argv);
  if (status) { return status; }

  //====================================================================================================//
  // Set up solver and get initial solution
  initializeSolver(solver);
  timer.start_timer(OverallTimeStats::INIT_SOLVE_TIME);
  solver->initialSolve();
  if (!checkSolverOptimality(solver, false)) {
    error_msg(errorstring, "Unable to solve initial LP relaxation.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  timer.end_timer(OverallTimeStats::INIT_SOLVE_TIME);
  boundInfo.lp_obj = solver->getObjValue();

  /** DEBUG TESTING BARRIER METHOD {
    OsiClpSolverInterface* interiorSolver = dynamic_cast<OsiClpSolverInterface*>(solver->clone());
    interiorSolver->getModelPtr()->barrier(false);

    if (interiorSolver->isProvenOptimal()) {
      std::cout << "Original\tBarrier\n";
      for (int i = 0; i < solver->getNumCols(); i++) {
        std::cout << solver->getColSolution()[i] << "\t";
        std::cout << interiorSolver->getColSolution()[i] << "\n";
      }
    } else {
      std::cerr << "Barrier method does not result in optimal solution." << std::endl;
    }
    exit(1);
  } **/

  //====================================================================================================//
  // Get IP solution if requested
  timer.end_timer(OverallTimeStats::TOTAL_TIME);
  std::vector<double> ip_solution;
  if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT)) {
    printf("\n## Try to find optimal IP solution to compare against. ##\n");
    BBInfo tmp_bb_info;
#ifdef USE_GUROBI
    int strategy = params.get(BB_STRATEGY);
    // If no solfile provided, check if we can find it in the same directory as the instance
    if (params.get(SOLFILE).empty()) {
      std::string f_name = dir + "/" + instname + "_gurobi.mst.gz";
      if (fexists(f_name.c_str())) {
        printf("No solution provided, but solution file found at %s.\n", f_name.c_str());
        params.set(SOLFILE, f_name);
      }
    }
    // If solfile provided, enable use_bound
    if (!params.get(SOLFILE).empty() && !use_bb_option(strategy, StrengtheningParameters::BB_Strategy_Options::use_best_bound)) {
      strategy = enable_bb_option(strategy, StrengtheningParameters::BB_Strategy_Options::use_best_bound);
    }
    doBranchAndBoundWithGurobi(params, strategy,
        params.get(stringParam::FILENAME).c_str(),
        tmp_bb_info, boundInfo.ip_obj, &ip_solution);
#endif
#ifdef TRACE
    /*{/// DEBUG
    bool first = true;
    for (int i = 0; i < (int) ip_solution.size(); i++) {
      const double val = ip_solution[i];
      if (isZero(val)) continue;
      if (!first) printf(", ");
      else first = false;
      printf("(%d,%g)", i, val);
    }
    printf("\n");
    }*/ /// DEBUG
#endif
    if (ip_solution.size() > 0) {
      if (isInfinity(std::abs(boundInfo.ip_obj))) {
        const double ip_obj = dotProduct(ip_solution.data(), solver->getObjCoefficients(), solver->getNumCols());
        boundInfo.ip_obj = ip_obj;
        params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
        fprintf(stdout, "Best known IP objective value is %s.\n", stringValue(boundInfo.ip_obj, "%g").c_str());
      }
    }
  } // get ip opt
  timer.start_timer(OverallTimeStats::TOTAL_TIME);

  //====================================================================================================//
  // Save original solver in case we wish to come back to it later
  origSolver = solver->clone();
  if (!origSolver->isProvenOptimal()) {
    origSolver->initialSolve();
    checkSolverOptimality(origSolver, false);
  }

  // Also save copies for calculating other objective values
  // We only need these if cuts other than mycuts are generated
  // solver         ::  stores the cuts we actually want to "count" (i.e., the base LP from which we generate cuts in each round)
  // GMICSolver     ::  only GMICs
  // strCutSolver   ::  if GMICs count, this is only mycuts; otherwise, it is both GMICs and mycuts
  // allCutSolver   ::  all generated cuts (same as solver if no GMICs or RCVMIP-strengthened cuts, in which case not generated)
  const int GOMORY_OPTION = params.get(intParam::GOMORY);
  const int SHOULD_ANALYZE_REGULARITY = params.get(StrengtheningParameters::intParam::ANALYZE_REGULARITY);
  if (GOMORY_OPTION != 0) {
    GMICSolver = solver->clone();
    strCutSolver = solver->clone();
    allCutSolver = solver->clone();
  }

  //====================================================================================================//
  // Now do rounds of cuts, until a limit is reached (e.g., time, number failures, number cuts, or all rounds are exhausted)
  boundInfo.num_mycut = 0, boundInfo.num_gmic = 0, boundInfo.num_str_affected_cuts = 0, boundInfo.num_rcvmip_str_affected_cuts = 0;
  int num_rounds = params.get(ROUNDS);
  std::vector<OsiCuts> mycuts_by_round(num_rounds);
  cutInfoVec.resize(num_rounds);
  boundInfoVec.resize(num_rounds);
  origCertInfoVec.resize(num_rounds);
  rcvmipCertInfoVec.resize(num_rounds);
  int round_ind = 0;
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind+1, num_rounds);
    }

    // Save a copy of the solver with all the cuts from previous rounds, but none of the cuts from this round
    OsiSolverInterface* roundOrigSolver = (round_ind > 0) ? solver->clone() : origSolver;
    
    // Initialize all the boundInfo entries for this round
    boundInfoVec[round_ind].lp_obj = roundOrigSolver->getObjValue();
    boundInfoVec[round_ind].num_mycut = 0;
    boundInfoVec[round_ind].num_gmic = 0;
    boundInfoVec[round_ind].num_str_affected_cuts = 0;
    boundInfoVec[round_ind].num_rcvmip_str_affected_cuts = 0;

    timer.start_timer(OverallTimeStats::CUT_TOTAL_TIME);

    //====================================================================================================//
    // Generate Gomory cuts
    // > 0: apply cuts to solver (affects disjunctive cuts generated); < 0: apply cuts to strCutSolver only
    // Option 1: GglGMI
    // Option 2: custom generate intersection cuts, calculate Farkas certificate, do strengthening
    // Option 3: custom generate intersection cuts, calculate Farkas certificate, do closed-form strengthening
    OsiCuts currGMICs;
    if (GOMORY_OPTION != 0) {
      timer.start_timer(OverallTimeStats::GOMORY_GEN_TIME);
      generateGomoryCuts(currGMICs, solver, GOMORY_OPTION, params.get(intParam::STRENGTHEN), params.get(doubleConst::AWAY), params.get(doubleConst::DIFFEPS), params.logfile);
      timer.end_timer(OverallTimeStats::GOMORY_GEN_TIME);

      // Apply GMICs to solvers
      boundInfo.num_gmic += currGMICs.sizeCuts();
      boundInfoVec[round_ind].num_gmic += currGMICs.sizeCuts();
      gmics.insert(currGMICs);

      // For timing, only count first application, to solver with only GMICs,
      // to get a sense of how a single LP relaxation is affected
      timer.start_timer(OverallTimeStats::GOMORY_APPLY_TIME);
      applyCutsCustom(GMICSolver, currGMICs);
      timer.end_timer(OverallTimeStats::GOMORY_APPLY_TIME);

      // Update bound
      boundInfo.gmic_obj = GMICSolver->getObjValue();
      boundInfoVec[round_ind].gmic_obj = boundInfo.gmic_obj;

      if (GOMORY_OPTION > 0) {
        applyCutsCustom(solver, currGMICs);
      }
      applyCutsCustom(allCutSolver, currGMICs, params.logfile);
      // else if (GOMORY_OPTION < 0) {
      //   applyCutsCustom(CutSolver, currGMICs);
      // }
    } // Gomory cut generation

    // { /// DEBUG
    //   // testGomory(solver, params);
    // } /// DEBUG

    //====================================================================================================//
    // Now for more general cuts
    OsiCuts& currCuts = mycuts_by_round[round_ind];

    const bool SHOULD_GENERATE_CUTS = (params.get(intParam::DISJ_TERMS) != 0);
    double CURR_EPS = params.get(StrengtheningParameters::doubleParam::EPS);
    Disjunction* disj = NULL;
    if (SHOULD_GENERATE_CUTS) {
      timer.start_timer(OverallTimeStats::CUT_GEN_TIME);

      // User can replace generateCuts method in CglAdvCut with whatever method is desired
      CglAdvCut gen(params);

      // Store the initial solve time in order to set a baseline for the PRLP resolve time
      gen.timer.add_value(
          CglAdvCut::CutTimeStatsName[static_cast<int>(CglAdvCut::CutTimeStats::INIT_SOLVE_TIME)],
          timer.get_value(OverallTimeStats::INIT_SOLVE_TIME));

      // Generate disjunctive cuts
      gen.generateCuts(*solver, currCuts);
      CURR_EPS = gen.probData.EPS;

      // Update timing from underlying generator
      CglVPC* vpc = &(gen.gen);
      TimeStats vpctimer = vpc->timer;
      std::vector<CglVPC::VPCTimeStats> vpc_stats = {
        CglVPC::VPCTimeStats::INIT_SOLVE_TIME,
        CglVPC::VPCTimeStats::DISJ_SETUP_TIME,
        CglVPC::VPCTimeStats::DISJ_GEN_TIME,
        CglVPC::VPCTimeStats::PRLP_SETUP_TIME,
        CglVPC::VPCTimeStats::PRLP_SOLVE_TIME,
        CglVPC::VPCTimeStats::GEN_CUTS_TIME,
      };
      std::vector<OverallTimeStats> overall_stats = {
        OverallTimeStats::VPC_INIT_SOLVE_TIME,
        OverallTimeStats::VPC_DISJ_SETUP_TIME,
        OverallTimeStats::VPC_DISJ_GEN_TIME,
        OverallTimeStats::VPC_PRLP_SETUP_TIME,
        OverallTimeStats::VPC_PRLP_SOLVE_TIME,
        OverallTimeStats::VPC_GEN_CUTS_TIME
      };
      for (int i = 0; i < (int) vpc_stats.size(); i++) {
        const CglVPC::VPCTimeStats stat = vpc_stats[i];
        const OverallTimeStats overall_stat = overall_stats[i];
        const clock_t currtimevalue = vpctimer.get_value(CglVPC::VPCTimeStatsName[static_cast<int>(stat)]);
        timer.add_value(overall_stat, currtimevalue);
      }

      // timer.add_value(
      //     CglAdvCut::CutTimeStatsName[static_cast<int>(CglVPC::VPCTimeStats::DISJ_SETUP_TIME)],
      //     vpc->timer.get_value(CglVPC::VPCTimeStats::DISJ_SETUP_TIME));

      // Update statistics about the disjunction objective value and cuts
      exitReason = gen.exitReason;
      if (gen.gen.disj()) {
        disj = (gen.gen.disj())->clone();
        disjInfo.num_disj++;
        boundInfo.num_mycut += gen.num_cuts;
        boundInfoVec[round_ind].num_mycut += gen.num_cuts;
        if (boundInfo.best_disj_obj < disj->best_obj) {
          boundInfo.best_disj_obj = disj->best_obj;
        }
        if (boundInfoVec[round_ind].best_disj_obj < disj->best_obj) {
          boundInfoVec[round_ind].best_disj_obj = disj->best_obj;
        }
        if (boundInfo.worst_disj_obj < disj->worst_obj) {
          boundInfo.worst_disj_obj = disj->worst_obj;
        }
        if (boundInfoVec[round_ind].worst_disj_obj < disj->worst_obj) {
          boundInfoVec[round_ind].best_disj_obj = disj->worst_obj;
        }
        if (round_ind == 0) {
          boundInfo.num_root_bounds_changed = disj->common_changed_var.size();
          boundInfo.root_obj = disj->root_obj;
        }
        boundInfoVec[round_ind].num_root_bounds_changed = disj->common_changed_var.size();
        boundInfoVec[round_ind].root_obj = disj->root_obj;
      }
      updateDisjInfo(disjInfo, disjInfo.num_disj, gen.gen);
      updateCutInfo(cutInfoVec[round_ind], &gen);

      // Get density of unstr cuts
      int total_support = 0;
      for (int cut_ind = 0; cut_ind < currCuts.sizeCuts(); cut_ind++) {
          OsiRowCut* disjCut = currCuts.rowCutPtr(cut_ind);
          const CoinPackedVector lhs = disjCut->row();
          const int num_elem = lhs.getNumElements();
          /*for (const double coeff : strCutCoeff) {
            if (!isZero(coeff)) num_elem++;
          }*/
          if (num_elem < cutInfoUnstr.min_support)
            cutInfoUnstr.min_support = num_elem;
          if (num_elem > cutInfoUnstr.max_support)
            cutInfoUnstr.max_support = num_elem;
          total_support += num_elem;
      }
      cutInfoUnstr.avg_support = (double) total_support / currCuts.sizeCuts();
      cutInfoUnstr.num_cuts += currCuts.sizeCuts();
      
      timer.end_timer(OverallTimeStats::CUT_GEN_TIME);
    } // SHOULD_GENERATE_CUTS
    else { // else update cutInfo with blanks
      setCutInfo(cutInfoVec[round_ind], 0, NULL);
    }

#if 0
    fprintf(stdout, "\n## Printing custom cuts ##\n");
    for (int cut_ind = 0; cut_ind < currCuts.sizeCuts(); cut_ind++) {
      printf("## Cut %d ##\n", cut_ind);
      const OsiRowCut* const cut = currCuts.rowCutPtr(cut_ind);
      cut->print();
    }
    fprintf(stdout, "Finished printing custom cuts.\n\n");
#endif

    // Check cuts against ip solution
    if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && !ip_solution.empty()) {
      const int num_violated = checkCutsAgainstFeasibleSolution(currCuts, ip_solution);
      printf("\n## Number of cuts violating the IP solution: %d ##\n", num_violated);
    }

    //====================================================================================================//
    // Get Farkas certificate and do strengthening
    OsiCuts unstrCurrCuts(currCuts); // save unstrengthened cuts
    std::vector<int> str_cut_ind;   // indices of cuts that were strengthened
    std::vector<CutCertificate> v;  // [cut][term][Farkas multiplier] in the end, per term, this will be of dimension rows + disj term ineqs + cols

    timer.start_timer(OverallTimeStats::STR_TOTAL_TIME);
    strengtheningHelper(currCuts, v, str_cut_ind, strInfo, boundInfoVec[round_ind], disj, solver, ip_solution);
    boundInfo.num_str_affected_cuts += str_cut_ind.size();
    setCertificateInfo(origCertInfoVec[round_ind], disj, v, solver->getNumRows(), solver->getNumCols(), str_cut_ind, CURR_EPS);
    timer.end_timer(OverallTimeStats::STR_TOTAL_TIME);

    mycuts.insert(currCuts);

    timer.end_timer(OverallTimeStats::CUT_TOTAL_TIME);

    //====================================================================================================//
    // Analyze regularity and irregularity
    timer.start_timer(OverallTimeStats::REG_TOTAL_TIME);

    // To start, obtain coefficient matrix (after adding globally-valid inequalities)
    CoinPackedMatrix Atilde;
    std::vector<double> btilde;
    if (SHOULD_ANALYZE_REGULARITY != 0 && disj != NULL) {
      const int mtilde = calculateNumRowsAtilde(disj, solver);
      const int num_orig_rows = solver->getNumRows();
      const int num_common_rows = disj->common_changed_var.size();
      const int num_bound_rows = calculateNumFiniteBounds(solver);
      printf("\n## Preparing Atilde matrix to analyze regularity (matrix will have %d rows = %d original constraints, %d globally-valid inequalities, %d finite lower+upper bounds). ##\n",
          mtilde, num_orig_rows, num_common_rows, num_bound_rows);
      
      timer.start_timer(OverallTimeStats::REG_GEN_ATILDE_TIME);
      prepareAtilde(Atilde, btilde, disj, solver, params.logfile);
      timer.end_timer(OverallTimeStats::REG_GEN_ATILDE_TIME);

      printf("Finished preparing Atilde matrix.\n");
    }

    // Compute rank of Atilde
    const bool SHOULD_COMPUTE_RANK = (params.get(StrengtheningParameters::intParam::ATILDE_COMPUTE_RANK) >= 1)
      && (Atilde.getNumRows() > 0);
    if (SHOULD_COMPUTE_RANK) {
      printf("\n## Computing rank of Atilde matrix (%d rows, %d columns, %1.4f sparsity). ##\n",
          Atilde.getNumRows(), Atilde.getNumCols(),
          Atilde.getNumElements() / (double) (Atilde.getNumRows() * Atilde.getNumCols()));
    }
    timer.start_timer(OverallTimeStats::REG_RANK_ATILDE_TIME);
    const int Atilderank = SHOULD_COMPUTE_RANK ?
        computeRank(&Atilde, std::vector<int>(), std::vector<int>()) :
        solver->getNumCols();
    timer.end_timer(OverallTimeStats::REG_RANK_ATILDE_TIME);
    if (SHOULD_COMPUTE_RANK) {
      printf("Finished computing rank of Atilde matrix in %s seconds.\n",
          stringValue(timer.get_value(OverallTimeStats::REG_RANK_ATILDE_TIME), "%1.2f").c_str());
    }

    // Analyze regularity of the existing certificate for each cut
    std::vector<RegularityStatus> orig_regularity_status;
    if (SHOULD_ANALYZE_REGULARITY >= 1 && currCuts.sizeCuts() > 0 && disj) {
      timer.start_timer(OverallTimeStats::REG_ANALYZE_ORIG_CERT_TIME);

      // Calculate the rank of the existing certificate for each cut
      // (Do not compute regularity of the cut overall)
      printf("\n## Analyzing regularity of certificate computed for each cut. ##\n");
      orig_regularity_status.resize(currCuts.sizeCuts(), RegularityStatus::UNKNOWN);

      // For purposes of linear independence,
      // when a bound is used in the certificate on a certain variable
      // it does not matter if it is a lower or upper bound
      // since we can multiply rows by -1 without changing the rank
      // std::vector<int> certificate_submx_rank(currCuts.sizeCuts(), 0);
      // std::vector<int> num_nonzero_multipliers(solver->getNumRows() + disj->common_changed_var.size() + solver->getNumCols(), 0);
      origCertInfoVec[round_ind].submx_rank.resize(currCuts.sizeCuts(), 0);
      origCertInfoVec[round_ind].num_nnz_mult.resize(currCuts.sizeCuts(), 0);
      // assert(Atilde.getNumRows() == solver->getNumRows() + disj->common_changed_var.size());
      // assert(Atilde.getNumRows() + disj->terms[0].changed_var.size() + solver->getNumCols() == v[0][0].size()); // dimension matches for cut 0, term 0

      for (int cut_ind = 0; cut_ind < currCuts.sizeCuts(); cut_ind++) {
        int curr_submx_rank = -1;
        int curr_num_nnz_mult = -1;
        const RegularityStatus status = analyzeCertificateRegularity(
            curr_submx_rank, curr_num_nnz_mult, v[cut_ind],
            disj, solver, Atilde, Atilderank, params);
        origCertInfoVec[round_ind].submx_rank[cut_ind] = curr_submx_rank;
        origCertInfoVec[round_ind].num_nnz_mult[cut_ind] = curr_num_nnz_mult;

        switch (static_cast<int>(status)) {
          case static_cast<int>(RegularityStatus::IRREG_LESS):
            origCertInfoVec[round_ind].num_irreg_less++;
            break;
          case static_cast<int>(RegularityStatus::REG):
            origCertInfoVec[round_ind].num_reg++;
            orig_regularity_status[cut_ind] = RegularityStatus::REG; // only change for regular because others may be inaccurate
            break;
          case static_cast<int>(RegularityStatus::IRREG_MORE):
            origCertInfoVec[round_ind].num_irreg_more++;
            break;
          default:
            error_msg(errorstring, "Invalid status %d from origCertInfoVec for round %d cut %d.\n", static_cast<int>(status), round_ind, cut_ind);
            writeErrorToLog(errorstring, params.logfile);
            exit(1);
        } // switch on status
        
    #ifdef TRACE
        fprintf(stdout, "Cut %d: rank = %d/%d,\tnum_nonzero_multipliers = %d,\tregularity status = % d (%s)\n",
            cut_ind,
            origCertInfoVec[round_ind].submx_rank[cut_ind],
            Atilderank,
            origCertInfoVec[round_ind].num_nnz_mult[cut_ind],
            static_cast<int>(status),
            getRegularityStatusName(status).c_str());
    #endif
      } // loop over certificates, analyzing each for regularity

      timer.end_timer(OverallTimeStats::REG_ANALYZE_ORIG_CERT_TIME);
    } // analyze regularity of *certificate* (not cut overall)

    // Next, analyze regularity of cut overall, using all possible certificates, with the RCVMIP by Serra and Balas (2020)
    // If this is performed, we can obtain new strengthened cuts to compare against the existing ones
    std::vector<CutCertificate> rcvmip_v; // [cut][term][Farkas multiplier] in the end, per term, this will be of dimension rows + disj term ineqs + cols
    OsiCuts rcvmipCurrCuts;
    std::vector<int> rcvmip_str_cut_ind; // indices of cuts that were strengthened
    std::vector<RegularityStatus> rcvmip_regularity_status;

    if (SHOULD_ANALYZE_REGULARITY >= 2 && currCuts.sizeCuts() > 0 && disj) {
      printf("\n## Analyzing regularity of cuts via RCVMIP of Serra and Balas (2020). ##\n");

      // Copy to rcvmip_v the contents of v
      rcvmip_v = v;
      rcvmip_regularity_status = orig_regularity_status;
      rcvmipCertInfoVec[round_ind].submx_rank = origCertInfoVec[round_ind].submx_rank;
      rcvmipCertInfoVec[round_ind].num_nnz_mult = origCertInfoVec[round_ind].num_nnz_mult;

      timer.start_timer(OverallTimeStats::REG_CALC_CERT_TIME);
      analyzeCutRegularity(rcvmip_v,
          rcvmipCertInfoVec[round_ind].submx_rank,
          rcvmipCertInfoVec[round_ind].num_nnz_mult,
          rcvmip_regularity_status,
          rcvmipCertInfoVec[round_ind].num_iterations,
          rcvmipCertInfoVec[round_ind].rcvmip_time,
          unstrCurrCuts, disj, solver,
          Atilde, Atilderank, params);
      timer.end_timer(OverallTimeStats::REG_CALC_CERT_TIME);
    
      for (int cut_ind = 0; cut_ind < currCuts.sizeCuts(); cut_ind++) {
        const RegularityStatus status = rcvmip_regularity_status[cut_ind];

        switch (static_cast<int>(status)) {
          case static_cast<int>(RegularityStatus::IRREG_LESS):
            rcvmipCertInfoVec[round_ind].num_irreg_less++;
            break;
          case static_cast<int>(RegularityStatus::REG):
            rcvmipCertInfoVec[round_ind].num_reg++;
            break;
          case static_cast<int>(RegularityStatus::IRREG_MORE):
            rcvmipCertInfoVec[round_ind].num_irreg_more++;
            break;
          case static_cast<int>(RegularityStatus::UNCONVERGED):
            rcvmipCertInfoVec[round_ind].num_unconverged++;
            break;
          default:
            error_msg(errorstring, "Invalid status %d from rcvmipCertInfoVec for round %d cut %d.\n", static_cast<int>(status), round_ind, cut_ind);
            writeErrorToLog(errorstring, params.logfile);
            exit(1);
        } // switch on status

    #ifdef TRACE
        fprintf(stdout, "Cut %d: rank = %d/%d,\tnum_nonzero_multipliers = %d,\tregularity status = % d (%s)\n",
            cut_ind,
            rcvmipCertInfoVec[round_ind].submx_rank[cut_ind],
            Atilderank,
            rcvmipCertInfoVec[round_ind].num_nnz_mult[cut_ind],
            static_cast<int>(status),
            getRegularityStatusName(status).c_str());
    #endif
      } // loop over cuts, analyzing each for regularity

      // Apply the new strengthening certificates to get new cuts
      rcvmipCurrCuts = unstrCurrCuts; // assignment operator essentially inserts each of the unstrCurrCuts into rcvmipCurrCuts
      strengtheningHelper(rcvmipCurrCuts, rcvmip_v, rcvmip_str_cut_ind, rcvmipStrInfo, boundInfoVec[round_ind], disj, solver, ip_solution, false);
      setCertificateInfo(rcvmipCertInfoVec[round_ind], disj, rcvmip_v, solver->getNumRows(), solver->getNumCols(), rcvmip_str_cut_ind, CURR_EPS);
      boundInfo.num_rcvmip_str_affected_cuts += rcvmip_str_cut_ind.size();
      
      rcvmip_cuts.insert(rcvmipCurrCuts);

    } // analyze regularity of *cut* (not just certificate)
    timer.end_timer(OverallTimeStats::REG_TOTAL_TIME);
    
    //====================================================================================================//
    // Apply cuts
    printf("\n## Applying disjunctive cuts (# cuts = %d). ##\n", (int) currCuts.sizeCuts());
    timer.start_timer(OverallTimeStats::TOTAL_APPLY_TIME);

    // To track time to apply cuts, need to be careful of which solver they are being added to
    // If GOMORY = 0, then no GMICs generated; solver is only one that exists
    //   (strCutSolver and allCutSolver are *not* created; but if they were, it would hold that solver = strCutSolver = allCutSolver)
    // If GOMORY = 1, then solver has GMICs from this round applied already (will be identical to allCutSolver)
    // If GOMORY = -1, then solver has only disjcuts from prior rounds (will be identical to strCutSolver)
    // Recall that "solver" is the one used for generating cuts for each round
    if (currCuts.sizeCuts() > 0) {
      if (GOMORY_OPTION != 0) { // GMICs generated, so need to track strVPC-only objective with strCutSolver
        // Apply to "strCutSolver", which tracks effect of strengthened cuts, without the effect of other cuts such as GMICs
        timer.start_timer(OverallTimeStats::CUT_APPLY_TIME);
        applyCutsCustom(strCutSolver, currCuts, params.logfile);
        timer.end_timer(OverallTimeStats::CUT_APPLY_TIME);

        // Apply to "solver" which has either only strengthened cuts, or also GMICs depending on GOMORY parameter
        // When GOMORY = -1, strCutSolver and solver are identical
        applyCutsCustom(solver, currCuts, params.logfile);

        // Apply to "allCutSolver", which has all cuts
        // When GOMORY = 1, allCutSolver and solver are identical
        applyCutsCustom(allCutSolver, currCuts, params.logfile);
      } // check if GMICs generated (params.get(GOMORY) != 0)
      else {
        // Else, no GMICs generated; only solver needed
        timer.start_timer(OverallTimeStats::CUT_APPLY_TIME);
        applyCutsCustom(solver, currCuts, params.logfile);
        timer.end_timer(OverallTimeStats::CUT_APPLY_TIME);
      }
    } // apply cuts if any were generated

    // Update bound info for strengthened cuts
    boundInfo.mycut_obj      = (GOMORY_OPTION != 0) ? strCutSolver->getObjValue() : solver->getObjValue();
    boundInfo.gmic_mycut_obj = (GOMORY_OPTION != 0) ? allCutSolver->getObjValue() : solver->getObjValue();

    // Currently, all_cuts_obj will only use strengthened cuts from one certificate
    // To get the objective when various certificates are used, see rcvmip_all_cuts_obj
    boundInfo.all_cuts_obj = boundInfo.gmic_mycut_obj;

    // In addition, we measure effect of unstrengthened cuts on their own
    const bool unstr_str_cuts_different = (unstrCurrCuts.sizeCuts() > 0 && strInfo.num_str_affected_cuts > 0); // if false, unstr info is same as str info
    if (unstr_str_cuts_different) {
      OsiSolverInterface* unstrCutSolver = roundOrigSolver->clone();
      
      // Add unstrengthened cuts
      applyCutsCustom(unstrCutSolver, unstrCurrCuts, params.logfile);
      boundInfo.unstr_mycut_obj = unstrCutSolver->getObjValue();

      // Add GMICs
      if (currGMICs.sizeCuts() > 0) { applyCutsCustom(unstrCutSolver, currGMICs, params.logfile); }
      boundInfo.unstr_gmic_mycut_obj = unstrCutSolver->getObjValue();

      boundInfo.unstr_all_cuts_obj = boundInfo.unstr_gmic_mycut_obj;

      if (unstrCutSolver) { delete unstrCutSolver; }
    } else {
      boundInfo.unstr_mycut_obj = boundInfo.mycut_obj;
      boundInfo.unstr_gmic_mycut_obj = boundInfo.gmic_mycut_obj;
      boundInfo.unstr_all_cuts_obj = boundInfo.unstr_gmic_mycut_obj;
    }

    // Repeat for RCVMIP-strengthened cuts
    const bool rcvmip_cuts_different = (SHOULD_ANALYZE_REGULARITY > 1 && rcvmipCurrCuts.sizeCuts() > 0 && boundInfo.num_rcvmip_str_affected_cuts > 0); // if false, rcvmip info is same as unstr or str info
    if (rcvmip_cuts_different) {
      OsiSolverInterface* rcvmipCutSolver = roundOrigSolver->clone();
      
      // Add RCVMIP-strengthened cuts
      applyCutsCustom(rcvmipCutSolver, rcvmipCurrCuts, params.logfile);
      boundInfo.rcvmip_mycut_obj = rcvmipCutSolver->getObjValue();

      // Add GMICs
      if (currGMICs.sizeCuts() > 0) { applyCutsCustom(rcvmipCutSolver, currGMICs, params.logfile); }
      boundInfo.rcvmip_gmic_mycut_obj = rcvmipCutSolver->getObjValue();

      // Add original-strengthened cuts
      applyCutsCustom(rcvmipCutSolver, currCuts, params.logfile);
      boundInfo.rcvmip_all_cuts_obj = rcvmipCutSolver->getObjValue();

      if (rcvmipCutSolver) { delete rcvmipCutSolver; }
    } else {
      if (SHOULD_ANALYZE_REGULARITY <= 1) {
        boundInfo.rcvmip_mycut_obj = boundInfo.mycut_obj;
        boundInfo.rcvmip_gmic_mycut_obj = boundInfo.gmic_mycut_obj;
      } else {
        boundInfo.rcvmip_mycut_obj = boundInfo.unstr_mycut_obj;
        boundInfo.rcvmip_gmic_mycut_obj = boundInfo.unstr_gmic_mycut_obj;
      }
      boundInfo.rcvmip_all_cuts_obj = boundInfo.all_cuts_obj; // include strengthened cuts from other certificates
    }

    // Set boundInfoVec too for str, unstr, rcvmip cuts
    boundInfoVec[round_ind].mycut_obj      = boundInfo.mycut_obj;
    boundInfoVec[round_ind].gmic_mycut_obj = boundInfo.gmic_mycut_obj;
    boundInfoVec[round_ind].all_cuts_obj   = boundInfo.all_cuts_obj;
    boundInfoVec[round_ind].unstr_mycut_obj      = boundInfo.unstr_mycut_obj;
    boundInfoVec[round_ind].unstr_gmic_mycut_obj = boundInfo.unstr_gmic_mycut_obj;
    boundInfoVec[round_ind].unstr_all_cuts_obj   = boundInfo.unstr_all_cuts_obj;
    boundInfoVec[round_ind].rcvmip_mycut_obj      = boundInfo.rcvmip_mycut_obj;
    boundInfoVec[round_ind].rcvmip_gmic_mycut_obj = boundInfo.rcvmip_gmic_mycut_obj;
    boundInfoVec[round_ind].rcvmip_all_cuts_obj   = boundInfo.rcvmip_all_cuts_obj;
    
    timer.end_timer(OverallTimeStats::TOTAL_APPLY_TIME);

    // Free memory from solvers/disj specific for this round
    if (roundOrigSolver && roundOrigSolver != origSolver) { delete roundOrigSolver; }
    if (disj) { delete disj; }
    
    // Print summary of bounds after this round
    printf("\n*** INFO:\n");
    printf(
        "Round %d/%d: Completed round of cut generation (exit reason: %s).\n", 
        round_ind + 1, params.get(ROUNDS),
        CglAdvCut::ExitReasonName[static_cast<int>(exitReason)].c_str());
    printf(
        "Cuts generated = %d (%d total). GMICs generated = %d.\n",
        boundInfoVec[round_ind].num_mycut,
        boundInfo.num_mycut,
        boundInfoVec[round_ind].num_gmic);
    fflush(stdout);
    printf("Obj value (unstrengthened): %s.\n", 
        stringValue(boundInfo.unstr_all_cuts_obj, "%1.6f").c_str());
    fflush(stdout);
    printf("Obj value (%d/%d strengthened): %s.\n",
        strInfo.num_str_affected_cuts, boundInfoVec[round_ind].num_mycut,
        stringValue(boundInfo.all_cuts_obj, "%1.6f").c_str());
    fflush(stdout);
    if (SHOULD_ANALYZE_REGULARITY > 1 && !isInfinity(boundInfo.rcvmip_all_cuts_obj)) {
      printf("Obj value (%d/%d RCVMIP strengthened): %s.\n",
          boundInfo.num_rcvmip_str_affected_cuts, boundInfoVec[round_ind].num_mycut,
          stringValue(boundInfo.rcvmip_all_cuts_obj, "%1.6f").c_str());
    }
    fflush(stdout);
    printf("Initial obj value: %s.\n", 
        stringValue(boundInfo.lp_obj, "%1.6f").c_str());
    if (round_ind > 0) {
      printf("Round start obj value: %s.\n",
        stringValue(roundOrigSolver->getObjValue(), "%1.6f").c_str());
    }
    printf("Disj lb: %s.\n",
        stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());
    fflush(stdout);
    printf("Elapsed time: %1.2f seconds. Time limit = %1.2f seconds.\n",
        timer.get_total_time(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME]),
        params.get(TIMELIMIT));
    printf("***\n");
    
    // Exit early if reached time limit
    if (timer.reachedTimeLimit(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME], params.get(TIMELIMIT))) {
      printf("\n*** INFO: Reached time limit (current time = %f, time limit = %f).\n",
          timer.get_total_time(OverallTimeStatsName[OverallTimeStats::TOTAL_TIME]), params.get(TIMELIMIT));
      break;
    }

    // Exit early from rounds of cuts if no cuts generated or solver is not optimal
    if ((SHOULD_GENERATE_CUTS && boundInfoVec[round_ind].num_mycut == 0)
        || (!SHOULD_GENERATE_CUTS && boundInfoVec[round_ind].num_gmic == 0)
        || !solver->isProvenOptimal()
        || isInfinity(std::abs(solver->getObjValue()))) {
      break;
    }
  } // loop over rounds of cuts
  if (round_ind < num_rounds) {
    num_rounds = round_ind+1;
  }

  printf(
      "\n## Finished with %d cuts and %d GMICs. Initial obj value: %s.",
      boundInfo.num_mycut, boundInfo.num_gmic,
      stringValue(boundInfo.lp_obj, "%1.6f").c_str());
  if (strInfo.num_str_affected_cuts > 0) printf(
      " Unstrengthened obj value: %s.",
      stringValue(boundInfo.unstr_mycut_obj, "%1.6f").c_str());
  if (boundInfo.num_mycut > 0) printf(" Final strcut obj value: %s.",
      stringValue(boundInfo.mycut_obj, "%1.6f").c_str());
  printf(" Final all cuts obj value: %s.",
      stringValue(boundInfo.all_cuts_obj, "%1.6f").c_str());
  if (SHOULD_ANALYZE_REGULARITY > 1) {
    printf(" Final RCVMIP-cuts obj value: %s.",
        stringValue(boundInfo.rcvmip_all_cuts_obj, "%1.6f").c_str());
  }
  printf(" ##\n");
  fflush(stdout);

  //====================================================================================================//
  // Do branch-and-bound experiments (if requested)
  if (params.get(BB_RUNS) != 0) {
    // Collect cuts from all rounds
    timer.start_timer(BB_TIME);
    runBBTests(params, &info_nocuts, &info_mycuts, &info_allcuts,
        params.get(stringParam::FILENAME), origSolver, boundInfo.ip_obj,
        &mycuts, NULL, &ip_solution);
    timer.end_timer(BB_TIME);
  }

  timer.end_timer(OverallTimeStats::TOTAL_TIME);

#ifdef PRINT_LP_WITH_CUTS
  if (boundInfo.num_mycut > 0) {
    std::string fileWithCuts = filename + "_cuts";
    solver->writeMps(fileWithCuts.c_str());
  }
#endif

  //====================================================================================================//
  // Do analyses in preparation for printing
  setCutInfo(cutInfo, num_rounds, cutInfoVec.data());
  analyzeStrength(params,
      GOMORY_OPTION == 0 ? NULL : GMICSolver,
      GOMORY_OPTION <= 0 ? solver : strCutSolver, 
      GOMORY_OPTION >= 0 ? solver : allCutSolver, 
      cutInfoGMICs, cutInfo,
      GOMORY_OPTION == 0 ? NULL : &gmics, 
      &mycuts,
      boundInfo, cut_output);
  analyzeBB(params, info_nocuts, info_mycuts, info_allcuts, bb_output);

  //====================================================================================================//
  // Finish up
  return wrapUp(0, argc, argv);
} /* main */

/**
 * Call this early to print welcome message, etc.
 */
int startUp(int argc, char** argv) {
  int status = 0;

  // Input handling
  printf("## Custom Cut Generator ##\n");
  printf("# Aleksandr M. Kazachkov\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  time(&start_time_t);
  struct tm* start_timeinfo = localtime(&start_time_t);
  snprintf(start_time_string, sizeof(start_time_string) / sizeof(char), "%s", asctime(start_timeinfo));
  printf("Start time: %s\n", start_time_string);

  /*
  printf("\n## Version Information ##\n");
#ifdef VPC_VERSION
  printf("# VPC Version %s\n", VPC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("# Cbc Version %s\n", CBC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("# Clp Version %s\n", CLP_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef USE_GUROBI
  printf("# Gurobi Version %s\n", GUROBI_VERSION_STRING.c_str());
#endif
#ifdef USE_CPLEX
  printf("# CPLEX Version %s\n", CPLEX_VERSION_STRING.c_str());
#endif
  printf("\n");
  */

  status = processArgs(argc, argv);
  if (status) { return status; }

  // Get instance file
  if (params.get(stringParam::FILENAME).empty()) {
    error_msg(errorstring,
        "No instance file provided. Use -f /path/to/instance.mps to specify the instance.\n");
    exit(1);
  }
  printf("Instance file: %s\n", params.get(stringParam::FILENAME).c_str());
  
  if (parseFilename(dir, instname, in_file_ext, params.get(stringParam::FILENAME), params.logfile) != 0) {
    error_msg(errorstring,
        "Unable to parse filename: %s. Found: dir=\"%s\", instname=\"%s\",ext=\"%s\".\n",
        params.get(stringParam::FILENAME).c_str(), dir.c_str(),
        instname.c_str(), in_file_ext.c_str());
    exit(1);
  }
  filename = dir + "/" + instname;

  // Prepare logfile
  const std::string logname = params.get(stringParam::LOGFILE);
  bool logexists = false; 
  if (!logname.empty()) {
    logexists = fexists(logname.c_str());
    params.logfile = fopen(logname.c_str(), "a");
  }

  // Read opt value (if not yet inputted)
  if (!isInfinity(params.get(doubleParam::IP_OBJ))) {
    boundInfo.ip_obj = params.get(doubleParam::IP_OBJ);
  }
  if (isInfinity(boundInfo.ip_obj) && !params.get(stringParam::OPTFILE).empty()) {
    std::string optfile = params.get(stringParam::OPTFILE);
    std::string csvext = ".csv";
    if (optfile.size() > csvext.size() && optfile.compare(optfile.size()-csvext.size(), csvext.size(), csvext) == 0) {
  #ifdef TRACE
      fprintf(stdout, "Reading objective information from \"%s\".\n", optfile.c_str());
  #endif
      boundInfo.ip_obj = getObjValueFromFile(optfile, params.get(stringParam::FILENAME), params.logfile);
      params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
      fprintf(stdout, "Best known IP objective value is %s.\n", stringValue(boundInfo.ip_obj, "%f").c_str());
      if (isInfinity(boundInfo.ip_obj)) {
        warning_msg(warnstring, "Did not find objective value.\n");
      }
    }
  }
  if (params.logfile != NULL) {
    if (!logexists) {
      printHeader(params, OverallTimeStatsName);
    }
    fprintf(params.logfile, "%s,", instname.c_str());
    printParams(params, params.logfile, 2); // only values
    fflush(params.logfile);
  }

  return status;
} /* startUp */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode, int argc, char** argv) {
  const int exitReasonInt = static_cast<int>(exitReason);

  time(&end_time_t);
  struct tm* end_timeinfo = localtime(&end_time_t);
  char end_time_string[25];
  snprintf(end_time_string, sizeof(end_time_string) / sizeof(char), "%s", asctime(end_timeinfo));

  FILE* logfile = params.logfile;
  if (logfile != NULL) {
    {
      // Bound and gap info
      printBoundAndGapInfo(boundInfo, params.logfile);
      // B&B info
      printSummaryBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(BB_RUNS) == 0);
      // Orig prob
      printOrigProbInfo(origSolver, params.logfile);
      // Post-cut prob
      printPostCutProbInfo(solver, cutInfoGMICs, cutInfo, params.logfile);
      // Disj info
      printDisjInfo(disjInfo, params.logfile);
      // Str info
      printStrInfo(strInfo, rcvmipStrInfo, params.logfile);
      // Regularity info
      printCertificateInfo(origCertInfoVec[0], rcvmipCertInfoVec[0], params.get(StrengtheningParameters::intParam::RCVMIP_MAX_ITERS), params.logfile);
      // Cut, obj, fail info
      printCutInfo(cutInfoGMICs, cutInfo, cutInfoUnstr, params.logfile);
      // Full B&B info
      printFullBBInfo({info_nocuts, info_mycuts}, params.logfile, params.get(BB_RUNS) == 0);
      // Print time info
      timer.print(params.logfile, 2); // only values
    }

    // Print version information
#ifdef CODE_VERSION
    fprintf(logfile, "%s,", CODE_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_VERSION
    fprintf(logfile, "%s,", VPC_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CBC_VERSION
    fprintf(logfile, "%s,", CBC_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef VPC_CLP_VERSION
    fprintf(logfile, "%s,", CLP_VERSION_STRING.substr(0,8).c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_GUROBI
    fprintf(logfile, "%s,", GUROBI_VERSION_STRING.c_str());
#else
    fprintf(logfile, ",");
#endif
#ifdef USE_CPLEX
    fprintf(logfile, "%s,", CPLEX_VERSION_STRING.c_str());
#else
    fprintf(logfile, ",");
#endif

    // Print exit reason and finish
    fprintf(logfile, "%s,", CglAdvCut::ExitReasonName[exitReasonInt].c_str());
    fprintf(logfile, "%s,", end_time_string);
    fprintf(logfile, "%.2f,", difftime(end_time_t, start_time_t));
    fprintf(logfile, "%s,", instname.c_str());
    fprintf(logfile, "DONE\n");
    fclose(logfile); // closes params.logfile
  }

#ifdef TRACE
  // Print parameters
  printf("\n## Parameter values ##\n");
  printParams(params, stdout, 1);
  printf("\n");
  printParams(params, stdout, 2);
  printf("\n");

  int NAME_WIDTH = 25;
  // Print time info
  printf("\n## Time information ##\n");
  for (int i = 0; i < NUM_TIME_STATS; i++) {
    std::string name = OverallTimeStatsName[i];
    if (params.get(intParam::BB_RUNS) == 0 && name.compare(0, 2, "BB") == 0)
      continue;
    printf("%-*.*s%s\n", NAME_WIDTH, NAME_WIDTH, (name).c_str(),
        stringValue(timer.get_time(name), "%.3f").c_str());
  }
#endif

  // Print results from adding cuts
  printf("%s", cut_output.c_str());

  // Print branch-and-bound results
  printf("%s", bb_output.c_str());

  printf("\n## Exiting cut generation with reason %s. ##\n", CglAdvCut::ExitReasonName[exitReasonInt].c_str());
#ifdef CODE_VERSION
  printf("Code Version: %s\n", CODE_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_VERSION
  printf("VPC Version: %s\n", VPC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CBC_VERSION
  printf("Cbc Version: %s\n", CBC_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef VPC_CLP_VERSION
  printf("Clp Version: %s\n", CLP_VERSION_STRING.substr(0,8).c_str());
#endif
#ifdef USE_GUROBI
  printf("Gurobi Version: %s\n", GUROBI_VERSION_STRING.c_str());
#endif
#ifdef USE_CPLEX
  printf("CPLEX Version: %s\n", CPLEX_VERSION_STRING.c_str());
#endif
  printf("Instance: %s\n", instname.c_str());
  if (!params.get(stringParam::LOGFILE).empty()) {
    printf("Log: %s\n", params.get(stringParam::LOGFILE).c_str());
  } else {
    printf("Log: stdout\n");
  }
  printf("Start time: %s\n", start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));
  { // Print command used to repeat this run
    if (argc > 0) printf("Command: ");
    for (int i = 0; i < argc; i++) {
      printf("%s ", argv[i]);
    }
    if (argc > 0) printf("\n");
  }

  if (solver) {
    delete solver;
  }
  if (origSolver) {
    delete origSolver;
  }
  if (GMICSolver) {
    delete GMICSolver;
  }
  if (strCutSolver) {
    delete strCutSolver;
  }
  if (allCutSolver) {
    delete allCutSolver;
  }
  return retCode;
} /* wrapUp */

void initializeSolver(OsiSolverInterface* &solver) {
  // Generate cuts
  solver = new SolverInterface;
  setLPSolverParameters(solver, params.get(VERBOSITY));

  int status = 0;
  if (in_file_ext.compare("lp") == 0) {
#ifdef TRACE
    printf("\n## Reading LP file. ##\n");
#endif
    status = solver->readLp(params.get(stringParam::FILENAME).c_str());
  } else {
    if (in_file_ext.compare("mps") == 0) {
#ifdef TRACE
      printf("\n## Reading MPS file. ##\n");
#endif
      status = solver->readMps(params.get(stringParam::FILENAME).c_str());
    } else {
      try {
#ifdef TRACE
        printf("\n## Reading MPS file. ##\n");
#endif
        status = solver->readMps(params.get(stringParam::FILENAME).c_str());
      } catch (std::exception& e) {
        error_msg(errorstring, "Unrecognized extension: %s.\n",
            in_file_ext.c_str());
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    }
  } // read file
  if (status < 0) {
    error_msg(errorstring, "Unable to read in file %s.\n",
        params.get(stringParam::FILENAME).c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Make sure we are doing a minimization problem; this is just to make later
  // comparisons simpler (i.e., a higher LP obj after adding the cut is better).
  if (solver->getObjSense() < 1e-3) {
    printf(
        "\n## Detected maximization problem. Negating objective function to make it minimization. ##\n");
    solver->setObjSense(1.0);
    const double* obj = solver->getObjCoefficients();
    for (int col = 0; col < solver->getNumCols(); col++) {
      solver->setObjCoeff(col, -1. * obj[col]);
    }
    double objOffset = getObjOffset(solver);
    if (objOffset != 0.) {
      solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
    }
  }
} /* initializeSolver */

/**
 * See params.hpp for descriptions of the parameters
 */
int processArgs(int argc, char** argv) {
  // Handle inputs
  // struct option declared in getopt.h
  // name: name of the long option
  // has_arg: 0,1,2 for none, required, or optional
  // *flag: how results are returned; if NULL, getopt_long() returns val (e.g., can be the equivalent short option character), and o/w getopt_long() returns 0, and flag points to a var which is set to val if the option is found, but left unchanged if the option is not found
  // val: value to return, or to load into the variable pointed to by flag
  const char* const short_opts = "a:b:B:c:d:f:g:hi:l:m:o:r:s:t:v:";
  const struct option long_opts[] =
  {
      {"analyze_regularity",    required_argument, 0, 'a'},
      {"atilde_compute_rank",   required_argument, 0, 'a'*'1'},
      {"bb_runs",               required_argument, 0, 'b'},
      {"bb_mode",               required_argument, 0, 'b'*'2'},
      {"bb_strategy",           required_argument, 0, 'B'},
      {"cutlimit",              required_argument, 0, 'c'},
      {"disj_terms",            required_argument, 0, 'd'},
      {"file",                  required_argument, 0, 'f'},
      {"gomory",                required_argument, 0, 'g'},
      {"help",                  no_argument,       0, 'h'},
      {"ip_obj",                required_argument, 0, 'i'},
      {"logfile",               required_argument, 0, 'l'},
      {"mode",                  required_argument, 0, 'm'},
      {"optfile",               required_argument, 0, 'o'},
      {"rounds",                required_argument, 0, 'r'},
      {"rcvmip_max_iters",      required_argument, 0, 'r'*'1'},
      {"rcvmip_total_timelimit",required_argument, 0, 'r'*'2'},
      {"rcvmip_cut_timelimit",  required_argument, 0, 'r'*'3'},
      {"solfile",               required_argument, 0, 's'*'1'},
      {"strengthen",            required_argument, 0, 's'},
      {"temp",                  required_argument, 0, 't'*'1'},
      {"timelimit",             required_argument, 0, 't'},
      {"verbosity",             required_argument, 0, 'v'},
      {nullptr, no_argument, nullptr, 0}
  };

  int inp;
  int status = 0;
  while ((inp = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (inp) {
      case 'a': {
                  int val;
                  intParam param = intParam::ANALYZE_REGULARITY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'a'*'1': {
                  int val;
                  intParam param = intParam::ATILDE_COMPUTE_RANK;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'b': {
                  int val;
                  intParam param = intParam::BB_RUNS;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'b'*'2': {
                  int val;
                  intParam param = intParam::BB_MODE;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'B': {
                  int val;
                  intParam param = intParam::BB_STRATEGY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'c': {
                  int val;
                  intParam param = intParam::CUTLIMIT;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'd': {
                  int val;
                  intParam param = intParam::DISJ_TERMS;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'f': {
                  params.set(stringParam::FILENAME, optarg);
                  break;
                }
      case 'g': {
                  int val;
                  intParam param = intParam::GOMORY;
                  if (!parseInt(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'i': {
                 double val;
                 doubleParam param = doubleParam::IP_OBJ;
                 if (!parseDouble(optarg, val)) {
                   error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                   exit(1);
                 }
                 params.set(param, val);
                 break;
               }
      case 'l': {
                  params.set(stringParam::LOGFILE, optarg);
                  break;
                }
      case 'm': {
                 int val;
                 intParam param = intParam::MODE;
                 if (!parseInt(optarg, val)) {
                   error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                   exit(1);
                 }
                 params.set(param, val);
                 break;
               }
      case 'o': {
                  params.set(stringParam::OPTFILE, optarg);
                  break;
                }
      case 'r': {
                   int val;
                   intParam param = intParam::ROUNDS;
                   if (!parseInt(optarg, val)) {
                     error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                     exit(1);
                   }
                   params.set(param, val);
                   break;
                 }
      case 'r' * '1': {
                        int val;
                        intParam param = intParam::RCVMIP_MAX_ITERS;
                        if (!parseInt(optarg, val)) {
                          error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                          exit(1);
                        }
                        params.set(param, val);
                        break;
                      }

      case 'r' * '2': {
                        double val;
                        doubleParam param = doubleParam::RCVMIP_TOTAL_TIMELIMIT;
                        if (!parseDouble(optarg, val)) {
                          error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                          exit(1);
                        }
                        params.set(param, val);
                        break;
                      }

      case 'r' * '3': {
                        double val;
                        doubleParam param = doubleParam::RCVMIP_CUT_TIMELIMIT;
                        if (!parseDouble(optarg, val)) {
                          error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                          exit(1);
                        }
                        params.set(param, val);
                        break;
                      }
      case 's': {
                   int val;
                   intParam param = intParam::STRENGTHEN;
                   if (!parseInt(optarg, val)) {
                     error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                     exit(1);
                   }
                   params.set(param, val);
                   break;
                 }
      case 's'*'1': {
                  params.set(stringParam::SOLFILE, optarg);
                  break;
                }
      case 't'*'1': {
                      int val;
                      intParam param = intParam::TEMP;
                      if (!parseInt(optarg, val)) {
                        error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                        exit(1);
                      }
                      params.set(param, val);
                      break;
                    }
      case 't': {
                  double val;
                  doubleParam param = doubleParam::TIMELIMIT;
                  if (!parseDouble(optarg, val)) {
                    error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                    exit(1);
                  }
                  params.set(param, val);
                  break;
                }
      case 'v': {
                   int val;
                   intParam param = intParam::VERBOSITY;
                   if (!parseInt(optarg, val)) {
                     error_msg(errorstring, "Error reading %s. Given value: %s.\n", params.name(param).c_str(), optarg);
                     exit(1);
                   }
                   params.set(param, val);
                   break;
                 }
      case 'h':
      case '?':
      default: {
                 // print help
                std::string helpstring;
                helpstring += "\n## DESCRIPTION ##\n";
                helpstring += "Code for generating custom cuts.\n";
                helpstring += "\n## OPTIONS ##\n";
                helpstring += "-h, --help\n\tPrint this help message.\n";
                helpstring += "--temp\n\tSet temporary options (e.g., value of 1 = do preprocessing on instance).\n";
                helpstring += "\n# Input/output #\n";
                helpstring += "-f file, --file=file\n\tFilename.\n";
                helpstring += "-i val, --ip_obj=val\n\tValue of integer optimum for this instance (takes precedence over -o/--optfile).\n";
                helpstring += "-l logfile, --logfile=logfile\n\tWhere to print log messages.\n";
                helpstring += "-o optfile, --optfile=optfile\n\tWhere to find integer optimum value information (a csv file formatted as \"instance_name,value\" on each row).\n";
                helpstring += "--solfile=solfile\n\tWhere to find integer optimum solution information (e.g., a mst/sol file produced by Gurobi/CPLEX/etc).\n";
                helpstring += "-v level, --verbosity=level\n\tVerbosity level (0: print little, 1: let solver output be visible).\n";
                helpstring += "\n# General cut options #\n";
                helpstring += "-c num cuts, --cutlimit=num cuts\n\tMaximum number of cuts to generate (0+ = as given, -k = k * # fractional variables at root).\n";
                helpstring += "-d num terms, --disj_terms=num terms\n\tMaximum number of disjunctive terms or disjunctions to generate (depending on mode).\n";
                helpstring += "-g +/- 0-3, --gomory=+/- 0-3\n\t0: do not use Gomory cuts, 1: generate Gomory cuts via CglGMI, 2: generate Gomory cuts via gmic.cpp, 3: try closed-form strengthening (<0: only gen, >0: also apply to LP).\n";
                helpstring += "-m mode, --mode=mode\n\tDescription needs to be entered.\n";
                helpstring += "-r num rounds, --rounds=num rounds\n\tNumber of rounds of cuts to apply.\n";
                helpstring += "-s 0/1/2, --strengthen=0/1/2\n\tWhether to strengthen cuts.\n";
                helpstring += "-t num seconds, --timelimit=num seconds\n\tTotal number of seconds allotted for cut generation.\n";
                helpstring += "\n# Branch-and-bound options #\n";
                helpstring += "-b 0+ --bb_runs=0+\n\tNumber of branch-and-bound repeats.\n";
                helpstring += "-B strategy --bb_strategy=strategy\n\tBranch-and-bound strategy (see BBHelper.hpp).\n";
                helpstring += "--bb_mode={0,1,10,11,100,...,111}\n\tWhich branch-and-bound experiments to run (ones = no cuts, tens = mycuts, hundreds = gmics).\n";
                helpstring += "\n# Regularity options #\n";
                helpstring += "-a 0/1/2, --analyze_regularity=0/1/2\n\t0: no, 1: yes, only first certificate 2: yes, use MIP to check for alternate certificates.\n";
                helpstring += "--rcvmip_max_iters=num iters\n\tMaximum number of iterations for RCVMIP.\n";
                helpstring += "--rcvmip_total_timelimit=num seconds\n\tTotal number of seconds allotted for RCVMIP (0: infinity). When specified, supercedes RCVMIP_CUT_TIMELIMIT.\n";
                helpstring += "--rcvmip_cut_timelimit=num seconds\n\tNumber of seconds allotted for generating certificate per cut with RCVMIP (0: infinity).\n";
                helpstring += "## END OF HELP ##\n";
                std::cout << helpstring << std::endl;
                status = 1;
               }
    } // switch statement for input
  } // process args
  //for (int i = 0; i < argc; i++) {
  //  std::cout << argv[i] << " ";
  //}
  //std::cout << std::endl;
  return status;
} /* processArgs */

void updateDisjInfo(SummaryDisjunctionInfo& disjInfo, const int num_disj, const CglVPC& gen) {
  if (num_disj <= 0)
    return;
  const Disjunction* const disj = gen.getDisjunction();
  const PRLP* const prlp = gen.getPRLP();
  if (!prlp)
    return;
  disjInfo.num_integer_sol += !(disj->integer_sol.empty());
  disjInfo.avg_num_terms = (disjInfo.avg_num_terms * (num_disj - 1) + disj->num_terms) / num_disj;
  disjInfo.avg_density_prlp = (disjInfo.avg_density_prlp * (num_disj - 1) + prlp->density) / num_disj;
  disjInfo.avg_num_rows_prlp += (disjInfo.avg_num_rows_prlp * (num_disj - 1) + prlp->getNumRows()) / num_disj;
  disjInfo.avg_num_cols_prlp += (disjInfo.avg_num_cols_prlp * (num_disj - 1) + prlp->getNumCols()) / num_disj;
  disjInfo.avg_num_points_prlp += (disjInfo.avg_num_points_prlp * (num_disj - 1) + prlp->numPoints) / num_disj;
  disjInfo.avg_num_rays_prlp += (disjInfo.avg_num_rays_prlp * (num_disj - 1) + prlp->numRays) / num_disj;
  try {
    const PartialBBDisjunction* const partialDisj =
        dynamic_cast<const PartialBBDisjunction* const >(disj);
    disjInfo.avg_explored_nodes += (disjInfo.avg_explored_nodes * (num_disj - 1) + partialDisj->data.num_nodes_on_tree) / num_disj;
    disjInfo.avg_pruned_nodes += (disjInfo.avg_pruned_nodes * (num_disj - 1) + partialDisj->data.num_pruned_nodes) / num_disj;
    disjInfo.avg_min_depth += (disjInfo.avg_min_depth * (num_disj - 1) + partialDisj->data.min_node_depth) / num_disj;
    disjInfo.avg_max_depth += (disjInfo.avg_max_depth * (num_disj - 1) + partialDisj->data.max_node_depth) / num_disj;
  } catch (std::exception& e) {

  }
} /* updateDisjInfo */

void strengtheningHelper(
    OsiCuts& currCuts,                      ///< [in/out] The cuts to be strengthened (in place)
    std::vector<CutCertificate>& v,         ///< [in/out] Certificate of cuts, per term, of dimension rows + disj term ineqs + cols with indices [cut][term][Farkas multiplier]
    std::vector<int>& str_cut_ind,          ///< [in/out] Indices of cuts that were strengthened
    SummaryStrengtheningInfo& strInfo,      ///< [in/out] The summary info for the strengthening
    SummaryBoundInfo& boundInfo,            ///< [in/out] The summary info for the bound
    const Disjunction* const disj,          ///< [in] The disjunction that generated the cuts
    const OsiSolverInterface* const solver, ///< [in] The solver that generated the cuts
    const std::vector<double>& ip_solution, ///< [in] Feasible integer solution
    const bool is_rcvmip                    ///< [in] Whether we are using the original (false) or RCVMIP (true) certificate --- if false, certificate \p v is calculated using #calcStrengtheningCertificate
) {
  // Check if cuts should be strengthened
  bool do_strengthening = params.get(intParam::STRENGTHEN) >= 1; // the strengthening parameter is set
  do_strengthening = do_strengthening && disj && disj->terms.size() > 0; // a disjunction exists, and it has terms
  do_strengthening = do_strengthening && currCuts.sizeCuts() > 0; // cuts have been generated
  do_strengthening = do_strengthening && disj->integer_sol.size() == 0; // TODO right now (2021-05-22) we cannot handle integer-feasible solutions found during branching

  if (!do_strengthening) {
    return;
  }

  const int num_cuts = currCuts.sizeCuts();
  printf("\n## Strengthening disjunctive cuts: (# cuts = %d). ##\n", num_cuts);

  // Retrieve the certificate
  if (!is_rcvmip) {
    calcStrengtheningCertificateHelper(currCuts, v, disj, solver);
  }

  // Apply the certificate
  applyStrengtheningCertificateHelper(currCuts, v, str_cut_ind, strInfo, boundInfo, disj, solver, ip_solution);

  // Print the results
  fprintf(stdout, "\nFinished strengthening (%d / %d cuts affected).\n", strInfo.num_str_affected_cuts, boundInfo.num_mycut);
  fprintf(stdout, "Number coeffs changed:\n");
  fprintf(stdout, "\ttotal: %g\n", strInfo.num_coeffs_strengthened[(int) Stat::total]);
  fprintf(stdout, "\tavg: %g\n", strInfo.num_coeffs_strengthened[(int) Stat::avg]);
  fprintf(stdout, "\tstddev: %g\n", strInfo.num_coeffs_strengthened[(int) Stat::stddev]);
  fprintf(stdout, "\tmin: %g\n", strInfo.num_coeffs_strengthened[(int) Stat::min]);
  fprintf(stdout, "\tmax: %g\n", strInfo.num_coeffs_strengthened[(int) Stat::max]);
  fprintf(stdout, "--------------------------------------------------\n");
#if 0
  fprintf(stdout, "\n## Printing strengthened custom cuts ##\n");
  for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
    printf("## Cut %d ##\n", cut_ind);
    const OsiRowCut* const cut = currCuts.rowCutPtr(cut_ind);
    cut->print();
  }
  fprintf(stdout, "Finished printing strengthened custom cuts.\n\n");
#endif
} /* strengtheningHelper */

void calcStrengtheningCertificateHelper(
    const OsiCuts& currCuts,                ///< [in] The cuts to be strengthened (in place)
    std::vector<CutCertificate>& v,         ///< [out] Certifcate of cuts that, in the end, per term, this will be of dimension rows + disj term ineqs + cols with indices [cut][term][Farkas multiplier]
    const Disjunction* const disj,          ///< [in] The disjunction that generated the cuts
    const OsiSolverInterface* const solver  ///< [in] The solver that generated the cuts
) {
  timer.start_timer(OverallTimeStats::STR_CALC_CERT_TIME);

  const int num_cuts = currCuts.sizeCuts();
  v.resize(num_cuts);
  for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
    v[cut_ind].resize(disj->num_terms);
  }
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    OsiSolverInterface* termSolver;
    disj->getSolverForTerm(termSolver, term_ind, solver, false, params.get(StrengtheningParameters::doubleConst::DIFFEPS), params.logfile);
    if (!termSolver) {
      printf("Disjunctive term %d/%d not created successfully.\n", term_ind+1, disj->num_terms);
      continue;
    }
    for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
      const OsiRowCut* disjCut = currCuts.rowCutPtr(cut_ind);

      // For each term of the disjunction,
      // we need to explicitly add the constraint(s) defining the disjunctive term
      const CoinPackedVector lhs = disjCut->row();

      getCertificate(v[cut_ind][term_ind], lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), termSolver, params.logfile);
    } // loop over cuts

    if (termSolver) { delete termSolver; }

  } // loop over disjunctive terms to retrieve certificate

  timer.end_timer(OverallTimeStats::STR_CALC_CERT_TIME);
} /* calcStrengtheningCertificateHelper */

void applyStrengtheningCertificateHelper(
    OsiCuts& currCuts,                      ///< [in/out] The cuts to be strengthened (in place)
    const std::vector<CutCertificate>& v,   ///< [in] Certifcate of cuts that, in the end, per term, this will be of dimension rows + disj term ineqs + cols with indices [cut][term][Farkas multiplier]
    std::vector<int>& str_cut_ind,          ///< [in/out] Indices of cuts that were strengthened
    SummaryStrengtheningInfo& strInfo,      ///< [in/out] The summary info for the strengthening
    SummaryBoundInfo& boundInfo,            ///< [in/out] The summary info for the bound
    const Disjunction* const disj,          ///< [in] The disjunction that generated the cuts
    const OsiSolverInterface* const solver, ///< [in] The solver that generated the cuts
    const std::vector<double>& ip_solution, ///< [in] Feasible integer solution
    const bool is_rcvmip                    ///< [in] Whether we are using the original (false) or RCVMIP (true) certificate
) {
  OverallTimeStats which_time = is_rcvmip ? OverallTimeStats::REG_APPLY_CERT_TIME : OverallTimeStats::STR_APPLY_CERT_TIME;
  timer.start_timer(which_time);

  const int num_cuts = currCuts.sizeCuts();
  strInfo.num_coeffs_strengthened.resize(static_cast<int>(Stat::num_stats), 0.); // total,avg,stddev,min,max
  strInfo.num_coeffs_strengthened[(int) Stat::min] = std::numeric_limits<int>::max();
  str_cut_ind.reserve(num_cuts);

  for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
    OsiRowCut* disjCut = currCuts.rowCutPtr(cut_ind);
    const CoinPackedVector lhs = disjCut->row();
    const double rhs = disjCut->rhs();
    std::vector<double> str_coeff;
    double str_rhs;
    const int curr_num_coeffs_str = strengthenCut(str_coeff, str_rhs,
        lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), rhs,
        disj, v[cut_ind], solver, params.logfile, ip_solution);

    // Update stats
    if (is_rcvmip) {
      boundInfo.num_rcvmip_str_affected_cuts += curr_num_coeffs_str > 0;
    } else {
      boundInfo.num_str_affected_cuts += curr_num_coeffs_str > 0;
    }
    strInfo.num_str_affected_cuts += curr_num_coeffs_str > 0;
    strInfo.num_coeffs_strengthened[(int) Stat::total] += curr_num_coeffs_str;
    strInfo.num_coeffs_strengthened[(int) Stat::avg] += (double) curr_num_coeffs_str / num_cuts;
    strInfo.num_coeffs_strengthened[(int) Stat::stddev] += (double) curr_num_coeffs_str * curr_num_coeffs_str / num_cuts;
    if (curr_num_coeffs_str < strInfo.num_coeffs_strengthened[(int) Stat::min]) {
      strInfo.num_coeffs_strengthened[(int) Stat::min] = curr_num_coeffs_str;
    }
    if (curr_num_coeffs_str > strInfo.num_coeffs_strengthened[(int) Stat::max]) {
      strInfo.num_coeffs_strengthened[(int) Stat::max] = curr_num_coeffs_str;
    }
    
    // Replace row if any coefficients were strengthened
    if (curr_num_coeffs_str > 0) {
      CoinPackedVector strCutCoeff(str_coeff.size(), str_coeff.data());
      disjCut->setRow(strCutCoeff);
      disjCut->setLb(str_rhs);
      str_cut_ind.push_back(cut_ind);
    }
  } // loop over cuts

  // Finish stddev calculation = sqrt(E[X^2] - E[X]^2)
  strInfo.num_coeffs_strengthened[(int) Stat::stddev] -= strInfo.num_coeffs_strengthened[(int) Stat::avg] * strInfo.num_coeffs_strengthened[(int) Stat::avg];
  strInfo.num_coeffs_strengthened[(int) Stat::stddev] = (strInfo.num_coeffs_strengthened[(int) Stat::stddev] > 0) ? std::sqrt(strInfo.num_coeffs_strengthened[(int) Stat::stddev]) : 0.;
  
  timer.end_timer(which_time);
} /* applyStrengtheningCertificateHelper */
