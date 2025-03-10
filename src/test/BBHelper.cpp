/**
 * @file BBHelper.cpp
 * @brief Helper functions for branch-and-bound
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */
#include <cstdio> // for tmpnam
#include <chrono>
#include <random> // random_device, mt19937 (mersenne twister engine), uniform_real_distribution

// Project files 
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCutsCustom
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;
#include "TimeStats.hpp"
#include "utility.hpp"

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>

#ifdef USE_CBC
// Cbc
#include <CbcModel.hpp>

// Variable selection
#include <CbcBranchDefaultDecision.hpp>
#include <OsiChooseVariable.hpp>

// General strategy
#include <CbcStrategy.hpp>
#endif /* USE_CBC */

#ifdef USE_GUROBI
#include "GurobiHelper.hpp"
#endif
#ifdef USE_CPLEX
#include "CplexHelper.hpp"
#endif

/**
 * Solver is changed unless we are doing row/col permutations
 */
void runBBTests(
    /// [in] Initial set of parameters (will be cloned to change random seeds)
    const Parameters& base_params,
    /// [out] Store information about branch-and-bound results with no (user-generated) cuts added
    SummaryBBInfo* const info_nocuts,
    /// [out] Store information about branch-and-bound results with only my user-generated cuts added
    SummaryBBInfo* const info_mycuts,
    /// [out] Store information about branch-and-bound results with all user-generated cuts added (including GMICs, for example)
    SummaryBBInfo* const info_allcuts,
    /// [in] Full path to instance
    const std::string fullfilename,
    /// [in] Original model
    const OsiSolverInterface* const solver,
    /// [in] Best bound (to provide to IP solver, if relevant option is selected)
    const double best_bound,
    /// [in] My cuts to add to solver
    const OsiCuts* mycuts,
    /// [in] (Optional) GMICs to test against
    const OsiCuts* const gmics,
    /// [in/out] IP solution to original problem (will usually be empty unless you are debugging; will be computed here if requested)
    std::vector<double>* const ip_solution) {
  Parameters params = base_params;
  const int num_mycuts = (mycuts != NULL) ? mycuts->sizeCuts() : 0;
  // TODO Set number of b&b runs to be zero if no cuts generated (need to enable though if we just want to do b&b)
  const int num_bb_runs = std::abs(params.get(intParam::BB_RUNS)); //* (num_mycuts > 0);
  if (num_bb_runs == 0)
    return;

  /***********************************************************************************
   * Perform branch-and-bound
   ***********************************************************************************/
  const int num_gmics = (gmics != NULL) ? gmics->sizeCuts() : 0;
  const int numCutsToAddPerRound = (num_gmics > 0) ? num_gmics : num_mycuts;

  // B&B mode: ones bit = no_cuts, tens bit = w/mycuts, hundreds bit = w/gmics
  const int mode_param = params.get(intParam::BB_MODE);
  const int mode_ones = mode_param % 10;
  const int mode_tens = (mode_param % 100 - (mode_param % 10)) / 10;
  const int mode_hundreds = (mode_param % 1000 - (mode_param % 100)) / 100;
  const bool branch_with_no_cuts = (mode_ones > 0);
  const bool branch_with_mycuts = (mode_tens > 0) && (num_mycuts > 0);
  const bool branch_with_gmics = (mode_hundreds > 0) && (num_gmics > 0);
  if (branch_with_no_cuts + branch_with_mycuts + branch_with_gmics == 0)
    return;

  const bool should_permute_rows_and_cols = (num_bb_runs >= 2)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)
      && !use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex);

#ifdef TRACE
  printf("\n## Performing branch-and-bound tests. ##\n");
#endif

  if (branch_with_no_cuts) {
    assert(info_nocuts != NULL);
    info_nocuts->vec_bb_info.resize(0);
    info_nocuts->vec_bb_info.resize(num_bb_runs);
  }
  if (branch_with_mycuts) {
    assert(info_mycuts != NULL);
    info_mycuts->vec_bb_info.resize(0);
    info_mycuts->num_cuts = num_mycuts;
    info_mycuts->vec_bb_info.resize(num_bb_runs);
  }
  if (branch_with_gmics) {
    assert(info_allcuts != NULL);
    info_allcuts->num_cuts = num_mycuts + num_gmics;
    info_allcuts->vec_bb_info.resize(num_bb_runs);
  }
  std::vector<int> row_permutation, col_permutation;
  if (should_permute_rows_and_cols) {
    row_permutation.resize(solver->getNumRows());
    col_permutation.resize(solver->getNumRows());
    for (int row = 0; row < solver->getNumRows(); row++) {
      row_permutation[row] = row;
    }
    for (int col = 0; col < solver->getNumCols(); col++) {
      col_permutation[col] = col;
    }
  }

  OsiCuts GMICsAndmycuts;
  if (branch_with_mycuts)
    GMICsAndmycuts.insert(*mycuts);
  if (branch_with_gmics) {
    GMICsAndmycuts.insert(*gmics);
  }

  OsiSolverInterface* runSolver = NULL;
  const int given_seed = params.get(intParam::RANDOM_SEED);
  const auto initial_random_seed = given_seed > 0 ? given_seed : std::chrono::system_clock::now().time_since_epoch().count();
  int random_seed;
  std::mt19937 rng;
  for (int run_ind = 0; run_ind < num_bb_runs; run_ind++) {
    // Change the random seed per run
    random_seed = initial_random_seed * (run_ind + 1);
    params.set(intParam::RANDOM_SEED, random_seed);

    // For Cbc, for every run after the first, we will randomize the rows and columns of the input
    if (should_permute_rows_and_cols) {
      //std::srand(run_ind + 2);  // sets the seed for pseudo-randomness in C++
      std::seed_seq curr_seed{ random_seed };
      rng.seed(curr_seed);
      std::vector<int> this_row_permutation(row_permutation);
      std::vector<int> this_col_permutation(col_permutation);
      std::shuffle(this_row_permutation.begin(), this_row_permutation.end(), rng);
      std::shuffle(this_col_permutation.begin(), this_col_permutation.end(), rng);

      // Create the permuted problem
      std::vector<double> rowLB(solver->getNumRows()), rowUB(solver->getNumRows());
      std::vector<double> colLB(solver->getNumCols()), colUB(solver->getNumCols());
      std::vector<double> obj(solver->getNumCols());
      const CoinPackedMatrix* mat = solver->getMatrixByRow();

      // Set the rows
      CoinPackedMatrix row_mat;
      row_mat.reverseOrdering();  // make it row-ordered
      row_mat.setMinorDim(solver->getNumCols());
      for (int row = 0; row < solver->getNumRows(); row++) {
        const int curr_row = this_row_permutation[row];
        row_mat.appendRow(mat->getVector(curr_row));
        rowLB[row] = solver->getRowLower()[curr_row];
        rowUB[row] = solver->getRowUpper()[curr_row];
        printf("New row: %d. Orig row: %d. Row upper: %s. Row lower: %s.\n",
            row, curr_row, stringValue(rowLB[row]).c_str(),
            stringValue(rowUB[row]).c_str());
      }
      row_mat.reverseOrdering(); // make it col-ordered for use in the next part

      // Set the cols
      CoinPackedMatrix col_mat;
      col_mat.setMinorDim(solver->getNumRows());
      for (int col = 0; col < solver->getNumCols(); col++) {
        const int curr_col = this_col_permutation[col];
        col_mat.appendCol(row_mat.getVector(curr_col));
        colLB[col] = solver->getColLower()[curr_col];
        colUB[col] = solver->getColUpper()[curr_col];
        obj[col] = solver->getObjCoefficients()[curr_col];
      }

      runSolver = new SolverInterface;
      runSolver->setObjSense(solver->getObjSense());
      double objOffset = 0.;
      solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
      runSolver->setDblParam(OsiDblParam::OsiObjOffset, objOffset);

      runSolver->loadProblem(col_mat, colLB.data(), colUB.data(), obj.data(),
          rowLB.data(), rowUB.data());

      // Set integers
      for (int col = 0; col < solver->getNumCols(); col++) {
        const int curr_col = this_col_permutation[col];
        if (solver->isInteger(curr_col)) {
          runSolver->setInteger(curr_col);
        }
      }
    } // test whether num_bb_runs >= 2 and permute row/col order

    // Do branch and bound
    if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::gurobi)) {
#ifdef USE_GUROBI
      if (params.get(TEMP) == static_cast<int>(TempOptions::CHECK_CUTS_AGAINST_BB_OPT) && num_mycuts > 0) {
        if (ip_solution && !ip_solution->empty()) {
          // Get the original solution
          BBInfo tmp_bb_info;
          doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
              fullfilename.c_str(), tmp_bb_info, best_bound, ip_solution);
        }

        if (!ip_solution || ip_solution->empty()) {
          error_msg(errorstring, "Unable to obain IP solution.\n");
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }

        // Check cuts
        for (int cut_ind = 0; cut_ind < num_mycuts; cut_ind++) {
          OsiRowCut currCut = mycuts->rowCut(cut_ind);
          const double rhs = currCut.rhs();
          const int num_el = currCut.row().getNumElements();
          const int* ind = currCut.row().getIndices();
          const double* el = currCut.row().getElements();
          const double activity = dotProduct(num_el, ind, el, ip_solution->data());

          if (lessThanVal(activity, rhs)) {
            warning_msg(warnstring, "Cut %d removes optimal solution. Activity: %.10f. Rhs: %.10f.\n", cut_ind, activity, rhs);
          }
        }
      } // checking cuts for violating the IP opt
      if (branch_with_no_cuts) {
        // NB: If user cuts is not explicitly set, this is WITHOUT user cuts
        doBranchAndBoundWithGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), info_nocuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_mycuts) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), mycuts, info_mycuts->vec_bb_info[run_ind],
            best_bound);
      }
      if (branch_with_gmics) {
        doBranchAndBoundWithUserCutsGurobi(params, params.get(BB_STRATEGY),
            fullfilename.c_str(), &GMICsAndmycuts,
            info_allcuts->vec_bb_info[run_ind], best_bound);
      }
#endif // USE_GUROBI
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cplex)) {
//#ifdef USE_CPLEX
//      doBranchAndBoundWithUserCutsCplexCallable(
//          GlobalVariables::in_f_name.c_str(), &structVPC, vec_bb_info_mycuts[run_ind]);
//      doBranchAndBoundWithUserCutsCplexCallable(
//          GlobalVariables::in_f_name.c_str(), &SICsAndmycuts,
//          vec_bb_info_allcuts[run_ind]);
//#endif // USE_CPLEX
    } else if (use_bb_option(params.get(BB_STRATEGY), BB_Strategy_Options::cbc)) {
      if (runSolver == NULL) {
        runSolver = const_cast<OsiSolverInterface*>(solver);
      }
      if (branch_with_no_cuts) {
        doBranchAndBoundNoCuts(params, runSolver, info_nocuts->vec_bb_info[run_ind]);
      }
      if (branch_with_mycuts) {
        doBranchAndBoundYesCuts(params, runSolver, info_mycuts->vec_bb_info[run_ind],
            *mycuts, false, numCutsToAddPerRound, 1, "\nBB with mycuts.\n");
      }
      if (branch_with_gmics) {
        doBranchAndBoundYesCuts(params, runSolver, info_allcuts->vec_bb_info[run_ind],
            GMICsAndmycuts, false, numCutsToAddPerRound, 1, "\nBB with VPC+Gomorys.\n");
      }
    }

    if (branch_with_no_cuts) updateBestBBInfo(info_nocuts->best_bb_info, info_nocuts->vec_bb_info[run_ind], (run_ind == 0));
    if (branch_with_mycuts) updateBestBBInfo(info_mycuts->best_bb_info, info_mycuts->vec_bb_info[run_ind], (run_ind == 0));
    if (branch_with_gmics) updateBestBBInfo(info_allcuts->best_bb_info, info_allcuts->vec_bb_info[run_ind], (run_ind == 0));

    // Free memory if necessary
    if (should_permute_rows_and_cols && runSolver && (runSolver != solver)) {
      delete runSolver;
    }
  } /* end iterating over runs */

  if (branch_with_no_cuts) {
    info_nocuts->first_bb_info = info_nocuts->vec_bb_info[0];
    averageBBInfo(info_nocuts->avg_bb_info, info_nocuts->vec_bb_info);
  }
  if (branch_with_mycuts) {
    info_mycuts->first_bb_info = info_mycuts->vec_bb_info[0];
    averageBBInfo(info_mycuts->avg_bb_info, info_mycuts->vec_bb_info);
  }
  if (branch_with_gmics) {
    info_allcuts->first_bb_info = info_allcuts->vec_bb_info[0];
    averageBBInfo(info_allcuts->avg_bb_info, info_allcuts->vec_bb_info);
  }
} /* runBBTests */

/** Methods related to BBInfo */
void updateBestBBInfo(BBInfo& best_info, const BBInfo& curr_info, const bool first) {
  best_info.obj = first ? curr_info.obj : CoinMin(best_info.obj, curr_info.obj);
  best_info.bound = first ? curr_info.bound : CoinMax(best_info.bound, curr_info.bound);
  best_info.iters = first ? curr_info.iters : CoinMin(best_info.iters, curr_info.iters);
  best_info.nodes = first ? curr_info.nodes : CoinMin(best_info.nodes, curr_info.nodes);
  best_info.root_passes = first ? curr_info.root_passes : CoinMin(best_info.root_passes, curr_info.root_passes);
  best_info.first_cut_pass = first ? curr_info.first_cut_pass : CoinMax(best_info.first_cut_pass, curr_info.first_cut_pass);
  best_info.last_cut_pass = first ? curr_info.last_cut_pass : CoinMax(best_info.last_cut_pass, curr_info.last_cut_pass);
  best_info.root_iters = first ? curr_info.root_iters : CoinMin(best_info.root_iters, curr_info.root_iters);
  best_info.root_time = first ? curr_info.root_time : CoinMin(best_info.root_time, curr_info.root_time);
  best_info.last_sol_time = first ? curr_info.last_sol_time : CoinMin(best_info.last_sol_time, curr_info.last_sol_time);
  best_info.time = first ? curr_info.time : CoinMin(best_info.time, curr_info.time);
} /* updateMinBBInfo */

void averageBBInfo(BBInfo& avg_info, const std::vector<BBInfo>& info) {
  for (const BBInfo& curr_info : info) {
    avg_info.obj += curr_info.obj;
    avg_info.bound += curr_info.bound;
    avg_info.iters += curr_info.iters;
    avg_info.nodes += curr_info.nodes;
    avg_info.root_passes += curr_info.root_passes;
    avg_info.first_cut_pass += curr_info.first_cut_pass;
    avg_info.last_cut_pass += curr_info.last_cut_pass;
    avg_info.root_iters += curr_info.root_iters;
    avg_info.root_time += curr_info.root_time;
    avg_info.last_sol_time += curr_info.last_sol_time;
    avg_info.time += curr_info.time;
  }
  const int num_bb_runs = info.size();
  avg_info.obj /= num_bb_runs;
  avg_info.bound /= num_bb_runs;
  avg_info.iters /= num_bb_runs;
  avg_info.nodes /= num_bb_runs;
  avg_info.root_passes /= num_bb_runs;
  avg_info.first_cut_pass /= num_bb_runs;
  avg_info.last_cut_pass /= num_bb_runs;
  avg_info.root_iters /= num_bb_runs;
  avg_info.root_time /= num_bb_runs;
  avg_info.last_sol_time /= num_bb_runs;
  avg_info.time /= num_bb_runs;
} /* averageBBInfo */

void createStringFromBBInfoVec(const std::vector<BBInfo>& vec_info,
    std::vector<std::string>& vec_str) {
  vec_str.resize(BB_INFO_CONTENTS.size());
  for (const BBInfo& info : vec_info) {
    vec_str[OBJ_BB_INFO_IND] += (!vec_str[OBJ_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.obj, "%.20f");
    vec_str[BOUND_BB_INFO_IND] += (!vec_str[BOUND_BB_INFO_IND].empty() ? ";" : "") + stringValue(info.bound, "%.20f");
    vec_str[ITERS_BB_INFO_IND] += (!vec_str[ITERS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.iters);
    vec_str[NODES_BB_INFO_IND] += (!vec_str[NODES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.nodes);
    vec_str[ROOT_PASSES_BB_INFO_IND] += (!vec_str[ROOT_PASSES_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_passes);
    vec_str[FIRST_CUT_PASS_BB_INFO_IND] += (!vec_str[FIRST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.first_cut_pass);
    vec_str[LAST_CUT_PASS_BB_INFO_IND] += (!vec_str[LAST_CUT_PASS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_cut_pass);
    vec_str[ROOT_ITERS_BB_INFO_IND] += (!vec_str[ROOT_ITERS_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_iters);
    vec_str[ROOT_TIME_BB_INFO_IND] += (!vec_str[ROOT_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.root_time);
    vec_str[LAST_SOL_TIME_BB_INFO_IND] += (!vec_str[LAST_SOL_TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.last_sol_time);
    vec_str[TIME_BB_INFO_IND] += (!vec_str[TIME_BB_INFO_IND].empty() ? ";" : "") + std::to_string(info.time);
  }
} /* createStringFromBBInfoVec */

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(
    FILE* const logfile,
    const OsiSolverInterface* const solver,
    std::string& f_name) {
  // Generate temporary file name
  char template_name[] = "/tmp/tmpmpsXXXXXX";

  int fd;
  if ((fd = mkstemp(template_name)) == -1) {
    error_msg(errorstring, "Could not generate temp file with error %s.\n", strerror(errno));
    writeErrorToLog(errorstring, logfile);
    throw std::logic_error(errorstring);
  }

  f_name = template_name;
  if (f_name.empty()) {
    error_msg(errorstring, "Could not generate temp file.\n");
    writeErrorToLog(errorstring, logfile);
    throw std::logic_error(errorstring);
  }
  solver->writeMps(template_name, "mps", solver->getObjSense());
  f_name += ".mps.gz"; // writeMps calls writeMpsNative, which invokes the CoinMpsIO writer with gzip option = 1
  //doBranchAndBoundWithCplex(f_name.c_str(), bb_opt, bb_iters, bb_nodes, bb_time);
  //remove(f_name.c_str());
} /* createTmpFileCopy (Osi) */

#ifdef USE_CBC
void setStrategyForBBTestCbc(const Parameters& params,
    CbcModel* const cbc_model,
    int seed = -1) {
  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  // Parameters that should always be set
  cbc_model->setMaximumSeconds(params.get(doubleParam::BB_TIMELIMIT)); // time limit
  cbc_model->setRandomSeed(seed); // random seed

  int strategy = params.get(intParam::BB_STRATEGY);
  if (strategy <= 0) {
    // Default strategy
    CbcStrategyDefault strategy;
    cbc_model->setStrategy(strategy);
    /*
    CbcStrategyDefault strategy(-1);
    strategy.setupPreProcessing(-1,0);
    cbc_model->setStrategy(strategy);
    cbc_model->setMaximumCutPassesAtRoot(0);
    cbc_model->setMaximumCutPasses(0);
    cbc_model->setWhenCuts(0);
    */
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      // Make sure dual reductions are off
    }

    // Turn off all cuts
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      cbc_model->setMaximumCutPassesAtRoot(0);
      cbc_model->setMaximumCutPasses(0);
      cbc_model->setWhenCuts(0);
    }

    // Presolve
    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      cbc_model->setTypePresolve(0);
    }
  }

  // Check if we should use strong branching
  // Not sure this works when using StrategyDefault as well...
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    CbcBranchDefaultDecision branch;
    OsiChooseStrong choose;
    choose.setNumberStrong(cbc_model->solver()->getNumCols());
    choose.setNumberBeforeTrusted(std::numeric_limits<int>::max());
    branch.setChooseMethod(choose);
    cbc_model->setBranchingMethod(branch);
  }
} /* setStrategyForBBTestCbc */

/************************************************************/
/**
 * Perform branch-and-bound without cuts
 */
void doBranchAndBoundNoCuts(const Parameters& params,
    const OsiSolverInterface* const solver, BBInfo& info) {
#ifdef TRACE
  printf("\nBB with no cuts.\n");
#endif
  const bool test_using_main = false;

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver);

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
#ifndef CBC_VERSION_210PLUS
    CbcMain0(*cbc_model); // need to add CbcSolverUsefulData to call: https://www.coin-or.org/Doxygen/Cbc/classCbcSolverUsefulData.html
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundNoCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
#endif
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundNoCuts */

/************************************************************/
/**
 * Perform branch-and-bound using the given cuts, perhaps doing cut selection
 */
void doBranchAndBoundYesCuts(const Parameters& params,
    const OsiSolverInterface* const solver, BBInfo& info, const OsiCuts& structCuts,
    const bool doCutSelection, const int numCutsToAddPerRound,
    const int maxRounds, const std::string logstring) {
//#ifdef TRACE
  printf("%s", logstring.c_str());
//#endif
  const bool test_using_main = false;

  // Check that there are cuts
  const int numCuts = structCuts.sizeCuts();
  if (numCuts == 0) {
    info.nodes = 0;
    info.time = 0.;
    return;
  }

  // Set up solver
  SolverInterface* BBSolver;
  BBSolver = dynamic_cast<SolverInterface*>(solver->clone());
  setLPSolverParameters(BBSolver, params.get(VERBOSITY));

  // Apply cuts
//  if (doCutSelection) {
//    applyCutsInRounds(BBSolver, structCuts, numCutsToAddPerRound, maxRounds);
//  } else {
    applyCutsCustom(BBSolver, structCuts);
//  }

  // Set up model
  CbcModel* cbc_model = new CbcModel;
  cbc_model->swapSolver(BBSolver);
  cbc_model->setModelOwnsSolver(true);
  setIPSolverParameters(cbc_model, params.get(VERBOSITY));

  // Set up options and run B&B
  if (!test_using_main) {
    setStrategyForBBTestCbc(params, cbc_model);
#ifdef TRACE
    cbc_model->branchAndBound(3);
#else
    cbc_model->branchAndBound(0);
#endif
  } else {
#ifndef CBC_VERSION_210PLUS
    CbcMain0(*cbc_model); // Need to add CbcSolverUsefulData to call: https://www.coin-or.org/Doxygen/Cbc/classCbcSolverUsefulData.html
    std::string name, logLevel, presolveOnOff, preprocessOnOff, cutsOnOff, heurOnOff, solveOption;
    name = "BBHelper_doBranchAndBoundYesCuts";
    presolveOnOff = "-presolve=off";
    preprocessOnOff = "-preprocess=off";
    cutsOnOff = "-cuts=off";
    heurOnOff = "-heur=off";
    solveOption = "-solve";
#ifdef TRACE
    logLevel = "-loglevel=3";
#else
    logLevel = "-loglevel=0";
#endif

    int argc = 0;
    const char** cbc_options = new const char*[20];
    cbc_options[argc++] = name.c_str();
    cbc_options[argc++] = logLevel.c_str();
    cbc_options[argc++] = presolveOnOff.c_str();
    cbc_options[argc++] = preprocessOnOff.c_str();
    cbc_options[argc++] = cutsOnOff.c_str();
    cbc_options[argc++] = heurOnOff.c_str();
    cbc_options[argc++] = solveOption.c_str();

    CbcMain1(argc, cbc_options, *cbc_model);
    delete[] cbc_options;
#endif
  }

  // Collect statistics
  info.time = CoinCpuTime()
      - cbc_model->getDblParam(CbcModel::CbcStartSeconds);
  info.obj = cbc_model->getObjValue();
  info.iters = cbc_model->getIterationCount();
  info.nodes = cbc_model->getNodeCount();

  // Free
  if (cbc_model) {
    delete cbc_model;
    cbc_model = NULL;
  }
} /* doBranchAndBoundYesCuts */
#endif // USE_CBC
