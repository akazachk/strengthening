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
#include "analysis.hpp"
#include "BBHelper.hpp"
#include "CglAdvCut.hpp"
#include "CutHelper.hpp"
#include "Disjunction.hpp" // DisjunctiveTerm, Disjunction
//#include "disjcuts.hpp"
#include "gmic.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "strengthen.hpp"
#include "utility.hpp"

enum OverallTimeStats {
  INIT_SOLVE_TIME,
  CUT_TIME,
  BB_TIME,
  APPLY_TIME,
  TOTAL_TIME,
  NUM_TIME_STATS
}; /* OverallTimeStats */
const std::vector<std::string> OverallTimeStatsName {
  "INIT_SOLVE_TIME",
  "CUT_TIME",
  "BB_TIME",
  "APPLY_TIME",
  "TOTAL_TIME"
};

// Main file variables
Parameters params;
OsiSolverInterface *solver, *origSolver;
OsiSolverInterface* GMICSolver = NULL;
OsiSolverInterface* CutSolver = NULL;
OsiCuts gmics, mycuts;
std::string dir = "", filename = "", instname = "", in_file_ext = "";
CglAdvCut::ExitReason exitReason;
TimeStats timer;
std::time_t start_time_t, end_time_t;
char start_time_string[25];

SummaryBoundInfo boundInfo;
SummaryBBInfo info_nocuts, info_mycuts, info_allcuts;
std::vector<SummaryCutInfo> cutInfoVec;
SummaryCutInfo cutInfo, cutInfoGMICs;

// For output
std::string cut_output = "", bb_output = "";

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

void startUp(int argc, char** argv);
void processArgs(int argc, char** argv);
void initializeSolver(OsiSolverInterface* &solver);
int wrapUp(int retCode);

/****************** MAIN FUNCTION **********************/
int main(int argc, char** argv) {
  // Do this early in your program's initialization
  std::signal(SIGABRT, signal_handler_with_error_msg);
  std::signal(SIGSEGV, signal_handler_with_error_msg);

  // Set up timing
  for (int t = 0; t < OverallTimeStats::NUM_TIME_STATS; t++) {
    timer.register_name(OverallTimeStatsName[t]);
  }

  // Print welcome message, set up logfile
  timer.start_timer(OverallTimeStats::TOTAL_TIME);
  startUp(argc, argv);

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

  // Save original solver in case we wish to come back to it later
  origSolver = solver->clone();
  if (!origSolver->isProvenOptimal()) {
    origSolver->initialSolve();
    checkSolverOptimality(origSolver, false);
  }

  // Also save copies for calculating other objective values
  // We only need these if cuts other than mycuts are generated
  // solver      ::  stores the cuts we actually want to "count"
  // GMICSolver  ::  only GMICs
  // CutSolver   ::  if GMICs count, this is only mycuts; otherwise, it is both GMICs and mycuts
  if (params.get(GOMORY) != 0) {
    GMICSolver = solver->clone();
    CutSolver = solver->clone();
  }

  // Now do rounds of cuts, until a limit is reached (e.g., time, number failures, number cuts, or all rounds are exhausted)
  boundInfo.num_mycuts = 0, boundInfo.num_gmic = 0;
  int num_rounds = params.get(ROUNDS);
  std::vector<OsiCuts> mycuts_by_round(num_rounds);
  cutInfoVec.resize(num_rounds);
  int round_ind = 0;
  for (round_ind = 0; round_ind < num_rounds; ++round_ind) {
    if (num_rounds > 1) {
      printf("\n## Starting round %d/%d. ##\n", round_ind+1, num_rounds);
    }
    timer.start_timer(OverallTimeStats::CUT_TIME);

    // Generate Gomory cuts
    // > 0: apply cuts to solver; < 0: apply cuts to CutSolver only
    // Option 1: GglGMI
    // Option 2: custom generate intersection cuts, calculate Farkas certificate, do strengthening
    if (std::abs(params.get(GOMORY)) == 1) {
      OsiCuts currGMICs;
      CglGMI GMIGen;
      GMIGen.getParam().setMAX_SUPPORT(solver->getNumCols());
      GMIGen.getParam().setMAX_SUPPORT_REL(0.5);
      GMIGen.generateCuts(*solver, currGMICs);
      gmics.insert(currGMICs);
      boundInfo.num_gmic += currGMICs.sizeCuts();
      applyCutsCustom(GMICSolver, currGMICs);
      boundInfo.gmic_obj = GMICSolver->getObjValue();
      if (params.get(GOMORY) == 1) {
        applyCutsCustom(solver, currGMICs);
      }
      if (params.get(GOMORY) == -1) {
        applyCutsCustom(CutSolver, currGMICs);
      }
    } // generate GMICs via CglGMIC

    if (std::abs(params.get(GOMORY)) == 2) {
      solver->enableFactorization();
      std::vector<int> varBasicInRow(solver->getNumRows(), -1);
      solver->getBasics(&varBasicInRow[0]);
      solver->disableFactorization();
      OsiCuts currGMICs;
      for (int var = 0; var < solver->getNumCols(); var++) {
        if (!solver->isInteger(var)) {
          continue;
        }

        const double val = solver->getColSolution()[var];
        if (isVal(val, std::floor(val)) || isVal(val, std::ceil(val))) {
          continue;
        }

        // Find row in which this variable is basic
        int splitVarRowIndex = -1;
        for (int row = 0; row < solver->getNumRows(); row++) {
          if (varBasicInRow[row] == var) {
            splitVarRowIndex = row;
            break;
          }
        }
        if (splitVarRowIndex == -1) {
          error_msg(errorstring, "Unable to find fractional variable %d in basis.\n", var);
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }

        // We generate the strengthened intersection cut (to compare against)
        OsiRowCut currGMIC;
        createMIG(currGMIC, solver, var, splitVarRowIndex, true);

        // Now we generate the unstrengthened intersection cut and see if we can strengthen it
        // The actual Farkas certificate for this cut can be derived in closed form from the basis
        OsiRowCut intCut;
        createMIG(intCut, solver, var, splitVarRowIndex, false);

        // For each term of the disjunction, 
        // we need to explicitly add the constraint(s) defining the disjunctive term
        const CoinPackedVector lhs = intCut.row();

        std::vector<double> v0, v1; // this will be of dimension rows + cols
        { // Check first side of the split
          // Calculate the certificate
          OsiSolverInterface* solver0 = solver->clone();
          const double el = -1.;
          solver0->addRow(1, &var, &el, -1. * std::floor(val), solver->getInfinity());

          getCertificate(v0, lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), solver0);

          if (solver0) delete solver0;
        } // check first side of the split

        { // Check second side of the split
          // Calculate the certificate
          OsiSolverInterface* solver1 = solver->clone();
          const double el = 1.;
          solver1->addRow(1, &var, &el, std::ceil(val), solver->getInfinity());

          getCertificate(v1, lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), solver1);

          if (solver1) delete solver1;
        } // check second side of the split

        // Double check that calculated value of Farkas certificate matches theoretical value
        const double delta0 = val - std::floor(val);
        const double delta1 = std::ceil(val) - val;
        if (!isVal(v0[solver->getNumRows()], 1.0 / delta0)) {
          error_msg(errorstring, "Value of disjunctive term multiplier does not match theory. u0: %.6g. theory: %.6g.\n", v0[solver->getNumRows()], 1.0 / delta0);
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }
        if (!isVal(v1[solver->getNumRows()], 1.0 / delta1)) {
          error_msg(errorstring, "Value of disjunctive term multiplier does not match theory. v0: %.6g. theory: %.6g.\n", v1[solver->getNumRows()], 1.0 / delta1);
          writeErrorToLog(errorstring, params.logfile);
          exit(1);
        }
      } // iterate over cols, generating GMICs
      gmics.insert(currGMICs);

      boundInfo.num_gmic += currGMICs.sizeCuts();
      applyCutsCustom(GMICSolver, currGMICs);
      boundInfo.gmic_obj = GMICSolver->getObjValue();
      if (params.get(GOMORY) == 2) {
        applyCutsCustom(solver, currGMICs);
      }
      if (params.get(GOMORY) == -2) {
        applyCutsCustom(CutSolver, currGMICs);
      }
    } // generate GMICs via gmic.hpp

    //==================================================//
    // Now for more general cuts
    // User can replace generateCuts method in CglAdvCut with whatever method is desired
    CglAdvCut gen(params);

    // Store the initial solve time in order to set a baseline for the PRLP resolve time
    gen.timer.add_value(
        CglAdvCut::CutTimeStatsName[static_cast<int>(CglAdvCut::CutTimeStats::INIT_SOLVE_TIME)],
        timer.get_value(OverallTimeStats::INIT_SOLVE_TIME));

    // Generate disjunctive cuts
    gen.generateCuts(*solver, mycuts_by_round[round_ind]); 
    Disjunction* disj = gen.gen.disj();
    exitReason = gen.exitReason;
    updateCutInfo(cutInfoVec[round_ind], &gen);
    boundInfo.num_mycuts += gen.num_cuts;

    // Get Farkas certificate
    if (disj && mycuts_by_round[round_ind].sizeCuts() > 0) {
      for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
        OsiSolverInterface* termSolver = disj->getSolverForTerm(term_ind, solver, params.get(StrengtheningParameters::doubleConst::DIFFEPS), params.logfile);
        if (!termSolver) {
          printf("Disjunctive term %d/%d not created successfully.\n", term_ind+1, disj->num_terms);
          continue;
        }
        for (int cut_ind = 0; cut_ind < mycuts_by_round[round_ind].sizeCuts(); cut_ind++) {
          OsiRowCut* disjCut = mycuts_by_round[round_ind].rowCutPtr(cut_ind);

          // For each term of the disjunction,
          // we need to explicitly add the constraint(s) defining the disjunctive term
          const CoinPackedVector lhs = disjCut->row();

          // Get dense cut coefficients so that we can set the objective vector
          std::vector<double> cut_coeff(solver->getNumCols(), 0.0);
          const int num_el = lhs.getNumElements();
          for (int i = 0; i < num_el; i++) {
            cut_coeff[lhs.getIndices()[i]] = lhs.getElements()[i];
          }

          std::vector<double> v; // this will be of dimension rows + cols
          getCertificate(v, lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), termSolver);
        } // loop over cuts
        if (termSolver) { delete termSolver; }
      } // loop over disjunctive terms
    } // check that disj exists and cuts were generated

    timer.end_timer(OverallTimeStats::CUT_TIME);

    timer.start_timer(OverallTimeStats::APPLY_TIME);
    applyCutsCustom(solver, mycuts_by_round[round_ind]);
    if (params.get(GOMORY) > 0) {
      boundInfo.gmic_mycuts_obj = solver->getObjValue();
      applyCutsCustom(CutSolver, mycuts_by_round[round_ind]);
      boundInfo.mycuts_obj = CutSolver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.gmic_mycuts_obj;
    }
    else if (params.get(GOMORY) < 0) {
      boundInfo.mycuts_obj = solver->getObjValue();
      applyCutsCustom(CutSolver, mycuts_by_round[round_ind]);
      boundInfo.gmic_mycuts_obj = CutSolver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.gmic_mycuts_obj;
    } else {
      boundInfo.mycuts_obj = solver->getObjValue();
      boundInfo.all_cuts_obj = boundInfo.mycuts_obj;
    }
    timer.end_timer(OverallTimeStats::APPLY_TIME);

    mycuts.insert(mycuts_by_round[round_ind]);

    printf(
        "\n## Round %d/%d: Completed round of cut generation (exit reason: %s). # cuts generated = %d.\n",
        round_ind + 1, params.get(ROUNDS),
        CglAdvCut::ExitReasonName[static_cast<int>(exitReason)].c_str(),
        mycuts_by_round[round_ind].sizeCuts());
    fflush(stdout);
    printf("Initial obj value: %1.6f. New obj value: %s. Disj lb: %s. ##\n",
        boundInfo.lp_obj, stringValue(solver->getObjValue(), "%1.6f").c_str(),
        stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());

    // Exit early from rounds of cuts if no cuts generated or solver is not optimal
    if (gen.num_cuts == 0 || !solver->isProvenOptimal()
        || isInfinity(std::abs(solver->getObjValue())))
      break;
  } // loop over rounds of cuts
  if (round_ind < num_rounds)
    num_rounds = round_ind+1;

  // Do branch-and-bound experiments (if requested)
  if (params.get(BB_RUNS) != 0) {
    // Collect cuts from all rounds
    timer.start_timer(BB_TIME);
    runBBTests(params, &info_nocuts, &info_mycuts, &info_allcuts,
        params.get(stringParam::FILENAME), solver, boundInfo.ip_obj, &mycuts, NULL);
    timer.end_timer(BB_TIME);
  }

  printf(
      "\n## Finished cut generation with %d cuts. Initial obj value: %s. Final obj value: %s. Disj lb: %s. ##\n",
      boundInfo.num_mycuts,
      stringValue(boundInfo.lp_obj, "%1.6f").c_str(),
      stringValue(boundInfo.mycuts_obj, "%1.6f").c_str(),
      stringValue(boundInfo.best_disj_obj, "%1.6f").c_str());
  timer.end_timer(OverallTimeStats::TOTAL_TIME);

#ifdef PRINT_LP_WITH_CUTS
  if (boundInfo.num_mycuts > 0) {
    std::string fileWithCuts = filename + "_cuts";
    solver->writeMps(fileWithCuts.c_str());
  }
#endif

  // Do analyses in preparation for printing
  setCutInfo(cutInfo, num_rounds, cutInfoVec.data());
  analyzeStrength(params, solver, cutInfoGMICs, cutInfo, &gmics, &mycuts,
      boundInfo, cut_output);
  analyzeBB(params, info_nocuts, info_mycuts, info_allcuts, bb_output);
  return wrapUp(0);
} /* main */

/**
 * Call this early to print welcome message, etc.
 */
void startUp(int argc, char** argv) {
  // Input handling
  printf("## Custom Cut Generator ##\n");
  printf("Aleksandr M. Kazachkov\n");
  for (int i = 0; i < argc; i++) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  time(&start_time_t);
  struct tm* start_timeinfo = localtime(&start_time_t);
  snprintf(start_time_string, sizeof(start_time_string) / sizeof(char), "%s", asctime(start_timeinfo));
  printf("Start time: %s\n", start_time_string);

  processArgs(argc, argv);

  // Get instance file
  printf("Instance file: %s\n", params.get(stringParam::FILENAME).c_str());
  
  parseFilename(dir, instname, in_file_ext, params.get(stringParam::FILENAME), params.logfile);
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
#ifdef TRACE
    std::cout << "Reading objective information from \"" + params.get(stringParam::OPTFILE) + "\"" << std::endl;
#endif
    boundInfo.ip_obj = getObjValueFromFile(params.get(stringParam::OPTFILE), params.get(stringParam::FILENAME), params.logfile);
    params.set(doubleParam::IP_OBJ, boundInfo.ip_obj);
#ifdef TRACE
    std::cout << "Best known objective value is " << boundInfo.ip_obj << std::endl;
#endif
    if (isInfinity(boundInfo.ip_obj)) {
      warning_msg(warnstring, "Did not find objective value.\n");
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
} /* startUp */

/**
 * Close the logfile and print to it
 */
int wrapUp(int retCode /*= 0*/) {
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
      printSummaryBBInfo({info_nocuts, info_mycuts}, params.logfile);
      // Orig prob
      printOrigProbInfo(origSolver, params.logfile);
      // Post-cut prob
      printPostCutProbInfo(solver, cutInfoGMICs, cutInfo, params.logfile);
      // Cut, obj, fail info
      printCutInfo(cutInfoGMICs, cutInfo, params.logfile);
      // Full B&B info
      printFullBBInfo({info_nocuts, info_mycuts}, params.logfile);
      // Print time info
      timer.print(params.logfile, 2); // only values
      // Print exit reason and finish
      fprintf(logfile, "%s,", CglAdvCut::ExitReasonName[exitReasonInt].c_str());
    }

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
  printf("Instance: %s\n", instname.c_str());
  printf("Start time: %s\n", start_time_string);
  printf("End time: %s\n", end_time_string);
  printf("Elapsed time: %.f seconds\n", difftime(end_time_t, start_time_t));

  if (solver) {
    delete solver;
  }
  if (origSolver) {
    delete origSolver;
  }
  if (GMICSolver) {
    delete GMICSolver;
  }
  if (CutSolver) {
    delete CutSolver;
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
    double objOffset = 0.;
    solver->getDblParam(OsiDblParam::OsiObjOffset, objOffset);
    if (objOffset != 0.) {
      solver->setDblParam(OsiDblParam::OsiObjOffset, -1. * objOffset);
    }
  }
} /* initializeSolver */

/**
 * See params.hpp for descriptions of the parameters
 */
void processArgs(int argc, char** argv) {
  // Handle inputs
  // struct option declared in getopt.h
  // name: name of the long option
  // has_arg: 0,1,2 for none, required, or optional
  // *flag: how results are returned; if NULL, getopt_long() returns val (e.g., can be the equivalent short option character), and o/w getopt_long() returns 0, and flag points to a var which is set to val if the option is found, but left unchanged if the option is not found
  // val: value to return, or to load into the variable pointed to by flag
  const char* const short_opts = "b:B:c:d:f:g:hi:l:m:o:r:s:t:v:";
  const struct option long_opts[] =
  {
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
      {"strengthen",            required_argument, 0, 's'},
      {"temp",                  required_argument, 0, 't'*'1'},
      {"timelimit",             required_argument, 0, 't'},
      {"verbosity",             required_argument, 0, 'v'},
      {nullptr, no_argument, nullptr, 0}
  };

  int inp;
  while ((inp = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
    switch (inp) {
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
                helpstring += "-v level, --verbosity=level\n\tVerbosity level (0: print little, 1: let solver output be visible).\n";
                helpstring += "\n# General cut options #\n";
                helpstring += "-c num cuts, --cutlimit=num cuts\n\tMaximum number of cuts to generate (0+ = as given, -k = k * # fractional variables at root).\n";
                helpstring += "-d num terms, --disj_terms=num terms\n\tMaximum number of disjunctive terms or disjunctions to generate (depending on mode).\n";
                helpstring += "-g -1/0/1, --gomory=-1/0/1\n\t0: do not use Gomory cuts before generating mycuts, +/-1: generate Gomory cuts before generating mycuts (-1: only gen, +1: also apply to LP).\n";
                helpstring += "-m mode, --mode=mode\n\tDescription needs to be entered.\n";
                helpstring += "-r num rounds, --rounds=num rounds\n\tNumber of rounds of cuts to apply.\n";
                helpstring += "-s 0/1/2, --strengthen=0/1/2\n\tWhether to strengthen cuts.\n";
                helpstring += "-t num seconds, --timelimit=num seconds\n\tTotal number of seconds allotted for cut generation.\n";
                helpstring += "\n# Branch-and-bound options #\n";
                helpstring += "-b 0+ --bb_runs=0+\n\tNumber of branch-and-bound repeats.\n";
                helpstring += "-B strategy --bb_strategy=strategy\n\tBranch-and-bound strategy (see BBHelper.hpp).\n";
                helpstring += "--bb_mode={0,1,10,11,100,...,111}\n\tWhich branch-and-bound experiments to run (ones = no cuts, tens = mycuts, hundreds = gmics).\n";
                helpstring += "## END OF HELP ##\n";
                std::cout << helpstring << std::endl;
                exit(1);
               }
    } // switch statement for input
  } // process args
  //for (int i = 0; i < argc; i++) {
  //  std::cout << argv[i] << " ";
  //}
  //std::cout << std::endl;
} /* processArgs */
