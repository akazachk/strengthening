// Name:     CplexHelper.cpp
// Author:   A. M. Kazachkov
// Date:     2019-Mar-01
//-----------------------------------------------------------------------------
#include "GurobiHelper.hpp"

// CPLEX
#ifdef VPC_USE_CPLEX
// C interface
#include <ilcplex/cplexx.h>

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(CPXENVptr& env, CPXLPptr& lp, std::string& f_name) {
  if (f_name.empty()) {
    // Generate temporary file name
    char template_name[] = "/tmp/tmpmpsXXXXXX";

    mkstemp(template_name);
    f_name = template_name;
    if (f_name.empty()) {
      error_msg(errorstring, "Could not generate temp file.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    f_name += ".mps.gz";
  }
  CPXXwriteprob(env, lp, f_name.c_str(), NULL);
} /* createTmpFileCopy (Cplex) */

void setStrategyForBBTestCplexCallable(CPXENVptr& env,
    const int seed = params.get(intConst::RANDOM_SEED),
    const double best_bound = GlobalVariables::bestObjValue) {
  int status = 0;

  // Parameters that should always be set
  status += CPXXsetdblparam(env, CPXPARAM_TimeLimit, GlobalVariables::param.getBB_TIMELIMIT()); // time limit
  status += CPXXsetlongparam(env, CPXPARAM_Threads, 1); // single-threaded
  status += CPXXsetintparam(env, CPXPARAM_RandomSeed, seed); // random seed

#ifndef TRACE
  status += CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
#endif
#ifdef TRACE
  status += CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
  status += CPXXsetintparam(env, CPXPARAM_MIP_Interval, 1);
#endif

  int strategy = param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  if (strategy <= 0) {
    // Default strategy
  } else {
    if (strategy & BB_Strategy_Options::user_cuts) {
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL); // disable dynamic search
      status += CPXXsetlongparam(env, CPXPARAM_Preprocessing_Linear, 0); // disable this so that presolve does not discard user cuts during preprocessing
      status += CPXXsetlongparam(env, CPXPARAM_Preprocessing_Reduce, CPX_PREREDUCE_PRIMALONLY); // disable dual reductions
    }

    // Turn off all cuts
    if (strategy & BB_Strategy_Options::all_cuts_off) {
      //status += CPXXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0);
      //status += CPXXsetlongparam(env, CPXPARAM_MIP_Limits_CutPasses, -1);
      //status += CPXXsetlongparam(env, CPXPARAM_MIP_Limits_EachCutLimit, 0); // all but Gomory
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_BQP, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_Cliques, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_Covers, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_Disjunctive, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_FlowCovers, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_PathCut, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_Gomory, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_GUBCovers, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_Implied, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_LocalImplied, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_LiftProj, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_MIRCut, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_MCFCut, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_RLT, -1);
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Cuts_ZeroHalfCut, -1);
      //    status += CPXXsetlongparam(env, CPXPARAM_MIP_Limits_GomoryPass, 0);
      //    status += CPXXsetlongparam(env, CPXPARAM_MIP_Limits_GomoryCand, 0);
    }

    // Presolve
    if (strategy & BB_Strategy_Options::presolve_off) {
      status += CPXXsetlongparam(env, CPXPARAM_Preprocessing_Presolve, 0); // turn off presolve overall
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_PresolveNode, 0); // turn off presolve at nodes (not turned off by above?)
      status += CPXXsetlongparam(env, CPXPARAM_Preprocessing_Relax, 0); // turn off LP presolve at root
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_Probe, -1); // turn off probing
    }

    // Heuristics
    if (strategy & BB_Strategy_Options::heuristics_off) {
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_HeuristicFreq, -1); // turn off the "periodic" heuristic
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_FPHeur, -1); // turn off feasibility pump
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_LBHeur, 0); // local branching heuristic (default is off anyway)
      status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_RINSHeur, 0); // turn off RINS
    }

    // Feed the solver the best bound provided
    if (strategy & BB_Strategy_Options::use_best_bound) {
      if (!isInfinity(best_bound)) {
        status += CPXXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
      }
    }
  } /* else, strategy > 0 */

  // Check if we should use strong branching
  if (std::abs(strategy) & BB_Strategy_Options::strong_branching_on) {
    status += CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
    status += CPXXsetintparam(env, CPXPARAM_MIP_Limits_StrongCand, CPX_BIGINT); // smaller than max int
    //std::numeric_limits<int>::max());
    status += CPXXsetlongparam(env, CPXPARAM_MIP_Limits_StrongIt, CPX_BIGLONG); // smaller than max long
        //std::numeric_limits<int>::max());
    status += CPXXsetlongparam(env, CPXPARAM_MIP_Strategy_BBInterval, 1); // always select best bound node
  }

  // Check that there are no errors
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred when setting parameters. Status: %d.\n", status);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* setStrategyForBBTestCplexCallable */

void readFileIntoCplexCallable(const char* f_name, CPXENVptr& env, CPXLPptr& lp) {
  int status = 0;

  /* Initialize the CPLEX environment */
  env = CPXXopenCPLEX (&status);

  /* If an error occurs, the status value indicates the reason for
     failure.  A call to CPXXgeterrorstring will produce the text of
     the error message.  Note that CPXXopenCPLEX produces no output,
     so the only way to see the cause of the error is to use
     CPXXgeterrorstring.  For other CPLEX routines, the errors will
     be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */
  if ( env == NULL ) {
    char errmsg[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring (env, status, errmsg);
    error_msg(errorstring, "CPLEX (C): Could not open CPLEX environment. Error: %s", errmsg);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  /* Create the problem, using the filename as the problem name */
  lp = CPXXcreateprob (env, &status, f_name);

  /* A returned pointer of NULL may mean that not enough memory
     was available or there was some other problem.  In the case of
     failure, an error message will have been written to the error
     channel from inside CPLEX.  In this example, the setting of
     the parameter CPXXPARAM_ScreenOutput causes the error message to
     appear on stdout.  Note that most CPLEX routines return
     an error code to indicate the reason for failure.   */
  if ( lp == NULL ) {
    error_msg(errorstring, "CPLEX (C): Failed to create the LP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXXreadcopyprob (env, lp, f_name, NULL);
  if ( status ) {
    error_msg(errorstring, "CPLEX (C): Failed to read and copy the problem data.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* readFileIntoCplexCallable */

void presolveModelWithCplexCallable(CPXENVptr& env, CPXLPptr& lp, double& presolved_opt, std::string& presolved_name) {
#ifdef TRACE
  printf("\n## CPLEX (C): Presolving model ##\n");
#endif
  //warning_msg(warnstring, "CPLEX presolve saving not currently working.\n");
  //return; // Currently not functioning
  const int strategy = param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND, BB_Strategy_Options::presolve_on);
  setStrategyForBBTestCplexCallable(env);
  param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND, strategy);

  int status = 0;

  // Save presolved model
  if (presolved_name.empty()) {
    char template_name[] = "/tmp/tmpmpsXXXXXX"; // generate temporary file name
    mktemp(template_name);
    presolved_name = template_name;
  }
  if (presolved_name.empty()) {
    error_msg(errorstring, "Could not generate temp file.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  double objoff; // objective offset between the original and presolved models
  status = CPXXpreslvwrite(env, lp, (presolved_name + ".pre").c_str(), &objoff);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during call to CPXXpreslvwrite.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  CPXENVptr new_env = NULL;
  CPXLPptr new_lp = NULL;
  // Initialize the CPLEX environment
  new_env = CPXXopenCPLEX(&status);
  if ( new_env == NULL ) {
    char errmsg[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring (new_env, status, errmsg);
    error_msg(errorstring, "CPLEX (C): Could not open CPLEX environment. Error: %s", errmsg);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Create the problem, using the correct problem name
  new_lp = CPXXcreateprob(new_env, &status, GlobalVariables::prob_name.c_str());
  if ( new_lp == NULL ) {
    error_msg(errorstring, "CPLEX (C): Failed to create the LP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Now read the file, and copy the data into the created lp
  status = CPXXreadcopyprob(new_env, new_lp, (presolved_name + ".pre").c_str(), NULL);
  if (status) {
    char errmsg[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring (new_env, status, errmsg);
    error_msg(errorstring, "CPLEX (C): Error occurred during reading in of presolved file. Error: %s", errmsg);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Write the problem to an mps file format
  presolved_name += ".mps.gz";
  status = CPXXwriteprob(new_env, new_lp, presolved_name.c_str(), NULL);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during writing of presolved file to mps file.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Change to an lp to get the presolved lp optimum
  status = CPXXchgprobtype(new_env, new_lp, CPXPROB_LP);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during changing problem type.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Optimize the lp
  status = CPXXlpopt(new_env, new_lp);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during solving presolved lp.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Ensure that we are optimal
  int optimstatus = CPXXgetstat(new_env, new_lp);
  if (optimstatus != CPX_STAT_OPTIMAL) {
    error_msg(errorstring, "CPLEX (C): Error occurred during solving presolved lp.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Save the presolved objective value
  status = CPXXgetobjval(new_env, new_lp, &presolved_opt);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Unable to get presolved objective value.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Free up the problem as allocated by CPXXcreateprob, if necessary
  if ( new_lp != NULL ) {
    status = CPXXfreeprob (new_env, &new_lp);
    if ( status ) {
      error_msg(errorstring, "CPLEX (C): CPXXfreeprob failed, error code %d.\n", status);
      writeErrorToLog(errorstring, params.logfile);
    }
  }

  // Free up the CPLEX new_environment, if necessary
  if ( new_env != NULL ) {
    status = CPXXcloseCPLEX (&new_env);
    if ( status ) {
      char errmsg[CPXMESSAGEBUFSIZE];
      CPXXgeterrorstring (new_env, status, errmsg);
      error_msg(errorstring, "CPLEX (C): Could not close CPLEX environment. Error: %s", errmsg);
      writeErrorToLog(errorstring, params.logfile);
    }
  }
  /* CPXCLPptr presolved;
  CPXXgetredlp(env, lp, &presolved);
  if (presolved) {
    //CPXXpreslvwrite
   //CPXwriteprob(env, presolved, ...);
  } else {
    error_msg(errorstring, "CPLEX (C): Not able to get presolved model from CPLEX.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }*/
  /*status = CPXXpresolve(env, lp, CPX_ALG_NONE); // CPX_ALG_NONE should be used for MIP
    if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during presolve.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }*/

  // Convert variables to continuous
  /*
  const CPXDIM num_vars = CPXXgetnumcols(env, lp);
  CPXDIM* indices = new CPXDIM[num_vars];
  char* xctype = new char[num_vars];
  for (CPXDIM i = 0; i < num_vars; i++) {
    indices[i] = i;
    xctype[i] = CPX_CONTINUOUS;
  }
  status = CPXXchgctype(env, lp, num_vars, indices, xctype);
  if (status) {
    error_msg(errorstring, "CPLEX (C): Error occurred during changing variable type.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  delete indices;
  delete xctype;
  */
} /* presolveModelWithCplexCallable (CPXENVptr, CPXLPptr) */

void presolveModelWithCplexCallable(const char* f_name, double& presolved_opt, std::string& presolved_name) {
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;
  readFileIntoCplexCallable(f_name, env, lp);
  presolveModelWithCplexCallable(env, lp, presolved_opt, presolved_name);

  /* Free up the problem as allocated by CPXXcreateprob, if necessary */
  int status = 0;
  if ( lp != NULL ) {
    status = CPXXfreeprob (env, &lp);
    if ( status ) {
      error_msg(errorstring, "CPLEX (C): CPXXfreeprob failed, error code %d.\n", status);
      writeErrorToLog(errorstring, params.logfile);
    }
  }

  /* Free up the CPLEX environment, if necessary */
  if ( env != NULL ) {
    status = CPXXcloseCPLEX (&env);

    /* Note that CPXXcloseCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPXXPARAM_ScreenOutput indicator is set to CPXX_ON. */

      if ( status ) {
        char errmsg[CPXMESSAGEBUFSIZE];
        CPXXgeterrorstring (env, status, errmsg);
        error_msg(errorstring, "CPLEX (C): Could not close CPLEX environment. Error: %s", errmsg);
        writeErrorToLog(errorstring, params.logfile);
      }
   }
} /* presolveModelWithCplexCallable (filename) */

void presolveModelWithCplexCallable(const OsiSolverInterface* const solver, double& presolved_opt, std::string& presolved_name) {
  std::string f_name;
  createTmpFileCopy(solver, f_name);
  presolveModelWithCplexCallable(f_name.c_str(), presolved_opt, presolved_name);
  remove(f_name.c_str()); // remove temporary file
} /* presolveModelWithCplexCallable (Osi) */

void doBranchAndBoundWithCplexCallable(CPXENVptr& env, CPXLPptr& lp, BBInfo& info) {
//#ifdef TRACE
  printf("\n## Running B&B with CPLEX (Callable). Strategy: %d. ##\n", param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND));
//#endif
  int status = 0;

  // Set CPLEX parameters
  setStrategyForBBTestCplexCallable(env);

  // Optimize the LP
  //status = CPXXlpopt (env, lp);

  /* Optimize the problem and obtain solution. */
  // Exceeding a user-specified CPLEX limit is not considered an error.
  // Proving the problem infeasible or unbounded is not considered an error.
  double start_time = 0., end_time = 0.;
  CPXXgettime(env, &start_time);
  status = CPXXmipopt (env, lp);
  CPXXgettime(env, &end_time);
  if ( status ) {
    error_msg(errorstring, "CPLEX (C): Failed to optimize MIP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  int optimstatus = CPXXgetstat (env, lp);
  /*
  status = CPXXsolution(env, lp, &optimstatus, &bb_opt, NULL, NULL, NULL, NULL);
  if ( status ) {
    error_msg(errorstring, "CPLEX (C): Failed to get solution.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  char* optimstatusstring = new char[CPXMESSAGEBUFSIZE];
  CPXXgetstatstring (env, optimstatus, optimstatusstring);
  printf("Sol stat string: %s\n", optimstatusstring);
  */

  // Possibly re-run
  if (optimstatus == CPXMIP_INForUNBD) {
    int presolve_flag = 0;
    CPXXgetintparam(env, CPXPARAM_Preprocessing_Presolve, &presolve_flag);
    if (presolve_flag) {
      status = CPXXsetlongparam(env, CPXPARAM_Preprocessing_Presolve, 0); // turn off presolve overall
      if (status) {
        error_msg(errorstring, "CPLEX (C): Failed to set presolve parameter.\n");
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }

      status = CPXXmipopt(env, lp);
      if ( status ) {
        error_msg(errorstring, "CPLEX (C): Failed to optimize MIP.\n");
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }

      optimstatus = CPXXgetstat (env, lp);
    }
  }

  if (optimstatus == CPXMIP_INForUNBD || optimstatus == CPXMIP_INFEASIBLE || optimstatus == CPXMIP_UNBOUNDED) {
    error_msg(errorstring, "CPLEX (C): Failed to optimize MIP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Get best objective value
  switch (optimstatus) {
    case CPXMIP_DETTIME_LIM_FEAS: {
    }
    case CPXMIP_DETTIME_LIM_INFEAS: {
    }
    case CPXMIP_NODE_LIM_FEAS: {
    }
    case CPXMIP_NODE_LIM_INFEAS: {
    }
    case CPXMIP_TIME_LIM_FEAS: {
    }
    case CPXMIP_TIME_LIM_INFEAS: {
      status = CPXXgetbestobjval (env, lp, &info.obj);
      break;
    }
    case CPXMIP_OPTIMAL_TOL: {
    }
    case CPXMIP_OPTIMAL: {
      status = CPXXgetobjval (env, lp, &info.obj);
      break;
    }
    default: {
      char* optimstatusstring = new char[CPXMESSAGEBUFSIZE];
      CPXXgetstatstring (env, optimstatus, optimstatusstring);
      error_msg(errorstring, "CPLEX (C): Other status after solve: %d, message: %s.\n", optimstatus, optimstatusstring);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  } /* switch optimstatus */

  if ( status ) {
    error_msg(errorstring, "CPLEX (C): Failed to obtain objective value.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  info.iters = CPXXgetmipitcnt(env, lp);
  info.nodes = CPXXgetnodecnt(env, lp);
  info.time = end_time - start_time;

#ifdef TRACE
  printf("CPLEX (C): Solution value: %1.6f.\n", info.obj);
  printf("CPLEX (C): Number iterations: %ld.\n", info.iters);
  printf("CPLEX (C): Number nodes: %ld.\n", info.nodes);
  printf("CPLEX (C): Time: %f.\n", info.time);
#endif

  /* Free up the problem as allocated by CPXXcreateprob, if necessary */
  if ( lp != NULL ) {
    status = CPXXfreeprob (env, &lp);
    if ( status ) {
      error_msg(errorstring, "CPLEX (C): CPXXfreeprob failed, error code %d.\n", status);
      writeErrorToLog(errorstring, params.logfile);
    }
  }

  /* Free up the CPLEX environment, if necessary */
  if ( env != NULL ) {
    status = CPXXcloseCPLEX (&env);

    /* Note that CPXXcloseCPLEX produces no output,
       so the only way to see the cause of the error is to use
       CPXXgeterrorstring.  For other CPLEX routines, the errors will
       be seen if the CPXXPARAM_ScreenOutput indicator is set to CPXX_ON. */

      if ( status ) {
        char errmsg[CPXMESSAGEBUFSIZE];
        CPXXgeterrorstring (env, status, errmsg);
        error_msg(errorstring, "CPLEX (C): Could not close CPLEX environment. Error: %s\n", errmsg);
        writeErrorToLog(errorstring, params.logfile);
      }
   }
} /* doBranchAndBoundWithCplexCallable (CPXENVptr, CPXLPptr) */

void doBranchAndBoundWithUserCutsCplexCallable(CPXENVptr& env, CPXLPptr& lp,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy) {
  // Ensure that user cuts setting is enabled
  const int strategy = GlobalVariables::param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  if (!(strategy & BB_Strategy_Options::user_cuts)) {
    warning_msg(warnstring, "Need to use user_cuts option; strategy currently: %d.\n", strategy);
    GlobalVariables::param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND,
        strategy | BB_Strategy_Options::user_cuts);
  }

  // Add user cuts
  if (cuts) {
    const int num_cuts = cuts.sizeCuts();
    std::vector<CPXNNZ> cutbeg(num_cuts); // index where each cut starts
    std::vector<CPXDIM> cutind; // variables involved in each cut
    std::vector<double> cutval; // coefficients for the variabels in each cut
    std::vector<double> cutrhs(num_cuts); // rhs of each cut
    std::string cutsens = ""; // sense of each cut
    for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
      const OsiRowCut* curr_cut = cuts.rowCutPtr(cut_ind);
      const int num_el = curr_cut->row().getNumElements();
      const int* ind = curr_cut->row().getIndices();
      const double* vals = curr_cut->row().getElements();

      cutbeg[cut_ind] = cutind.size();
      cutind.insert(cutind.end(), ind, ind + num_el);
      cutval.insert(cutval.end(), vals, vals + num_el);
      cutrhs[cut_ind] = curr_cut->rhs();
      cutsens += "G";
    } /* end iterating over cuts */

    int status = CPXXaddusercuts(env, lp, cutbeg.size(), cutind.size(), cutrhs.data(), cutsens.c_str(), cutbeg.data(), cutind.data(), cutval.data(), NULL);
    if (status) {
      error_msg(errorstring, "CPLEX (C): Error adding user cuts. Error: %d.\n", status);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    if (addAsLazy) {
      status = CPXXaddlazyconstraints(env, lp, cutbeg.size(), cutind.size(), cutrhs.data(), cutsens.c_str(), cutbeg.data(), cutind.data(), cutval.data(), NULL);
      if (status) {
        error_msg(errorstring, "CPLEX (C): Error adding lazy cuts. Error: %d.\n", status);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    }
  } /* ensure cuts is not NULL */
  // Continue in normal routine
  doBranchAndBoundWithCplexCallable(env, lp, info); // does freeing

  GlobalVariables::param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND, strategy);
} /* doBranchAndBoundWithUserCutsCplexCallable (CPXenvptr) */

void doBranchAndBoundWithCplexCallable(const char* f_name, BBInfo& info) {
#ifdef TRACE
  printf("\n## Reading from file into CPLEX (C). ##\n");
#endif
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;
  readFileIntoCplexCallable(f_name, env, lp);
  const int strategy = GlobalVariables::param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  if (strategy & BB_Strategy_Options::user_cuts) {
    doBranchAndBoundWithUserCutsCplexCallable(env, lp, NULL, info, false); // add as lazy switch not available here
  } else {
    doBranchAndBoundWithCplexCallable(env, lp, info); // does freeing
  }
} /* doBranchAndBoundWithCplexCallable (file) */

void doBranchAndBoundWithCplexCallable(const OsiSolverInterface* const solver,
    BBInfo& info) {
  std::string f_name;
  createTmpFileCopy(solver, f_name);
  doBranchAndBoundWithCplexCallable(f_name.c_str(), info);
  remove(f_name.c_str());
} /* doBranchAndBoundWithCplexCallable (Osi) */

void doBranchAndBoundWithUserCutsCplexCallable(const char* f_name,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy) {
#ifdef TRACE
  printf("\n## Reading from file into CPLEX (C) and adding user cuts. ##\n");
#endif
  CPXENVptr env = NULL;
  CPXLPptr lp = NULL;
  int status = 0;

  /* Initialize the CPLEX environment */
  env = CPXXopenCPLEX (&status);

  /* If an error occurs, the status value indicates the reason for
     failure.  A call to CPXXgeterrorstring will produce the text of
     the error message.  Note that CPXXopenCPLEX produces no output,
     so the only way to see the cause of the error is to use
     CPXXgeterrorstring.  For other CPLEX routines, the errors will
     be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */
  if ( env == NULL ) {
    char errmsg[CPXMESSAGEBUFSIZE];
    CPXXgeterrorstring (env, status, errmsg);
    error_msg(errorstring, "CPLEX (C): Could not open CPLEX environment. Error: %s", errmsg);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  /* Create the problem, using the filename as the problem name */
  lp = CPXXcreateprob (env, &status, f_name);

  /* A returned pointer of NULL may mean that not enough memory
     was available or there was some other problem.  In the case of
     failure, an error message will have been written to the error
     channel from inside CPLEX.  In this example, the setting of
     the parameter CPXXPARAM_ScreenOutput causes the error message to
     appear on stdout.  Note that most CPLEX routines return
     an error code to indicate the reason for failure.   */
  if ( lp == NULL ) {
    error_msg(errorstring, "CPLEX (C): Failed to create the LP.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  /* Now read the file, and copy the data into the created lp */
  status = CPXXreadcopyprob (env, lp, f_name, NULL);
  if ( status ) {
    error_msg(errorstring, "CPLEX (C): Failed to read and copy the problem data.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  doBranchAndBoundWithUserCutsCplexCallable(env, lp, cuts, info, addAsLazy); // does freeing
} /* doBranchAndBoundWithUserCutsCplexCallable (filename) */

void doBranchAndBoundWithUserCutsCplexCallable(const OsiSolverInterface* const solver,
    const OsiCuts& cuts, BBInfo& info, const bool addAsLazy) {
  std::string f_name;
  createTmpFileCopy(solver, f_name);
  doBranchAndBoundWithUserCutsCplexCallable(f_name.c_str(), cuts, info, addAsLazy);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithUserCutsCplexCallable (Osi) */

// C++ interface
#ifdef VPC_USE_CPLEX_CONCERT
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

void setStrategyForBBTestCplexConcert(IloCplex& cplex,
    const int seed = params.get(intConst::RANDOM_SEED),
    const double best_bound = GlobalVariables::bestObjValue) {
  // Parameters that should always be set
  cplex.setParam(IloCplex::Param::TimeLimit, GlobalVariables::param.getBB_TIMELIMIT()); // time limit
  cplex.setParam(IloCplex::Param::Threads, 1); // single-threaded
  cplex.setParam(IloCplex::Param::RandomSeed, seed); // random seed

#ifndef TRACE
    cplex.setOut(cplex.getEnv().getNullStream());
#endif

  int strategy = param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  if (strategy <= 0) {
    // Default strategy
  } else {
    if (strategy & BB_Strategy_Options::user_cuts) {
      cplex.setParam(IloCplex::Param::MIP::Strategy::Search, CPX_MIPSEARCH_TRADITIONAL); // disable dynamic search
      cplex.setParam(IloCplex::Param::Preprocessing::Linear, 0); // presolved model vars linear comb of original
      cplex.setParam(IloCplex::Param::Preprocessing::Reduce, CPX_PREREDUCE_PRIMALONLY); // disable dual reductions
    }

    // Turn off cuts
    if (strategy & BB_Strategy_Options::all_cuts_off) {
      //cplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 0.);
      //    cplex.setParam(IloCplex::Param::MIP::Limits::CutPasses, -1);
      //cplex.setParam(IloCplex::Param::MIP::Limits::EachCutLimit, 0); // disable all but Gomory
      cplex.setParam(IloCplex::Param::MIP::Cuts::BQP, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::Cliques, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::Covers, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::Disjunctive, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::PathCut, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::GUBCovers, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::Implied, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::LocalImplied, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::LiftProj, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::MCFCut, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::RLT, -1);
      cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);
      //    cplex.setParam(IloCplex::Param::MIP::Limits::GomoryPass, 0);
      //    cplex.setParam(IloCplex::Param::MIP::Limits::GomoryCand, 0);
    }

    // Presolve
    if (strategy & BB_Strategy_Options::presolve_off) {
      cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0); // turn off presolve overall
      cplex.setParam(IloCplex::Param::MIP::Strategy::PresolveNode, 0); // turn off presolve at nodes
      cplex.setParam(IloCplex::Param::Preprocessing::Relax, 0); // turn off LP presolve at root
      cplex.setParam(IloCplex::Param::MIP::Strategy::Probe, -1); // turn off probing
    }

    // Heuristics
    if (strategy & BB_Strategy_Options::heuristics_off) {
      cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1); // turn off the "periodic" heuristic
      cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, -1); // turn off feasibility pump
      cplex.setParam(IloCplex::Param::MIP::Strategy::LBHeur, 0); // local branching heuristic (default is off anyway)
      cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, 0); // turn off RINS
    }

    // Feed the solver the best bound provided
    if (strategy & BB_Strategy_Options::use_best_bound) {
      if (!isInfinity(best_bound)) {
        cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
      }
    }
  }

  // Check if we should use strong branching
  if (std::abs(strategy) & BB_Strategy_Options::strong_branching_on) {
    cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect,
        IloCplex::VariableSelect::Strong);
    cplex.setParam(IloCplex::Param::MIP::Limits::StrongCand, CPX_BIGINT); // smaller than max int
        //std::numeric_limits<int>::max());
    cplex.setParam(IloCplex::Param::MIP::Limits::StrongIt, CPX_BIGLONG); // smaller than max long
        //std::numeric_limits<int>::max());
    cplex.setParam(IloCplex::Param::MIP::Strategy::BBInterval, 1); // always select best bound node
  }
} /* setStrategyForBBTestCplexConcert */

void doBranchAndBoundWithCplexConcert(IloEnv& env, IloModel& model, IloCplex& cplex,
    BBInfo& info) {
#ifdef TRACE
  printf("\n## Running B&B with CPLEX (Concert). Strategy: %d. ##\n", param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND));
#endif
  setStrategyForBBTestCplexConcert(cplex); // set parameters

  try {
    IloObjective obj;
    IloNumVarArray var(env);
    IloRangeArray rng(env);
    //model.add(IloConversion(env, var, ILOFLOAT));

    const IloNum start_time = cplex.getCplexTime();
    bool isFeasibleSolutionFound = cplex.solve();
    IloNum end_time = cplex.getCplexTime();

    IloAlgorithm::Status optimstatus = cplex.getStatus();
    if (optimstatus == IloAlgorithm::Status::InfeasibleOrUnbounded) {
      if (presolve_flag) {
        isFeasibleSolutionFound = cplex.solve();
        end_time = cplex.getCplexTime();
        optimstatus = cplex.getStatus();
      }
    }

    /*
    if (!isFeasibleSolutionFound) {
      error_msg(errorstring, "CPLEX (C): No feasible solution found using Cplex. Status: %d.\n", optimstatus);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    */

    if (optimstatus == IloAlgorithm::Status::Optimal) {
      info.obj = cplex.getObjValue();
    } else if (optimstatus == IloCplex::CplexStatus::AbortTimeLim) {
      // Abort due to time limit
      info.obj = cplex.getBestObjValue();
    } else if (optimstatus == IloCplex::CplexStatus::AbortItLim) {
      // Abort due to iteration limit
      info.obj = cplex.getBestObjValue();
    } else {
      error_msg(errorstring, "CPLEX (C): Other status after solve: %d.\n", optimstatus);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    info.iters = cplex.getNiterations();
    info.nodes = cplex.getNnodes();
    info.time = end_time - start_time;

#ifdef TRACE
    printf("CPLEX (Concert): Solution value: %e.\n", info.obj);
    printf("CPLEX (Concert): Number iterations: %ld.\n", info.iters);
    printf("CPLEX (Concert): Number nodes: %ld.\n", info.nodes);
    printf("CPLEX (Concert): Time: %f.\n", info.time);
#endif
    //IloNumArray vals(env);
    //cplex.getValues(vals, var);
    //printf("Solution status: %s.\n", cplex.getStatus());

    //env.out() << "Solution vector = " << vals << endl;

    /*
    try {     // basis may not exist
      IloCplex::BasisStatusArray cstat(env);
      cplex.getBasisStatuses(cstat, var);
      env.out() << "Basis statuses  = " << cstat << endl;
    } catch (...) {
    }
    */
    //env.out() << "Maximum bound violation = "
    //    << cplex.getQuality(IloCplex::MaxPrimalInfeas) << endl;
//    IloNumVarArray vars(env);
//    vars.add(IloNumVar(env, 0.0, 40.0));
//    vars.add(IloNumVar(env));
//    vars.add(IloNumVar(env));
//    model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
//    model.add(-vars[0] + vars[1] + vars[2] <= 20);
//    model.add(vars[0] - 3 * vars[1] + vars[2] <= 30);
//    IloCplex cplex(model);
//    if (!cplex.solve()) {
//      env.error() << "Failed to optimize LP." << std::endl;
//      throw(-1);
//    }
//    IloNumArray vals(env);
//    env.out() << "Solution status = " << cplex.getStatus() << std::endl;
//    env.out() << "Solution value = " << cplex.getObjValue() << std::endl;
//    cplex.getValues(vals, vars);
//    env.out() << "Values = " << vals << std::endl;
  } catch (IloException& e) {
    error_msg(errorstring, "CPLEX (C++): Exception caught: %s.\n", e.getMessage());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "CPLEX (C++): Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  env.end();
} /* doBranchAndBoundWithCplexConcert (IloEnv, IloModel, IloCplex) */

void doBranchAndBoundWithUserCutsCplexConcert(IloEnv& env, IloModel& model, IloCplex& cplex,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy) {
  // Ensure that user cuts setting is enabled
  const int strategy = GlobalVariables::param.getParamVal(ParamIndices::BB_STRATEGY_PARAM_IND);
  if (!(strategy & BB_Strategy_Options::user_cuts)) {
    warning_msg(warnstring,
        "Need to use user_cuts option; strategy currently: %d.\n", strategy);
    GlobalVariables::param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND,
        strategy | BB_Strategy_Options::user_cuts);
  }

#ifdef TRACE
  printf("\n## Reading from file into CPLEX (Concert) and adding user cuts. ##\n");
#endif
  try {
    if (cuts) {
      // Add user cuts
      IloRangeArray cplex_cuts(env);
      for (int cut_ind = 0; cut_ind < cuts.sizeCuts(); cut_ind++) {
        const OsiRowCut* curr_cut = cuts.rowCutPtr(cut_ind);
  //      const double* curr_coeff = curr_cut->row().denseVector(n);
        const int* ind = curr_cut->row().getIndices();
        const double* vals = curr_cut->row().getElements();

        // Set IloConstraint
        IloRange curr_cplex_cut(env, curr_cut->rhs(), IloInfinity, "ObjCut");
        for (int i = 0; i < curr_cut->row().getNumElements(); i++) {
          const IloNumVar curr_var = var[ind[i]];
          const IloNum curr_value = vals[i];
          curr_cplex_cut.setLinearCoef(curr_var, curr_value);
        }
        cplex_cuts.add(curr_cplex_cut);
      } /* finish setting cplex version of given cuts */
      cplex.addUserCuts(cplex_cuts);
      if (addAsLazy) {
        cplex.addLazyConstraints(cplex_cuts);
      }
      cplex_cuts.endElements();
      cplex_cuts.end();
    } /* check that cuts is not NULL */

    doBranchAndBoundWithCplexConcert(env, model, cplex, info); // does freeing
  } catch (IloException& e) {
    error_msg(errorstring, "CPLEX (C++): Exception caught: %s.\n", e.getMessage());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "CPLEX (C++): Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  GlobalVariables::param.setParamVal(ParamIndices::BB_STRATEGY_PARAM_IND, strategy);
} /* doBranchAndBoundWithUserCutsCplexConcert (IloEnv, IloModel, IloCplex) */

void doBranchAndBoundWithCplexConcert(const char* f_name, BBInfo& info) {
#ifdef TRACE
  printf("\n## Reading from file into CPLEX (Concert). ##\n");
#endif
  IloEnv env;
  try {
    IloModel model(env);
    IloCplex cplex(model);
    cplex.importModel(model, f_name, obj, var, rng);
    if (strategy & BB_Strategy_Options::user_cuts) {
      doBranchAndBoundWithUserCutsCplexConcert(env, model, cplex, NULL, info, false); // add as lazy is turned off here because no option is passed
    } else {
      doBranchAndBoundWithCplexConcert(env, model, cplex, info); // does freeing
    }
  } catch (IloException& e) {
    error_msg(errorstring, "CPLEX (C++): Exception caught: %s.\n", e.getMessage());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "CPLEX (C++): Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithCplexConcert (filename) */

void doBranchAndBoundWithCplexConcert(const OsiSolverInterface* const solver, BBInfo& info) {
  std::string f_name;
  createTmpFileCopy(solver, f_name);
  doBranchAndBoundWithCplexConcert(f_name.c_str(), info);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithCplexConcert (Osi, creates temporary file from solver) */

void doBranchAndBoundWithUserCutsCplexConcert(const char* f_name,
    const OsiCuts* cuts, BBInfo& info, const bool addAsLazy) {
  IloEnv env;
  try {
    IloModel model(env);
    IloCplex cplex(model);
    cplex.importModel(model, f_name, obj, var, rng);

    doBranchAndBoundWithUserCutsCplexConcert(env, model, cplex, cuts, info, addAsLazy);
  } catch (IloException& e) {
    error_msg(errorstring, "CPLEX (C++): Exception caught: %s.\n", e.getMessage());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "CPLEX (C++): Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithUserCutsCplexConcert (f_name) */

void doBranchAndBoundWithUserCutsCplexConcert(const OsiSolverInterface* const solver,
    const OsiCuts& cuts, BBInfo& info, const bool addAsLazy) {
  std::string f_name;
  createTmpFileCopy(solver, f_name);
  doBranchAndBoundWithUserCutsCplexConcert(f_name.c_str(), cuts, info, addAsLazy);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithUserCutsCplexConcert (Osi) */
#endif /* VPC_USE_CPLEX_CONCERT */
#endif /* VPC_USE_CPLEX */
