/**
 * @file GurobiHelper.cpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-28
 */
#include "GurobiHelper.hpp"

//#include <cstdio> // for tmpnam

// Project files
#include "BBHelper.hpp"
#include "CutHelper.hpp" // applyCuts
#include "SolverHelper.hpp"
#include "Parameters.hpp"
#include "utility.hpp" // createTmpFilename
using namespace StrengtheningParameters;

// COIN-OR
#include <CoinTime.hpp>
#include <OsiCuts.hpp>
#include <CoinPackedMatrix.hpp>
#include <OsiSolverInterface.hpp>

// Gurobi
#ifdef USE_GUROBI
#include <gurobi_c++.h>

/**
 * Creates temporary file (in /tmp) so that it can be read by a different solver
 * It does not delete the file
 */
void createTmpFileCopy(FILE* const logfile, GRBModel& model, std::string& f_name, const std::string add_ext = ".mps.gz") {
  if (f_name.empty()) {
    try {
      createTmpFilename(f_name, add_ext);
    } catch (const std::exception &e) {
      error_msg(errorstring, "Could not generate temp file: %s.\n", e.what());
      writeErrorToLog(errorstring, logfile);
      exit(1);
    }
  } else {
    // Ensure f_name has proper extension
    if (!add_ext.empty() && f_name.compare(f_name.size()-add_ext.size(), add_ext.size(), add_ext) != 0) {
      f_name += add_ext;
    }
  }
  model.write(f_name.c_str());
} /* createTmpFileCopy (Gurobi) */

GRBModel* buildGRBModelFromOsi(const OsiSolverInterface* const solver, FILE* const logfile) {
  // Below DOES NOT WORK because when you make a copy of the file,
  // order of the variables may not be preserved
  // (Could fix this with ensuring variable names are checked, but also some precision may get lost in writing to a file)
  // std::string f_name;
  // createTmpFileCopy(logfile, solver, f_name);
  GRBEnv env = GRBEnv();
  GRBModel* model = new GRBModel(env);
  // remove(f_name.c_str()); // remove temporary file

  // Create variables
  for (int col = 0; col < solver->getNumCols(); col++) {
    const double lb = solver->getColLower()[col];
    const double ub = solver->getColUpper()[col];
    const double obj = solver->getObjCoefficients()[col];
    const std::string name = solver->getColName(col);
    const bool is_int = solver->isInteger(col);
    const bool is_bin = is_int && (lb == 0 && ub == 1);
    if (is_bin) {
      model->addVar(lb, ub, obj, GRB_BINARY, name);
    } else if (is_int) {
      model->addVar(lb, ub, obj, GRB_INTEGER, name);
    } else {
      model->addVar(lb, ub, obj, GRB_CONTINUOUS, name);
    }
  }
  model->update();

  // Copy constraints
  GRBVar* var = model->getVars();
  const CoinPackedMatrix* matrix = solver->getMatrixByRow();
  for (int row = 0; row < solver->getNumRows(); row++) {
    const double lb = solver->getRowLower()[row];
    const double ub = solver->getRowUpper()[row];
    const std::string name = solver->getRowName(row);
    GRBLinExpr expr;
    for (CoinBigIndex j = matrix->getVectorStarts()[row]; j < matrix->getVectorStarts()[row] + matrix->getVectorLengths()[row]; j++) {
      const int col = matrix->getIndices()[j];
      const double coeff = matrix->getElements()[j];
      expr += coeff * var[col];
    }
    if (lb == ub) {
      model->addConstr(expr, GRB_EQUAL, lb, name);
    } else if (isInfinity(std::abs(ub)) && !isInfinity(std::abs(lb))) {
      model->addConstr(expr, GRB_GREATER_EQUAL, lb, name);
    } else if (isInfinity(std::abs(lb)) && !isInfinity(std::abs(ub))) {
      model->addConstr(expr, GRB_LESS_EQUAL, ub, name);
    } else if (!isInfinity(std::abs(lb)) && !isInfinity(std::abs(ub))) {
      model->addConstr(expr, GRB_LESS_EQUAL, ub, name + "_ub");
      model->addConstr(expr, GRB_GREATER_EQUAL, lb, name + "_lb");
    } else {
      error_msg(errorstring, "buildGRBModelFromOsi: Constraint %s has no finite bounds.\n", name.c_str());
      writeErrorToLog(errorstring, logfile);
      exit(1);
    }
  }
  model->update();

  if (var) { delete[] var; }
  return model;
} /* buildGRBModelFromOsi */

/// @details This function is used to solve a model with Gurobi
///
/// @return GRB_IntAttr_Status value
///
/// The possible return statuses are:
/// \li GRB_LOADED = 1: Model is loaded, but no solution information is available.
/// \li GRB_OPTIMAL: Model was solved to optimality (subject to tolerances), and an optimal solution is available.
/// \li GRB_INFEASIBLE: Model was proven to be infeasible.
/// \li GRB_INF_OR_UNBD: Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.
/// \li GRB_UNBOUNDED: Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.
/// \li GRB_CUTOFF: Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.
/// \li GRB_ITERATION_LIMIT: Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.
/// \li GRB_NODE_LIMIT: Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.
/// \li GRB_TIME_LIMIT: Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.
/// \li GRB_SOLUTION_LIMIT: Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.
/// \li GRB_INTERRUPTED: Optimization was terminated by the user.
/// \li GRB_NUMERIC: Optimization was terminated due to unrecoverable numerical difficulties.
/// \li GRB_SUBOPTIMAL: Unable to satisfy optimality tolerances; a sub-optimal solution is available.
/// \li GRB_INPROGRESS: A non-blocking optimization call was made (by setting the NonBlocking parameter to 1 in a Gurobi Compute Server environment), but the associated optimization run is not yet complete.
/// \li GRB_USER_OBJ_LIMIT: Optimization terminated because the user objective function has reached the value specified in the ObjLimit parameter and no better solution has been found.
/// \li GRB_WORK_LIMIT: Optimization terminated because the total work expended exceeded the value specified in the WorkLimit parameter.
int solveGRBModel(GRBModel& model, FILE* const logfile) {
  int optimstatus = -1;
  try {
    model.optimize();

    optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_INF_OR_UNBD) {
      const int presolve_flag = model.get(GRB_IntParam_Presolve);
      if (presolve_flag) {
        model.set(GRB_IntParam_Presolve, 0);
        model.optimize();
        optimstatus = model.get(GRB_IntAttr_Status);
      }
    }
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  return optimstatus;
} /* solveGRBModel */

void saveSolution(
    std::vector<double>& solution,
    const GRBModel& model) {
  GRBVar* vars = model.getVars(); // An array of all variables in the model. Note that this array is heap-allocated, and must be returned to the heap by the user.
  const int num_vars = model.get(GRB_IntAttr_NumVars);
  // double* values = model.get(GRB_DoubleAttr_X, vars, num_vars); // error due to model being const
  // solution.assign(values, values + num_vars);
  solution.resize(num_vars);
  for (int i = 0; i < num_vars; i++) {
    solution[i] = vars[i].get(GRB_DoubleAttr_X);
  }
  if (vars) {
    delete[] vars;
  }
} /* saveSolution */

/// @details Sets model parameters for branch-and-bound tests
void setStrategyForBBTestGurobi(const Parameters& params, const int strategy,
    GRBModel& model, const double best_bound, int seed) {
  if (seed < 0) seed = params.get(intParam::RANDOM_SEED);
  // Parameters that should always be set
  model.set(GRB_DoubleParam_TimeLimit, params.get(doubleConst::BB_TIMELIMIT)); // time limit
  model.set(GRB_IntParam_Threads, 1); // single-threaded
  if (seed >= 0) {
    model.set(GRB_IntParam_Seed, seed); // random seed
  }
//  model.set(GRB_DoubleParam_MIPGap, param.getEPS()); // I guess the default 1e-4 is okay, though it messes up for large objective values

  if (params.get(VERBOSITY) == 0) {
    model.set(GRB_IntParam_OutputFlag, 0); // turn off output
  } else {
    model.set(GRB_IntParam_OutputFlag, 1); // turn on output
    model.set(GRB_IntParam_DisplayInterval, 1); // set display frequency
  }

  if (strategy <= 0) {
    // Default strategy
  } else {
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      //model.set(GRB_IntParam_DualReductions, 0); // disable dual reductions
      model.set(GRB_IntParam_PreCrush, 1); // must be enabled when using user cuts
    }

    // Turn off all cuts
    if (use_bb_option(strategy, BB_Strategy_Options::all_cuts_off)) {
      model.set(GRB_IntParam_Cuts, 0); // turn off all cuts
    }

    // Presolve
    if (use_bb_option(strategy, BB_Strategy_Options::presolve_off)) {
      model.set(GRB_IntParam_Presolve, 0); // disable presolve
    }

    // Heuristics
    if (use_bb_option(strategy, BB_Strategy_Options::heuristics_off)) {
      model.set(GRB_DoubleParam_Heuristics, 0.); // disable heuristics
      //model.set(GRB_IntParam_PumpPasses, 0); // disable feasibility pump; is this turned off by setting heuristic frequency to 0.?
      //model.set(GRB_IntParam_RINS, 0); // disable RINS; is this turned off by setting heuristic frequency to 0.?
    }

    // Feed the solver the best bound provided
    if (use_bb_option(strategy, BB_Strategy_Options::use_best_bound)) {
      if (!isInfinity(std::abs(best_bound))) {
        // BestObjStop: stop when primal bound <= z
        // BestBdStop: stop when dual bound >= z
        // Cutoff: prune subtrees with objective value > z
        //model.set(GRB_DoubleParam_BestObjStop, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
        model.set(GRB_DoubleParam_BestBdStop, best_bound - 1e-7); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
        //model.set(GRB_DoubleParam_Cutoff, best_bound + 1e-3); // give the solver the best IP objective value (it is a minimization problem) with a tolerance
      }
      // Check if user provides mip start or solution file
      std::string solfile = params.get(stringParam::SOLFILE);
      std::string ext1 = "_gurobi.sol.gz";
      std::string ext2 = "_gurobi.sol";
      std::string ext3 = "_gurobi.mst.gz";
      std::string ext4 = "_gurobi.mst";
      bool user_provides_start = false;
      user_provides_start |= (solfile.size() > ext1.size()) && (solfile.compare(solfile.size() - ext1.size(), ext1.size(), ext1) == 0);
      user_provides_start |= (solfile.size() > ext2.size()) && (solfile.compare(solfile.size() - ext2.size(), ext2.size(), ext2) == 0);
      user_provides_start |= (solfile.size() > ext3.size()) && (solfile.compare(solfile.size() - ext3.size(), ext3.size(), ext3) == 0);
      user_provides_start |= (solfile.size() > ext4.size()) && (solfile.compare(solfile.size() - ext4.size(), ext4.size(), ext4) == 0);
      if (user_provides_start) {
        printf("Gurobi: Reading optimal solution from %s.\n", solfile.c_str());
        model.read(solfile);
      }
    }
  } /* strategy > 0 */

  // Check if we should use strong branching
  if (use_bb_option(std::abs(strategy), BB_Strategy_Options::strong_branching_on)) {
    model.set(GRB_IntParam_VarBranch, 3); // turn on strong branching
  }
} /* setStrategyForBBTestGurobi */

/**
 * User cut callback
 */
class GurobiUserCutCallback : public GRBCallback {
  public:
    Parameters params;
    int num_vars;
    double obj_offset;
    GRBVar* grb_vars;
    const OsiCuts* cuts;
    bool addAsLazy;
    double first_lp_opt;
    BBInfo info;

    GurobiUserCutCallback(const Parameters& params, int num_vars, double obj_offset,
        GRBVar* vars, const OsiCuts* cutsToAdd, const bool lazyFlag = false)
        : params(params), num_vars(num_vars), obj_offset(obj_offset), grb_vars(vars),
        cuts(cutsToAdd), addAsLazy(lazyFlag) {
      this->info.obj = params.get(doubleConst::INF);
      this->first_lp_opt = -1 * params.get(doubleConst::INF);
      this->info.root_passes = 0;
      this->info.first_cut_pass = -1 * params.get(doubleConst::INF);
      this->info.last_cut_pass = -1 * params.get(doubleConst::INF);
      this->info.root_iters = 0;
      this->info.root_time = 0.;
      this->info.last_sol_time = 0.;
//      this->numTimesApplied = 0;
    } /* constructor */

    /**
     * Destructor
     */
    ~GurobiUserCutCallback() {
      if (this->grb_vars) {
        delete[] grb_vars;
      }
    } /* destructor */

  protected:
    double getObjValue() {
      double* vals = GRBCallback::getNodeRel(grb_vars, num_vars);
      double obj_val = obj_offset;
      for (int i = 0; i < num_vars; i++) {
        obj_val += vals[i] * grb_vars[i].get(GRB_DoubleAttr_Obj);
      }
      if (vals) {
        delete[] vals;
      }
      return obj_val;
    } /* getObjValue */

    void callback() {
      try {
        if (where == GRB_CB_MIPSOL) {
          const double obj = GRBCallback::getDoubleInfo(GRB_CB_MIPSOL_OBJ);
          if (obj < this->info.obj) { // integer-feasible solution with better solution has been found
            this->info.last_sol_time = GRBCallback::getDoubleInfo(GRB_CB_RUNTIME);
            this->info.obj = obj;
          }
        }
        else if (where == GRB_CB_MIP) {
          const int num_nodes = GRBCallback::getDoubleInfo(GRB_CB_MIP_NODCNT);
          if (num_nodes > 0) {
            return;
          }
          this->info.root_iters = GRBCallback::getDoubleInfo(GRB_CB_MIP_ITRCNT);
        }
        else if (where == GRB_CB_MIPNODE) {
          const int num_nodes = GRBCallback::getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
          if (num_nodes > 0) {
            return;
          }
          this->info.root_passes++;
          this->info.root_time = GRBCallback::getDoubleInfo(GRB_CB_RUNTIME);

          if (GRBCallback::getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) {
            return;
          }

          const double objValue = getObjValue();
          this->info.last_cut_pass = objValue;

          // Make sure this is our first entry into the root
          if (this->info.root_passes > 1) {
            if (this->info.root_passes == 2) {
              this->info.first_cut_pass = objValue;
            }
            return;
          }

          this->first_lp_opt = objValue;
          this->info.first_cut_pass  = objValue;

          // Add cuts to the model, one at a time
          if (cuts) {
            for (int cut_ind = 0; cut_ind < cuts->sizeCuts(); cut_ind++) {
              const OsiRowCut* curr_cut = cuts->rowCutPtr(cut_ind);
              const int num_el = curr_cut->row().getNumElements();
              const int* ind = curr_cut->row().getIndices();
              const double* vals = curr_cut->row().getElements();

              // Cannot add it all at once due to the C++ interface being worse than the C one
              //GRBaddconstr(model, num_el, ind, vals, GRB_GREATER_EQUAL, curr_cut->rhs(), NULL);
              GRBLinExpr lhs = 0;
              for (int i = 0; i < num_el; i++) {
                const int curr_var = ind[i];
                const double curr_value = vals[i];
                lhs += curr_value * grb_vars[curr_var];
              }
              // model.addConstr(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              addCut(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              if (addAsLazy) {
                addCut(lhs, GRB_GREATER_EQUAL, curr_cut->rhs());
              }
            } /* add cuts to the model */
          } /* check that cuts is not null */
        } /* where == MIPNODE (make sure we are in the right "where") */
//        else if (where == GRB_CB_MIP) {
//          // Get the number of times cuts were applied
//          num_total_cuts_applied = GRBCallback::getIntInfo(GRB_CB_MIP_CUTCNT);
//        } /* where == MIP (make sure we are in the right "where") */
//        else if (where == GRB_CB_MESSAGE) {
//          // Number user cuts
//          std::string msg = GRBCallback::getStringInfo(GRB_CB_MSG_STRING);
//        } /* where == MESSAGE */
      } catch (GRBException& e) {
        error_msg(errorstring, "Gurobi: Error during callback: %s\n", e.getMessage().c_str());
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } catch (...) {
        error_msg(errorstring, "Gurobi: Error during callback.\n");
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    } /* callback */
}; /* class GurobiUserCutCallback */

void presolveModelWithGurobi(const Parameters& params, int strategy, 
    GRBModel& model, double& presolved_lp_opt, std::string& presolved_name,
    const double best_bound) {
//#ifdef TRACE
  printf("\n## Gurobi: Presolving model ##\n");
//#endif
  try {
    GRBModel presolved_model = model.presolve();
    GRBModel presolved_model_mip = presolved_model;
    GRBVar* vars = presolved_model.getVars();
    const int num_vars = presolved_model.get(GRB_IntAttr_NumVars);
    for (int i = 0; i < num_vars; i++) {
      vars[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }

    strategy = static_cast<int>(BB_Strategy_Options::presolve_on);
    setStrategyForBBTestGurobi(params, strategy, presolved_model, best_bound);

    printf("Solving Gurobi-presolved model to get new LP optimum value.\n");
    presolved_model.optimize();

    // Save optimal value
    int optimstatus = presolved_model.get(GRB_IntAttr_Status);
    if (optimstatus == GRB_OPTIMAL) {
        presolved_lp_opt = presolved_model.get(GRB_DoubleAttr_ObjVal);
        size_t slashindex = presolved_name.find_last_of("/\\");
        presolved_model_mip.set(GRB_StringAttr_ModelName, presolved_name.substr(slashindex+1));
        printf("Saving Gurobi-presolved model to \"%s.mps.gz\".\n", presolved_name.c_str());
        createTmpFileCopy(params.logfile, presolved_model_mip, presolved_name, ".mps.gz"); // adds .mps.gz ext
        if (vars) {
          delete[] vars;
        }
    } else {
      error_msg(errorstring, "Gurobi: Was not able to solve presolved model to optimality. Status: %d.\n", optimstatus);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* presolveModelWithGurobi (GRBModel) */

void presolveModelWithGurobi(const Parameters& params, int strategy, const char* f_name,
    double& presolved_lp_opt, std::string& presolved_name, const double best_bound) {
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    presolveModelWithGurobi(params, strategy, model, presolved_lp_opt, presolved_name, best_bound);
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* presolveModelWithGurobi (filename) */

void presolveModelWithGurobi(const Parameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_lp_opt,
    std::string& presolved_name, const double best_bound) {
  std::string f_name = "";
  createTmpFileCopy(params.logfile, solver, f_name);
  presolveModelWithGurobi(params, strategy, f_name.c_str(), presolved_lp_opt, presolved_name, best_bound);
  remove(f_name.c_str()); // remove temporary file
} /* presolveModelWithGurobi (Osi) */

void doBranchAndBoundWithGurobi(const Parameters& params, int strategy,
    GRBModel& model, BBInfo& info, const double best_bound, 
    std::vector<double>* const solution = NULL) {
//#ifdef TRACE
  printf("\tRunning B&B with Gurobi. Strategy: %d. Random seed: %d.\n",
      strategy, params.get(intParam::RANDOM_SEED));
//#endif
  try {
    setStrategyForBBTestGurobi(params, strategy, model, best_bound);

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_INF_OR_UNBD) {
      const int presolve_flag = model.get(GRB_IntParam_Presolve);
      if (presolve_flag) {
        model.set(GRB_IntParam_Presolve, 0);
        model.optimize();
        optimstatus = model.get(GRB_IntAttr_Status);
      }
    }

    if (optimstatus == GRB_INFEASIBLE || optimstatus == GRB_UNBOUNDED || optimstatus == GRB_INF_OR_UNBD) {
      error_msg(errorstring, "Gurobi: Failed to optimize MIP.\n");
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }

    switch (optimstatus) {
      case GRB_CUTOFF:
      case GRB_ITERATION_LIMIT:
      case GRB_NODE_LIMIT:
      case GRB_TIME_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_USER_OBJ_LIMIT: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      case GRB_OPTIMAL: {
        info.obj = model.get(GRB_DoubleAttr_ObjVal);
        info.bound = model.get(GRB_DoubleAttr_ObjBound);
        break;
      }
      default: {
        error_msg(errorstring, "Gurobi: Other status after solve: %d.\n", optimstatus);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      }
    } // switch optimistatus
    info.iters = (long) model.get(GRB_DoubleAttr_IterCount);
    info.nodes = (long) model.get(GRB_DoubleAttr_NodeCount);
    info.time = model.get(GRB_DoubleAttr_Runtime);

#ifdef TRACE
    printf("Gurobi: Solution value: %s.\n", stringValue(info.obj, "%1.6f").c_str());
    printf("Gurobi: Best bound: %s.\n", stringValue(info.bound, "%1.6f").c_str());
    printf("Gurobi: Number iterations: %ld.\n", info.iters);
    printf("Gurobi: Number nodes: %ld.\n", info.nodes);
    printf("Gurobi: Time: %f.\n", info.time);
#endif

   // Save the solution if needed
   if (solution) {
     // Get variables and double-check objective matches
     double obj = model.get(GRB_DoubleAttr_ObjCon);
     GRBVar* vars = model.getVars();
     const int num_vars = model.get(GRB_IntAttr_NumVars);
     (*solution).resize(num_vars);
     for (int i = 0; i < num_vars; i++) {
       (*solution)[i] = vars[i].get(GRB_DoubleAttr_X);
       obj += (*solution)[i] * vars[i].get(GRB_DoubleAttr_Obj);
     }
     if (vars) {
       delete[] vars;
     }
     const bool add_same = isVal(obj, info.obj, 1e-3);
     const bool mult_same = ( (isZero(obj) || isZero(info.obj)) && add_same )
       || ( (obj/info.obj > 0) && isVal(std::abs(obj/info.obj), 1., .01) ); 
     if (!add_same && !mult_same) {
       error_msg(errorstring, "Gurobi: %g, computed objective value from solution, does not match solver's obj value %g.\n", obj, info.obj);
       writeErrorToLog(errorstring, params.logfile);
       exit(1);
     }
     if (!isVal(obj, info.obj)) {
       warning_msg(warnstring, "Gurobi: %g, computed objective value from solution, does not match solver's obj value %g.\n", obj, info.obj);
     }
   } // save ip solution

   if (use_temp_option(params.get(intParam::TEMP), TempOptions::SAVE_IP_OPT) && params.get(stringParam::SOLFILE).empty()) {
     std::string dir, instname, in_file_ext;
     std::string filename = params.get(stringParam::FILENAME);
     parseFilename(dir, instname, in_file_ext, filename, params.logfile);
     std::string f_name = dir + "/" + instname + "_gurobi.mst.gz";
     if (!fexists(f_name.c_str())) {
       printf("Saving Gurobi MIP start file to \"%s\".\n", f_name.c_str());
       model.write(f_name);
     }
   }
  } catch(GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (GRBModel) */

void doBranchAndBoundWithUserCutsGurobi(const Parameters& params,
    int strategy, GRBModel& model, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy,
    std::vector<double>* const solution = NULL) {
  // Ensure that user cuts setting is enabled
  if (!use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
    warning_msg(warnstring,
        "Need to use user_cuts option; strategy currently: %d.\n", strategy);
    strategy = enable_bb_option(strategy, BB_Strategy_Options::user_cuts);
  }

  try {
    // Due to using the C++ interface, we need to access the variable list first
    GRBVar* grb_vars = model.getVars();
    GurobiUserCutCallback cb = GurobiUserCutCallback(params,
        model.get(GRB_IntAttr_NumVars), model.get(GRB_DoubleAttr_ObjCon),
        grb_vars, cuts, addAsLazy);
    model.setCallback(&cb);
    if (grb_vars) { delete[] grb_vars; }

    // Update the model
    model.update();
    /*
    const int retcode = GRBupdatemodel(&model);
    if (retcode) {
      error_msg(errorstring, "Gurobi: Error updating the model; error code: %d.\n", retcode);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
    */

    // Finally, run B&B
    doBranchAndBoundWithGurobi(params, strategy, model, info, best_bound, solution);

    // Save information
    info.root_passes = cb.info.root_passes;
    if (info.root_passes > 0) {
      info.first_cut_pass = cb.info.first_cut_pass; // second because first is lp opt val
      info.last_cut_pass = cb.info.last_cut_pass;
      info.root_iters = cb.info.root_iters;
      info.root_time = cb.info.root_time;
      info.last_sol_time = cb.info.last_sol_time;
    } else {
      // I guess this can happen if we solve during presolve,
      // or if we do no root passes of cuts
      if (info.nodes == 0) {
        info.first_cut_pass = info.obj;
        info.last_cut_pass = info.obj;
        info.root_iters = info.iters; // all iters were at the root
        info.root_time = info.time; // all time was spent at the root
        info.last_sol_time = cb.info.last_sol_time; // roughly the same as total time in this case
      } else {
        info.first_cut_pass = cb.first_lp_opt;
        info.last_cut_pass = cb.first_lp_opt;
        info.root_iters = cb.info.root_iters;
        info.root_time = cb.info.root_time; // all time was spent at the root
        info.last_sol_time = cb.info.last_sol_time; // roughly the same as total time in this case
      }
    }
//#ifdef TRACE
//    printf("Gurobi: Times cuts applied: %d.\n", cb.numTimesApplied);
//#endif
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithUserCutsGurobi (GRBModel) */

void doBranchAndBoundWithGurobi(const Parameters& params, int strategy,
    const char* f_name, BBInfo& info, const double best_bound,
    std::vector<double>* const solution) {
#ifdef TRACE
  printf("\n## Reading from file into Gurobi. ##\n");
#endif
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    if (use_bb_option(strategy, BB_Strategy_Options::user_cuts)) {
      doBranchAndBoundWithUserCutsGurobi(params, strategy, model, NULL, info, best_bound, false, solution);
    } else {
      doBranchAndBoundWithGurobi(params, strategy, model, info, best_bound, solution);
    }
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithGurobi (filename) */

void doBranchAndBoundWithGurobi(const Parameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound, std::vector<double>* const solution) {
  std::string f_name;
  createTmpFileCopy(params.logfile, solver, f_name);
  doBranchAndBoundWithGurobi(params, strategy, f_name.c_str(), info, best_bound,
      solution);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithGurobi (Osi) */

void doBranchAndBoundWithUserCutsGurobi(const Parameters& params,
    int strategy, const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
#ifdef TRACE
  printf("\n## Reading from file into Gurobi and adding user cuts. ##\n");
#endif
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env, f_name);
    doBranchAndBoundWithUserCutsGurobi(params, strategy, model, cuts, info, best_bound, addAsLazy);
  } catch (GRBException& e) {
    error_msg(errorstring, "Gurobi: Exception caught: %s\n", e.getMessage().c_str());
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  } catch (...) {
    error_msg(errorstring, "Gurobi: Unknown exception caught.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
} /* doBranchAndBoundWithUserCutsGurobi (filename) */

void doBranchAndBoundWithUserCutsGurobi(const Parameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy) {
  std::string f_name;
  createTmpFileCopy(params.logfile, solver, f_name);
  doBranchAndBoundWithUserCutsGurobi(params, strategy, f_name.c_str(), cuts, info, best_bound, addAsLazy);
  remove(f_name.c_str()); // remove temporary file
} /* doBranchAndBoundWithUserCutsGurobi (Osi) */
#endif /* USE_GUROBI */
