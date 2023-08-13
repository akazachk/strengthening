/**
 * @file GurobiHelper.hpp
 * @author A. M. Kazachkov
 * @date 2019-Feb-28
 */
#pragma once

#include <limits>
#include <string>
#include <vector>

class OsiSolverInterface;
class OsiCuts;

struct BBInfo;
namespace StrengtheningParameters {
  struct Parameters;
}

#ifdef USE_GUROBI
class GRBModel;

/// @brief Use OsiSolverInterface to build GRBModel instance
GRBModel* buildGRBModelFromOsi(const OsiSolverInterface* const solver, FILE* const logfile);

/// @brief Solve GRBModel \p model and return status
int solveGRBModel(GRBModel& model, FILE* const logfile);

/// @brief Save model solution into double vector
void saveSolution(std::vector<double>& solution, const GRBModel& model);

/// @brief Set \p model parameters
void setStrategyForBBTestGurobi(
    const StrengtheningParameters::Parameters& params,
    const int strategy,
    GRBModel& model,
    const double best_bound = std::numeric_limits<double>::max(),
    int seed = -1);

void presolveModelWithGurobi(const StrengtheningParameters::Parameters& params, int strategy,
    const char* f_name, double& presolved_lp_opt, std::string& presolved_name,
    const double best_bound);
void presolveModelWithGurobi(const StrengtheningParameters::Parameters& params, int strategy,
    const OsiSolverInterface* const solver, double& presolved_lp_opt, std::string& presolved_name);
void doBranchAndBoundWithGurobi(const StrengtheningParameters::Parameters& params, int strategy, 
    const char* f_name, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithGurobi(const StrengtheningParameters::Parameters& params, int strategy,
    const OsiSolverInterface* const solver, BBInfo& info,
    const double best_bound = std::numeric_limits<double>::max(),
    std::vector<double>* const solution = NULL);
void doBranchAndBoundWithUserCutsGurobi(const StrengtheningParameters::Parameters& params, int strategy,
    const char* f_name, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
void doBranchAndBoundWithUserCutsGurobi(const StrengtheningParameters::Parameters& params, int strategy,
    const OsiSolverInterface* const solver, const OsiCuts* cuts, BBInfo& info,
    const double best_bound, const bool addAsLazy = false);
#endif /* USE_GUROBI */
