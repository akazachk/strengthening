/**
 * @file disjcuts.cpp
 * @author A. M. Kazachkov
 * @date 2019-12-27
 */
#include "disjcuts.hpp"

// COIN-OR includes
#include <OsiSolverInterface.hpp>
#include <OsiCuts.hpp>

// VPC includes
//#include "OsiProblemData.hpp"
#include "Disjunction.hpp"
#include "VPCSolverInterface.hpp"
#include "VPCParameters.hpp"

/**
 * @brief Generate disjunctive cuts
 *
 * @return Generated cuts and disjunction used to generate them
 */
void genDisjCuts(
    /// [in]
    const OsiSolverInterface* const si,
    /// [in]
    const int numdisjterms,
    /// [in]
    const int maxcuts,
    /// [in/out]
    OsiCuts* const cuts,
    /// [out] disjunction used to generate cuts (**user must free**)
    Disjunction* disj) {
   std::unique_ptr<VPCSolverInterface> solver = std::make_unique<VPCSolverInterface>();
   solver->load(si);
   solver->solver->initialSolve();
   const double init_obj = solver->getObjValue();

   VPCParametersNamespace::VPCParameters* params = solver->params;
   params->set(VPCParametersNamespace::DISJ_TERMS, numdisjterms);
   params->set(VPCParametersNamespace::CUTLIMIT, maxcuts);

   solver->generateCuts();
   solver->applyCuts();
   solver->solver->resolve();
   const double final_obj = solver->getObjValue();
   fprintf(stdout, "Initial obj value: %g. Final obj value: %g.\n", init_obj, final_obj);

   cuts->insert(*(solver->cuts));
   if (solver->disj)
     disj = solver->disj->clone();
} /* genDisjCuts */

