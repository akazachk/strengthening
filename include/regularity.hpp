/**
 * @file cglp.hpp
 * @author A. M. Kazachkov
 * @date 2023-08-02
 */
#pragma once

#include <cstdio> // FILE

class OsiSolverInterface;
class OsiRowCut;
class Disjunction;

/**
 * @brief Generate a CGLP from a cut
 */
void genRCVMILPFromCut(OsiSolverInterface* liftingSolver,
    const OsiRowCut* const cut, const Disjunction* const disj, 
    const OsiSolverInterface* const solver, FILE* const logfile,
    const bool use_min_sum_delta = true);
