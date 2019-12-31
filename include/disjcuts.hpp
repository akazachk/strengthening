/**
 * @file disjcuts.hpp
 * @author A. M. Kazachkov
 * @date 2019-12-27
 */
#pragma once

class OsiSolverInterface;
class OsiCuts;
class Disjunction;

void genDisjCuts(const OsiSolverInterface* const si,
    const int numdisjterms,
    const int maxcuts,
    OsiCuts* const cuts,
    Disjunction* disj);
