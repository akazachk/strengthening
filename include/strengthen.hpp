/**
 * @file strengthen.hpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#pragma once

#include <vector>
#include <cstdio>

// COIN-OR files
class OsiSolverInterface;
class CoinPackedMatrix;

// Project files
#include "CutCertificate.hpp" // TermCutCertificate and CutCertificate
class Disjunction;
class DisjunctiveTerm;

/// @brief Return Farkas certificate for a single term, given a valid cut (in sparse form) 
/// and a linear program describing the feasible region of a polyhedron
void getCertificate(TermCutCertificate& v, const int num_elem, const int* const ind, 
    const double* const coeff, OsiSolverInterface* const solver, FILE* logfile);

/// @brief Not tested
void getCertificateForTerm(
    TermCutCertificate& v, 
    const int num_elem, 
    const int* const ind, 
    const double* const coeff,
    const OsiSolverInterface* const si,
    const DisjunctiveTerm* const term,
    const double DIFFEPS,
    FILE* logfile);

/// @brief Not implemented
void getCertificateTrivial(TermCutCertificate& v, const int num_elem, const int* const ind, 
    const double* const coeff, const OsiSolverInterface* const solver, const Disjunction* const disj);

/// @brief Attempt to strengthen coefficients of given cut
int strengthenCut(
    std::vector<double>& str_coeff,
    double& str_rhs,
    const int num_elem, 
    const int* const ind, 
    const double* const coeff,
    const double rhs,
    const Disjunction* const disj,
    const CutCertificate& v, 
    const OsiSolverInterface* const solver,
    FILE* logfile,
    const std::vector<double>& ip_solution = {});

/// @brief Calculate new cut coefficient
bool strengthenCutCoefficient(
    double& str_coeff,
    double& str_rhs,
    const int var,
    const double coeff,
    const Disjunction* const disj,
    const std::vector<double>& lb_term,
    const CutCertificate& v, 
    const OsiSolverInterface* const solver,
    OsiSolverInterface* const mono,
    FILE* logfile);

/// @brief Creates monoidal strengthening IP
void setupMonoidalIP(
    OsiSolverInterface* const mono,
    const int var,
    const Disjunction* const disj,
    const std::vector<double>& lb_term,
    const CutCertificate& v, 
    const OsiSolverInterface* const solver,
    const double mult);
    
/// @brief Change relevant rhs for rows of monoidal IP set up as in setupMonoidalIP
void updateMonoidalIP(
    OsiSolverInterface* const mono,
    const int var,
    const Disjunction* const disj,
    const CutCertificate& v, 
    const OsiSolverInterface* const solver,
    const double mult);
