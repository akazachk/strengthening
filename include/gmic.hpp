/**
 * @file gmic.hpp
 * @author A. M. Kazachkov
 * @date 2019-11-19
 */
#pragma once

#include <vector>
#include <cstdio>

class OsiCuts;
class OsiRowCut;
class OsiSolverInterface;

/// @brief enum for the different types of Gomory cuts
enum class GomoryType {
  CglGMI = 1, ///< Generate GMICs via CglGMI
  CreateMIG_CustomStrengthen = 2, ///< Generate GMICs via createMIG and apply custom strengthening via strengthen.hpp
  CreateMIG_ClosedFormStrengthen = 3, ///< Generate GMICs via createMIG and apply closed-form strengthening
};

/// @brief Generate Gomory cuts based on \p option (0: CglGMI, 1: ) 
void generateGomoryCuts(
    OsiCuts& currGMICs,
    const OsiSolverInterface* const si,
    const int option,
    const int strengthen_option,
    const int MIN_SUPPORT_THRESHOLD,
    const double MAX_SUPPORT_REL,
    const double AWAY,
    const double DIFFEPS,
    FILE* logfile);

/// @brief Adapted from CglLandP: create Gomory mixed-integer cut
void createMIG(OsiRowCut &cut, const OsiSolverInterface* const solver,
    const int splitVarIndex, const int splitVarRowIndex, const bool strengthen = true,
    FILE* logfile = NULL);

/// @brief From CglLandPSimplex: put a vector into structural space
void eliminate_slacks(std::vector<double>& vec, const OsiSolverInterface* const solver);

/// @brief From CglLandP: return the coefficients of the intersection cut
double unstrengthenedIntersectionCutCoeff(double abar, double f0);

/// @brief From CglLandP: compute the modularized row coefficient for an integer variable
double modularizedCoeff(double abar, double f0);

/// @brief Adapted from CglLandP: return the coefficients of the unstrengthened intersection cut
double unstrengthenedIntersectionCutCoeff(double abar, double f0);

/// @brief Adapted from CglLandP: return the coefficients of the strengthened intersection cut
double strengthenedIntersectionCutCoeff(double abar, double f0);

/// @brief Return #strengthenedIntersectionCutCoeff if \p strengthen is true; otherwise #unstrengthenedIntersectionCutCoeff
double intersectionCutCoeff(double abar, double f0, const bool strengthen = false);

/// @brief Return #strengthenedIntersectionCutCoeff if \p strengthen is true and \p currFracVar is structural integer; otherwise #unstrengthenedIntersectionCutCoeff
double intersectionCutCoeff(double abar, double f0,
    const OsiSolverInterface* const solver = NULL, 
    const int currFracVar = -1,
    const bool strengthen = false);

