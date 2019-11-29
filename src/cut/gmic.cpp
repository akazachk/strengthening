/**
 * @file gmic.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-19
 */
#include "gmic.hpp"

#include <OsiRowCut.hpp>

#include "Parameters.hpp" // SolverInterface
#include "SolverHelper.hpp" // isBasicVar
#include "utility.hpp"

/** Adapted from CglLandP **/
void createMIG(
    OsiRowCut &cut, 
    const OsiSolverInterface* const solver,
    const int splitVarIndex,
    const int splitVarRowIndex,
    const bool strengthen) {
  const int numcols = solver->getNumCols();
  const int numrows = solver->getNumRows();

  double splitVarVal = solver->getColSolution()[splitVarIndex];

  std::vector<double> coeff(solver->getNumCols());

  const double* colLower = solver->getColLower();
  const double* rowLower = solver->getRowLower();
  const double* colUpper = solver->getColUpper();
  const double* rowUpper = solver->getRowUpper();

  bool mustDelete = false;
  const CoinWarmStartBasis* basis_;
  try {
    const SolverInterface* si = dynamic_cast<const SolverInterface*>(solver);
    basis_ = dynamic_cast<const CoinWarmStartBasis*>(si->getConstPointerToWarmStart());
  } catch (std::exception& e) {
    error_msg(errorstring, "Failed to cast to SolverInterface.\n");
    exit(1);
  }

  std::vector<double> basisRowStruct(numcols), basisRowSlack(numrows);
  solver->enableFactorization();
  solver->getBInvARow(splitVarRowIndex, &basisRowStruct[0], &basisRowSlack[0]);
  solver->disableFactorization();

  double f0 = splitVarVal - std::floor(splitVarVal);

  cut.setUb(COIN_DBL_MAX);
  std::vector<double> vec(numcols + numrows, 0.);
  double cutRhs = 1.;
  assert(std::abs(cutRhs) < 1e100); // TODO assert gracefully
  for (int var = 0; var < solver->getNumCols() + solver->getNumRows(); var++) {
    if (isBasicVar(solver, var)) { continue; }
    double value = 0.0;
    if (var < solver->getNumCols()) {
      const CoinWarmStartBasis::Status status = basis_->getStructStatus(var);
      if (status == CoinWarmStartBasis::atUpperBound) { // TODO I think this might be wrong for ub vars
        value = -intersectionCutCoeff(-basisRowStruct[var], f0, solver, var, strengthen);
        cutRhs += value * colUpper[var];
      } else if (status == CoinWarmStartBasis::atLowerBound) {
        value = intersectionCutCoeff(basisRowStruct[var], f0, solver, var, strengthen);
        cutRhs += value * colLower[var];
      } else {
        error_msg(errorstring, "createMIG: Invalid basis\n");
        if (mustDelete && basis_) {
          delete basis_;
        }
        exit(1); // probably better to throw than exit
      }
    } else {
      const int row = var - solver->getNumCols();
      if (lessThanVal(std::abs(basisRowSlack[row]), 0.0)) {
        continue;
      }
      if (!isInfinity(rowUpper[row])) {
        value = intersectionCutCoeff(basisRowSlack[row], f0, solver, var, strengthen);
        cutRhs -= value * rowUpper[row];
      } else {
        value = -intersectionCutCoeff(-basisRowSlack[row], f0, solver, var, strengthen);
        cutRhs -= value * rowLower[row];
        assert(
            basis_->getArtifStatus(row) == CoinWarmStartBasis::atUpperBound
            || !isInfinity(rowUpper[row])); // TODO assert gracefully
      }
    }
    assert(std::abs(cutRhs) < 1e100); // TODO assert gracefully
    vec[var] = value; // had "original index" but we do not work in subspace
  }

  //Eliminate slacks
  eliminate_slacks(vec, solver);

  //Pack vec into the cut
  std::vector<int> inds(numcols);
  int nelem = 0;
  for (int i = 0; i < numcols; i++) {
    if (std::abs(vec[i]) > COIN_INDEXED_TINY_ELEMENT) {
      vec[nelem] = vec[i];
      inds[nelem++] = i;
    }
  }

  cut.setLb(cutRhs);
  cut.setRow(nelem, inds.data(), vec.data(), false);
  if (mustDelete && basis_) {
    delete basis_;
  }
} /* createMIG */

void
eliminate_slacks(std::vector<double>& vec, const OsiSolverInterface* const solver)
{
    const CoinPackedMatrix * mat = solver->getMatrixByCol();
    const CoinBigIndex * starts = mat->getVectorStarts();
    const int * lengths = mat->getVectorLengths();
    const double * values = mat->getElements();
    const CoinBigIndex * indices = mat->getIndices();
    const double * vecSlacks = vec.data() + solver->getNumCols();
    for (int j = 0 ; j < solver->getNumCols() ; j++)
    {
        const CoinBigIndex& start = starts[j];
        CoinBigIndex end = start + lengths[j];
        double & val = vec[j]; //vec[original_index_[j]];
        for (CoinBigIndex k = start ; k < end ; k++)
        {
            val -= vecSlacks[indices[k]] * values[k];
        }
    }
} /* eliminate_slacks */

/** From CglLandP: return the coefficients of the intersection cut */
double unstrengthenedIntersectionCutCoeff(double abar, double f0) {
  if (abar > 0) {
    //return alpha_i * (1 - beta);
    return (abar / f0);
  } else {
    //return -alpha_i * beta;
    return (-abar / (1 - f0));
  }
} /* unstrengthenedIntersectionCutCoeff */

/** From CglLandP compute the modularized row coefficient for an integer variable */
double modularizedCoeff(double abar, double f0) {
  double f_i = abar - floor(abar);
  if (f_i <= f0) {
    return f_i;
  } else {
    return f_i - 1;
  }
} /* modularizedCoeff */

/** Adapted from CglLandP: return the coefficients of the strengthened intersection cut */
double strengthenedIntersectionCutCoeff(double abar, double f0) {
    //const OsiSolverInterface* const solver, 
    //const int currFracVar) {
  //if ((currFracVar >= 0)
  //    && ((currFracVar >= solver->getNumCols())
  //        || !(solver->isInteger(currFracVar)))) {
  //  return unstrengthenedIntersectionCutCoeff(abar, f0);
  //} else {
//    return modularizedCoeff(abar, f0);
    double f_i = abar - std::floor(abar);
    if (f_i < f0) {
//      return f_i * (1 - f0);
      return f_i / f0;
    } else {
//      return (1 - f_i) * f0;
      return (1 - f_i) / (1 - f0);
    }
  //}
} /* strengthenedIntersectionCutCoeff */

double intersectionCutCoeff(double abar, double f0, const bool strengthen) {
  if (!strengthen) {
    return unstrengthenedIntersectionCutCoeff(abar, f0);
  } else {
    return strengthenedIntersectionCutCoeff(abar, f0);
  }
} /* intersectionCutCoeff */

double intersectionCutCoeff(double abar, double f0,
    const OsiSolverInterface* const solver, const int currFracVar,
    const bool strengthen) {
  if (!strengthen || (currFracVar >= solver->getNumCols())
      || !(solver->isInteger(currFracVar))) {
    return unstrengthenedIntersectionCutCoeff(abar, f0);
  } else {
    return strengthenedIntersectionCutCoeff(abar, f0);
  }
} /* intersectionCutCoeff */
