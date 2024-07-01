/**
 * @file gmic.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-19
 */
#include "gmic.hpp"

#include <CglGMI.hpp>
#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>

#include "SplitDisjunction.hpp"
#include "SolverHelper.hpp" // isBasicVar, isRayVar, isNonBasicLBVar, isNonBasicUBVar
#include "SolverInterface.hpp"
#include "strengthen.hpp"
#include "utility.hpp"
#include "verify.hpp" // getCutFromCertificate

/** 
 * @details Generate Gomory cuts based on selected \p option (see #GomoryType) and \p strengthen_option
 *
 * Option 1: GglGMI
 * Option 2: custom generate intersection cuts, calculate Farkas certificate, do strengthening
 * Option 3: custom generate intersection cuts, calculate Farkas certificate, do closed-form strengthening
 */
void generateGomoryCuts(
    // [out] Gomory cuts that we generate
    OsiCuts& currGMICs,
    /// [in/out] solver from which we generate cuts; assumed optimal; non-const because we need to enable factorization
    OsiSolverInterface* const solver,
    /// [in] which GMIC generation method to use (see #GomoryType)
    const int option,
    /// [in] how to strengthen the cut
    const int strengthen_option,
    /// [in] CglGMI's sparsity parameter rejects with bad support when # cut coefficient nonzeros >= MAX_SUPPORT + MAX_SUPPORT_REL * # cols; for us MIN_SUPPORT == MIN_SUPPORT_THRESHOLD
    const int MIN_SUPPORT_THRESHOLD,
    /// [in] CglGMI's sparsity parameter rejects with bad support when # cut coefficient nonzeros >= MAX_SUPPORT + MAX_SUPPORT_REL * # cols
    const double MAX_SUPPORT_REL,
    /// [in] epsilon for deciding whether a value is fractional
    const double AWAY,
    /// [in] tolerance for reconstructing disjunctive term solution
    const double DIFFEPS,
    FILE* logfile) {
  fprintf(stdout, "\n## Starting GMIC generation with mode %d and strengthen option %d. ##\n", option, strengthen_option);

  StatVector num_coeff_str_stats;
  int num_cuts_strengthened = 0;
  initializeStats(num_coeff_str_stats);

  // Use CglGMI
  if (std::abs(option) == static_cast<int>(GomoryType::CglGMI)) {
    CglGMI GMIGen;
    // Set parameters so that many GMIs are generated
    // CglGMI's MAX_SUPPORT parameter is equivalent to the MIN_SUPPORT_THRESHOLD value in this code
    // CglGMI's sparsity parameter rejects with bad support when # cut coefficient nonzeros >= MAX_SUPPORT + MAX_SUPPORT_REL * # cols
    // This means that there will be a constant offset on the max allowable support for GMICs compared to VPCs (in the favor of GMICs being denser)
    GMIGen.getParam().setMAX_SUPPORT(MIN_SUPPORT_THRESHOLD);
    GMIGen.getParam().setMAX_SUPPORT_REL(MAX_SUPPORT_REL);
    GMIGen.getParam().setMAXDYN(solver->getInfinity());
    GMIGen.getParam().setMINVIOL(0.0);
    GMIGen.getParam().setAWAY(AWAY);
    GMIGen.generateCuts(*solver, currGMICs);
  } // generate GMICs via CglGMI

  // Generate GMICs via createMIG and apply custom strengthening via strengthen.hpp
  else if (std::abs(option) == static_cast<int>(GomoryType::CreateMIG_CustomStrengthen)) {
    solver->enableFactorization();

    // Get variables basic in which row
    std::vector<int> varBasicInRow(solver->getNumRows(), -1);
    solver->getBasics(&varBasicInRow[0]);
    solver->disableFactorization();

    // Loop over variables and generate cuts for any that are integer-restricted but take a fractional value
    for (int var = 0; var < solver->getNumCols(); var++) {
      if (!solver->isInteger(var)) {
        continue;
      }

      const double val = solver->getColSolution()[var];
      if (isVal(val, std::floor(val), AWAY) || isVal(val, std::ceil(val), AWAY)) {
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
        writeErrorToLog(errorstring, logfile);
        exit(1);
      }

      // Now we generate the unstrengthened intersection cut and see if we can strengthen it
      // The actual Farkas certificate for this cut can be derived in closed form from the basis
      OsiRowCut intCut;
      createMIG(intCut, solver, var, splitVarRowIndex, false, logfile);

      if (strengthen_option > 0) {
        // For each term of the disjunction, 
        // we need to explicitly add the constraint(s) defining the disjunctive term
        const CoinPackedVector lhs = intCut.row();

        std::vector<std::vector<double> > v(2); // [term][val] for each term, this will be of dimension rows + cols
        { // Check first side of the split
          // Calculate the certificate
          OsiSolverInterface* solver0 = solver->clone();
          const double el = -1.;
          solver0->addRow(1, &var, &el, -1. * std::floor(val), solver->getInfinity());

          getCertificate(v[0], lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), solver0, logfile);

          if (solver0) delete solver0;
        } // check first side of the split

        { // Check second side of the split
          // Calculate the certificate
          OsiSolverInterface* solver1 = solver->clone();
          const double el = 1.;
          solver1->addRow(1, &var, &el, std::ceil(val), solver->getInfinity());

          getCertificate(v[1], lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), solver1, logfile);

          if (solver1) delete solver1;
        } // check second side of the split

        // Do strengthening; for this we first need to setup a split disjunction
        SplitDisjunction disj;
        disj.var = var;
        disj.prepareDisjunction(solver);
        const double rhs = intCut.rhs();
        std::vector<double> str_coeff;
        double str_rhs;
        const int curr_num_coeffs_str = strengthenCut(str_coeff, str_rhs, lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), rhs, &disj, v, solver, logfile);

        // Update stats
        num_cuts_strengthened += (curr_num_coeffs_str > 0);
        updateStatsBeforeFinalize(num_coeff_str_stats, curr_num_coeffs_str);

        // Replace row if any coefficients were strengthened
        if (curr_num_coeffs_str > 0) {
          CoinPackedVector strCutCoeff(str_coeff.size(), str_coeff.data());
          intCut.setRow(strCutCoeff);
          intCut.setLb(str_rhs);
        }

        // We generate the strengthened intersection cut (to compare against)
        //OsiRowCut strIntCut;
        //createMIG(strIntCut, solver, var, splitVarRowIndex, true, logfile);
      } // do strengthening

      // Insert new cut into currGMICs
      currGMICs.insert(intCut);
    } // iterate over cols, generating GMICs
    finalizeStats(num_coeff_str_stats, currGMICs.sizeCuts());
  } // Generate GMICs via createMIG and apply custom strengthening via strengthen.hpp

  // Generate GMICs via createMIG and apply closed-form strengthening
  else if (std::abs(option) == static_cast<int>(GomoryType::CreateMIG_ClosedFormStrengthen)) {
    solver->enableFactorization();

    // Get variables basic in which row
    std::vector<int> varBasicInRow(solver->getNumRows(), -1);
    solver->getBasics(&varBasicInRow[0]);

    // Get nonbasic variables
    std::vector<int> NBVarIndex;
    NBVarIndex.reserve(solver->getNumCols());
    for (int var = 0; var < solver->getNumCols() + solver->getNumRows(); var++) {
      if (isRayVar(solver, var)) {
        NBVarIndex.push_back(var);
      }
    } // loop over vars, collecting which are nb

    // Collect basis inverse
    std::vector<std::vector<double> > currRay(NBVarIndex.size());
    for (int ray_ind = 0; ray_ind < (int) NBVarIndex.size(); ray_ind++) {
      const int NBVar = NBVarIndex[ray_ind];
      currRay[ray_ind].resize(solver->getNumRows());
      solver->getBInvACol(NBVar, &(currRay[ray_ind][0]));
      if (isNonBasicLBVar(solver, NBVar)) {
        for (int row = 0; row < solver->getNumRows(); row++) {
          currRay[ray_ind][row] *= -1;
        }
      }
    }
    solver->disableFactorization();

    for (int var = 0; var < solver->getNumCols(); var++) {
      if (!solver->isInteger(var)) {
        continue;
      }

      const double val = solver->getColSolution()[var];
      if (isVal(val, std::floor(val), AWAY) || isVal(val, std::ceil(val), AWAY)) {
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
        writeErrorToLog(errorstring, logfile);
        exit(1);
      }

      // Now we generate the unstrengthened intersection cut and see if we can strengthen it
      // The actual Farkas certificate for this cut can be derived in closed form from the basis
      OsiRowCut intCut;
      createMIG(intCut, solver, var, splitVarRowIndex, false, logfile);

      if (strengthen_option > 0) {
        std::vector<std::vector<double> > v(2); // [term][val] for each term, this will be of dimension rows + cols
        v[0].resize(solver->getNumRows() + 1 + solver->getNumCols(), 0.0);
        v[1].resize(solver->getNumRows() + 1 + solver->getNumCols(), 0.0);

        // Set values using closed-form formula
        const double delta0 = val - std::floor(val); // f0
        const double delta1 = std::ceil(val) - val; // 1-f0
        const double inv_delta_sum = 1 / ((1/delta0) + (1/delta1)); // f0 (1-f0) in split case
        v[0][solver->getNumRows()] = inv_delta_sum / delta0;
        v[1][solver->getNumRows()] = inv_delta_sum / delta1;

        // Compute lambda^t_i = (d^t A_N^-1)_i / delta_t and lambda_i = max_t lambda^t_i
        // where the disjunction is (d^0 x >= d^0_0) \vee (d^1 x \ge d^1_0)
        // and here (d^0, d^0_0) = (-e_k, -std::floor(val))
        // and (d^1, d^1_0) = (e_k, std::ceil(val))
        double v_rhs = 1;
        for (int ray_ind = 0; ray_ind < (int) NBVarIndex.size(); ray_ind++) {
          const int NBVar = NBVarIndex[ray_ind];
          const int tmp_row = NBVar - solver->getNumCols();
          double mult = 1.;
          if (NBVar < solver->getNumCols()) {
            mult = isNonBasicUBVar(solver, NBVar) ? -1. : 1.;
          } else {
            mult = (solver->getRowSense()[tmp_row] == 'L') ? -1. : 1.;
          }
          const double lambda0 = -currRay[ray_ind][splitVarRowIndex] / delta0;
          const double lambda1 = currRay[ray_ind][splitVarRowIndex] / delta1;
          const double lambda[2] = { lambda0, lambda1 };
          const int t_i = (lambda0 > lambda1) ? 0 : 1;
          const int v_row = (NBVar < solver->getNumCols()) ? (solver->getNumRows() + 1 + NBVar) : tmp_row;
          for (int t = 0; t < 2; t++) {
            v[t][v_row] = mult * inv_delta_sum * (lambda[t_i] - lambda[t]);
          }

          double curr_rhs = 0.;;
          if (NBVar < solver->getNumCols()) {
            curr_rhs = (mult < 0) ? -1. * solver->getColUpper()[NBVar] : solver->getColLower()[NBVar];
          } else {
            curr_rhs = mult * solver->getRightHandSide()[tmp_row];
          }
          v_rhs += lambda[t_i] * curr_rhs;
        }

        v_rhs *= inv_delta_sum;

        // Scale v as appropriate based on intCut rhs vs calculated rhs
        const double scale = std::abs(intCut.rhs() / v_rhs);
        for (int t = 0; t < 2; t++) {
          for (auto& vti : v[t]) {
            vti *= scale;
          }
        }

        const CoinPackedVector lhs = intCut.row();
#ifdef TRACE
        for (int t = 0; t < 2; t++) {
          SplitDisjunction disj;
          disj.var = var;
          disj.prepareDisjunction(solver);

          OsiSolverInterface* termSolver;
          disj.getSolverForTerm(termSolver, t, solver, false, DIFFEPS, logfile);
          std::vector<double> new_coeff(solver->getNumCols());
          getCutFromCertificate(new_coeff, v[t], termSolver);

          std::vector<double> cut_coeff(solver->getNumCols(), 0.0);

          const int num_el = lhs.getNumElements();
          for (int i = 0; i < num_el; i++) {
            cut_coeff[lhs.getIndices()[i]] = lhs.getElements()[i];
          }

          int num_errors = 0;
          double total_diff = 0.;
          for (int i = 0; i < solver->getNumCols(); i++) {
            const double diff = cut_coeff[i] - new_coeff[i];
            if (greaterThanVal(std::abs(diff), 0.0)) {
              fprintf(stderr, "%d: cut: %.6f\tcalc: %.6f\tdiff: %g\n", i, cut_coeff[i], new_coeff[i], diff);
              num_errors++;
              total_diff += std::abs(diff);
            }
          }
          if (num_errors > 0) printf("Number of differences between true and calculated cuts: %d. Total difference: %g.\n", num_errors, total_diff);

          if (termSolver) { delete termSolver; }
        }
#endif

        // Do strengthening; for this we first need to setup a split disjunction
        SplitDisjunction disj;
        disj.var = var;
        disj.prepareDisjunction(solver);
        const double rhs = intCut.rhs();
        std::vector<double> str_coeff;
        double str_rhs;
        const int curr_num_coeffs_str = strengthenCut(str_coeff, str_rhs, lhs.getNumElements(), lhs.getIndices(), lhs.getElements(), rhs, &disj, v, solver, logfile);

        // Update stats
        num_cuts_strengthened += (curr_num_coeffs_str > 0);
        updateStatsBeforeFinalize(num_coeff_str_stats, curr_num_coeffs_str);

        // Replace row
        if (curr_num_coeffs_str > 0) {
          CoinPackedVector strCutCoeff(str_coeff.size(), str_coeff.data());
          intCut.setRow(strCutCoeff);
          intCut.setLb(str_rhs);
        }
      } // do strengthening

      // Insert new cut into currGMICs
      currGMICs.insert(intCut);
    } // iterate over cols, generating GMICs

    finalizeStats(num_coeff_str_stats, currGMICs.sizeCuts());
  } // Generate GMICs via createMIG and apply closed-form strengthening

  fprintf(stdout, "Finished generating %d GMICs with mode %d and strengthen option %d", currGMICs.sizeCuts(), option, strengthen_option);
  if (strengthen_option != 0) { fprintf(stdout, " (%d / %d cuts strengthened)", num_cuts_strengthened, currGMICs.sizeCuts()); }
  fprintf(stdout, ".\n");

  if (strengthen_option != 0) {
    fprintf(stdout, "Number coeffs changed:\n");
    printStats(num_coeff_str_stats, true, '\n', stdout);
  }
  fprintf(stdout, "--------------------------------------------------\n");
#if 0
  fprintf(stdout, "\n## Printing GMICs ##\n");
  for (int cut_ind = 0; cut_ind < currGMICs.sizeCuts(); cut_ind++) {
    printf("## Cut %d ##\n", cut_ind);
    const OsiRowCut* const cut = currGMICs.rowCutPtr(cut_ind);
    //cut->print();
    const double rhs_mult = (isInfinity(std::abs(cut->lb()))) ? -1. : 1.;
    double mult = rhs_mult;
    if (!isZero(cut->rhs())) {
      mult /= std::abs(cut->rhs());
    }
    const CoinPackedVector row = cut->row();
    const int num_elem = row.getNumElements();
    const int* ind = row.getIndices();
    const double* coeff = row.getElements();
    for (int i = 0; i < num_elem; i++) {
      fprintf(stdout, "(%d, %.6g)\n", ind[i], mult * coeff[i]);
    }
    fprintf(stdout, "(unnormalized) rhs: %g\n\n", rhs_mult * cut->rhs());
    //if (cut_ind > -1) exit(1);
  }
  fprintf(stdout, "Finished printing GMICs.\n\n");
#endif
} /* generateGomoryCuts */

void createMIG(
    /// [out] Cut to generate
    OsiRowCut &cut, 
    /// [in] Solver from which we generate cuts; assumed optimal
    const OsiSolverInterface* const solver,
    /// [in] Index of the variable to split
    const int splitVarIndex,
    /// [in] Index of the row in which the splitVarIndex is basic
    const int splitVarRowIndex,
    /// [in] Whether to strengthen the cut
    const bool strengthen,
    /// [in] Where to write log messages
    FILE* logfile) {
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
    // Skip basic variables (be wary of difference between basis_ and isBasicVar)
    // 2024-06-21: sometimes, such as with bc_presolved, var is basic with basis_,
    // but not isBasicVar (maybe due to free variables)
    // if (isBasicVar(solver, var)) { continue; }
    const int row = var - solver->getNumCols();
    const CoinWarmStartBasis::Status status = (row < 0) ? basis_->getStructStatus(var) : basis_->getArtifStatus(row);
    if (status == CoinWarmStartBasis::basic) {
        continue;
    }

    double value = 0.0;
    if (var < solver->getNumCols()) {  
      if (status == CoinWarmStartBasis::atUpperBound) { // TODO I think this might be wrong for ub vars
        value = -intersectionCutCoeff(-basisRowStruct[var], f0, solver, var, strengthen);
        cutRhs += value * colUpper[var];
      } else if (status == CoinWarmStartBasis::atLowerBound) {
        value = intersectionCutCoeff(basisRowStruct[var], f0, solver, var, strengthen);
        cutRhs += value * colLower[var];
      } else {
        error_msg(errorstring, "createMIG: Invalid basis\n");
        writeErrorToLog(errorstring, logfile);
        if (mustDelete && basis_) {
          delete basis_;
        }
        exit(1); // probably better to throw than exit
      }
    } else {
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
            (status == CoinWarmStartBasis::atUpperBound)
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

double unstrengthenedIntersectionCutCoeff(double abar, double f0) {
  if (abar > 0) {
    //return alpha_i * (1 - beta);
    return (abar / f0);
  } else {
    //return -alpha_i * beta;
    return (-abar / (1 - f0));
  }
} /* unstrengthenedIntersectionCutCoeff */

double modularizedCoeff(double abar, double f0) {
  double f_i = abar - floor(abar);
  if (f_i <= f0) {
    return f_i;
  } else {
    return f_i - 1;
  }
} /* modularizedCoeff */

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
