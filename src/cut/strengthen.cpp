/**
 * @file strengthen.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#include "strengthen.hpp"

// COIN-OR files
#include <OsiRowCut.hpp> // for setting a disjunctive term
#include <OsiCuts.hpp>
#include <CbcModel.hpp>

// Project files
#include "CutHelper.hpp" // for setObjective and addToObjective
#include "Disjunction.hpp"
//#include "Parameters.hpp"
#include "SolverInterface.hpp"
#include "SolverHelper.hpp" // isBasicVar, checkSolverOptimality, isNonBasicUBVar
#include "utility.hpp"
#include "verify.hpp"

#ifdef USE_EIGEN
#include "eigen.hpp"
#endif

/// @details Given problem is (with A an [m x n] matrix)
///   Ax - s = b
///   x >= 0
///   s >= 0.
/// Suppose the basis is B (m columns); note that some of these might be slack variables.
/// Denote by N the remaining n columns (the cobasis).
///
/// If solver is proven optimal, and the cut is valid for the basis cone,
/// then it suffices to take the inverse of the optimal basis.
/// In particular, we want alpha^T = v^T A.
/// 
/// Note that alpha is n-dimensional, v is m-dimensional.
/// The components of v corresponding to rows in which the slack is basic will be 0.
/// We need to get access to the inverse of the nxn submatrix of A corresponding to nonbasic slacks.
///
/// * Invertible case
/// If A is invertible,
/// then we can take v^T = alpha^T A^{-1}
/// or, equivalently, v = (A^{-1})^T alpha.
/// This means that v_j can be calcuated as the dot product of the j-th column of A^{-1} with alpha.
///
/// * General case
/// If A is not invertible, then we need only consider the relaxed corner polyhedron.
/// That is, we only need the basis inverse applied to the matrix.
void getCertificate(
    /// [out] Farkas multipliers (vector of length m + m_t + n)
    TermCutCertificate& v,
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] LP solver corresponding to disjunctive term
    OsiSolverInterface* const term_solver,
    /// [in] logfile for error printing
    FILE* logfile) {
  v.clear();
  v.resize(term_solver->getNumCols() + term_solver->getNumRows(), 0.0);

  // Get dense cut coefficients so that we can set the objective vector
  std::vector<double> cut_coeff(term_solver->getNumCols(), 0.0);
  for (int i = 0; i < num_elem; i++) {
    cut_coeff[ind[i]] = coeff[i];
  }
  term_solver->setObjective(cut_coeff.data());
  term_solver->resolve();

  if (!term_solver->isProvenOptimal()) {
    return;
  }

  term_solver->enableFactorization();

  // Collect nonbasic variables
  // TODO Can / should we do this with getBInv calls instead?
  std::vector<int> rows, cols;
  //std::vector<int> NBVarIndex;
  //NBVarIndex.reserve(term_solver->getNumCols());
  for (int var = 0; var < term_solver->getNumCols() + term_solver->getNumRows(); var++) {
    if (!isBasicVar(term_solver, var)) {
      //NBVarIndex.push_back(var);
      if (var < term_solver->getNumCols()) {
        cols.push_back(var);
      } else {
        rows.push_back(var - term_solver->getNumCols());
      }
    }
  }

#ifdef USE_EIGEN
  calculateCertificateEigen(v, num_elem, ind, coeff, term_solver, rows, cols, logfile);
#else
  error_msg(
      errorstring,
      "Calculating certificate without Eigen is not yet implemented. Exiting.\n");
  writeErrorToLog(errorstring, logfile);
  throw std::logic_error(errorstring);
#endif // USE_EIGEN

  term_solver->disableFactorization();

  // Verify certificate results in same cut
  checkCutHelper(cut_coeff, v, -1, NULL, term_solver, logfile);
} /* getCertificate */

void getCertificateForTerm(
    /// [out] Farkas multipliers
    TermCutCertificate& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] disjunctive term specification
    const DisjunctiveTerm* const term,
    /// [in] tolerance for reconstructing disjunctive term solution
    const double DIFFEPS,
    /// [in] logfile for error printing
    FILE* logfile) {
  error_msg(errorstring, "Currently implemented incorrectly (missing common_changed_var effect).\n");
  writeErrorToLog(errorstring, logfile);
  exit(1);

  OsiSolverInterface* termSolver = solver->clone();
  const int curr_num_changed_bounds = term->changed_var.size();
  std::vector < std::vector<int> > commonTermIndices(curr_num_changed_bounds);
  std::vector < std::vector<double> > commonTermCoeff(curr_num_changed_bounds);
  std::vector<double> commonTermRHS(curr_num_changed_bounds);

  for (int i = 0; i < curr_num_changed_bounds; i++) {
    const int col = term->changed_var[i];
    const double coeff = (term->changed_bound[i] <= 0) ? 1. : -1.;
    const double val = term->changed_value[i];
    commonTermIndices[i].resize(1, col);
    commonTermCoeff[i].resize(1, coeff);
    commonTermRHS[i] = coeff * val;
    if (term->changed_bound[i] <= 0) {
      termSolver->setColLower(col, val);
    } else {
      termSolver->setColUpper(col, val);
    }
  }

  const int curr_num_added_ineqs = term->ineqs.size();
  for (int i = 0; i < curr_num_added_ineqs; i++) {
    const OsiRowCut* currCut = &term->ineqs[i];
    termSolver->applyRowCuts(1, currCut); // hopefully this works
  }

  // Set the warm start
  if (term->basis && !(termSolver->setWarmStart(term->basis))) {
    error_msg(errorstring,
        "Warm start information not accepted for term.\n");
    writeErrorToLog(errorstring, logfile);
    exit(1);
  }

  // Resolve and check the objective matches
#ifdef TRACE
  printf("\n## Solving for term ##\n");
#endif
  termSolver->resolve();
  const bool calcAndFeasTerm = checkSolverOptimality(termSolver, true);

  // If something went wrong
  if (!calcAndFeasTerm) {
    printf("\n## Term is not proven optimal. Exiting from this term. ##\n");
    delete termSolver;
    return;
  }

  // Sometimes we run into a few issues getting the ``right'' value
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    termSolver->resolve();
  }
  if (!isVal(termSolver->getObjValue(), term->obj, DIFFEPS)) {
    double ratio = termSolver->getObjValue() / term->obj;
    if (ratio < 1.) {
      ratio = 1. / ratio;
    }
    // Allow it to be up to 3% off without causing an error
    if (greaterThanVal(ratio, 1.03)) {
      error_msg(errorstring,
          "Objective at disjunctive term is incorrect. Before, it was %s, now it is %s.\n",
          stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
      writeErrorToLog(errorstring, logfile);
      exit(1);
    } else {
      warning_msg(warnstring,
          "Objective at disjunctive term is incorrect. Before, it was %s, now it is %s.\n",
          stringValue(term->obj, "%1.3f").c_str(),
          stringValue(termSolver->getObjValue(), "%1.3f").c_str());
    }
#ifdef TRACE
    std::string commonName;
    Disjunction::setCgsName(commonName, curr_num_changed_bounds, commonTermIndices,
        commonTermCoeff, commonTermRHS, false);
    printf("Bounds changed: %s.\n", commonName.c_str());
#endif
  } // check that objective value matches

  getCertificate(v, num_elem, ind, coeff, termSolver, logfile);
  if (termSolver) { delete termSolver; }
} /* getCertificateForTerm */

/// Use Theorem B.5 of Kazachkov 2018 dissertation to get certificate with trivial normalization
void getCertificateTrivial(
    /// [out] Farkas multipliers
    TermCutCertificate& v, 
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    // [in] disjunction
    const Disjunction* const disj,
    // [in] initial LP solver
    const OsiSolverInterface* const solver) {
  // TODO
} /* getCertificateTrivial */

/**
 * @brief Attempt to strengthen coefficients of given cuts
 *
 * Required:
 * (1) original cuts
 * (2) disjunction for which the cut is valid
 * (3) Farkas certificate for the cut for this disjunction, 
 *      where the first m rows correspond to the original constraint matrix, 
 *      the next m_t rows correspond to the multipliers on the disjunctive term inequalities
 *      the last n rows are bounds on the variables, 
 * (4) globally valid lower bounds on the disjunctive term inequalties
 *
 * With those, the strengthening requires finding values m_1,...,m_T
 *  str_coeff[k] = coeff[k] + max_t { -u^t_k + u^t_0 (D^t_0 - \ell^t) m_t }
 *
 * @return Number of coefficients strengthened
 *
 * TODO If k is nonbasic at a nonzero bound, what do we do?
 * TODO Get this to work with general disjunctive inequalities, not just bound changes (see below TODO)
 */
int strengthenCutS(
    /// [out] strengthened cuts
    OsiCuts& str_cuts,
    /// [in] original cuts
    const OsiCuts& cuts,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const CutCertificate& v, 
    /// [in] original solver (used to get globally-valid lower bounds for the disjunctive terms)
    const OsiSolverInterface* const solver,
    /// [in] logfile for error printing
    FILE* logfile,
    /// [in] IP solution to original problem (will usually be empty unless you are debugging)
    const std::vector<double>& ip_solution) {
  const int num_cuts = cuts.sizeCuts();
  if (num_cuts == 0) { return 0; }
  if (!disj) { return 0; }
  const int num_terms = disj->num_terms;
  if (num_terms < 2) { return 0; }
  for (int i = 0; i < num_cuts; i++) {
    str_cuts.insert(cuts.rowCut(i));
  }

  // Get globally lower bound for each of the disjunctive terms
  // and use this to set up term-by-term coefficients for strengthening
  // TODO get this to work with generic inequality description of disjunction, not just the bound changes version
  std::vector<std::vector<double> > disj_lb_diff(num_terms);
  std::vector<double> lb_term(num_terms, 0.0); // u^t_0 (D^t_0 - \ell^t)
  for (int term_ind = 0; term_ind < num_terms; term_ind++) {
    const DisjunctiveTerm& term = disj->terms[term_ind];
    const int num_common = (int) disj->common_changed_var.size() + disj->common_ineqs.size();
    const int num_disj_ineqs = num_common + term.changed_var.size() + term.ineqs.size();
    // Resize current term vector
    disj_lb_diff[term_ind].resize(num_disj_ineqs, 0.0);

    for (int bound_ind = 0; bound_ind < num_disj_ineqs; bound_ind++) {
      const std::vector<int>& changed_var = (bound_ind < num_common) ? disj->common_changed_var : term.changed_var;
      const std::vector<int>& changed_bound = (bound_ind < num_common) ? disj->common_changed_bound : term.changed_bound;
      const std::vector<double>& changed_value = (bound_ind < num_common) ? disj->common_changed_value : term.changed_value;
      const int real_bound_ind = (bound_ind < num_common) ? bound_ind : bound_ind - num_common;
      const int var = changed_var[real_bound_ind];
      const double mult = (changed_bound[real_bound_ind] <= 0) ? 1.0 : -1.0;
      const double bd = (changed_bound[real_bound_ind] <= 0) ? solver->getColLower()[var] : solver->getColUpper()[var];
      const double val = changed_value[real_bound_ind];
      double lb = mult * bd;
      if (isInfinity(std::abs(lb))) {
        {
          // No bound readily available; try to optimize to find one
          OsiSolverInterface* tmpSolver = solver->clone();
          addToObjectiveFromPackedVector(tmpSolver, NULL, true);
          setConstantObjectiveFromPackedVector(tmpSolver, mult, 1, &var);
          tmpSolver->resolve();
          checkSolverOptimality(tmpSolver, false);
          if (tmpSolver->isProvenOptimal()) {
            lb = tmpSolver->getObjValue();
          }

          if (tmpSolver) { delete tmpSolver; }
        }
        if (isInfinity(std::abs(lb))) {
          warning_msg(warnstring,
              "Cannot strengthen cut using this disjunction because variable %d is missing its %s bound.\n",
              var, mult > 0 ? "lower" : "upper");
          return 0;
        }
      } // check if lb is infinite

      disj_lb_diff[term_ind][bound_ind] = mult * val - lb;
#ifdef TRACE
      assert(disj_lb_diff[term_ind][bound_ind] > -1e-3);
#endif
      if (!isZero(v[term_ind][solver->getNumRows() + bound_ind]) && !isZero(disj_lb_diff[term_ind][bound_ind])) {
        lb_term[term_ind] += v[term_ind][solver->getNumRows() + bound_ind] * disj_lb_diff[term_ind][bound_ind];
        if (isZero(lb_term[term_ind])) {
          lb_term[term_ind] = 0.;
        }
      }
    } // loop over bounds
  } // loop over terms

  // Strengthen the coefficient on each cut
  SolverInterface* mono = NULL;
  int num_coeffs_changed = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
    if (!solver->isInteger(col)) continue;

    // Check if both lower and upper bound have positive multipliers across all terms
    int lt_zero_ind = -1, gt_zero_ind = -1;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const DisjunctiveTerm& term = disj->terms[term_ind];
      const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
      const double ukt = v[term_ind][solver->getNumRows() + num_disj_ineqs + col];
      if (lessThanVal(ukt, 0)) {
        if (lt_zero_ind == -1) lt_zero_ind = term_ind;
      } else if (greaterThanVal(ukt, 0)) {
        if (gt_zero_ind == -1) gt_zero_ind = term_ind;
      }
      if (lt_zero_ind != -1 && gt_zero_ind != -1) break;
    } // loop over terms
    const double mult = (gt_zero_ind == -1) ? -1 : 1.;

    if (!mono && num_terms > 2) {
      mono = new SolverInterface;
      setupMonoidalIP(mono, col, disj, lb_term, v, solver, mult);
    } else if (mono != NULL) {
      updateMonoidalIP(mono, col, disj, v, solver, mult);
    }

    // Now try to strengthen the coefficients on the integer-restricted variables
    for (int cut_ind = 0; cut_ind < num_cuts; cut_ind++) {
      OsiRowCut& cut = str_cuts.rowCut(cut_ind);
      double str_rhs = cut.rhs();

      CoinPackedVector& row = cut.mutableRow();
      const int num_el = row.getNumElements();
      int* ind = row.getIndices();
      double* el = row.getElements();
      int ind_of_col = 0;
      for (ind_of_col = 0; ind_of_col < num_el; ind_of_col++) {
        if (ind[ind_of_col] >= col) {
          break;
        }
      }

      // Check if this index already existed in the cut
      if (ind_of_col < num_el && ind[ind_of_col] == col) {
        double& str_coeff = el[ind_of_col];
        num_coeffs_changed += strengthenCutCoefficient(str_coeff, str_rhs, col, str_coeff, disj, lb_term, v, solver, mono, logfile);
      } else {
        // The coefficient on col is currently 0
        double str_coeff = 0.;
        num_coeffs_changed += strengthenCutCoefficient(str_coeff, str_rhs, col, str_coeff, disj, lb_term, v, solver, mono, logfile);
        if (!isZero(str_coeff)) {
          // We need to insert this into the cut
          row.insert(ind_of_col, str_coeff);
        }
      }

      cut.setLb(str_rhs);
    } // loop over cuts
  } // loop over cols

  if (mono) { delete mono; }

  return num_coeffs_changed;
} /* strengthenCutS */

/**
 * @details Required:
 * (1) original cut
 * (2) disjunction for which the cut is valid
 * (3) Farkas certificate for the cut for this disjunction, 
 *      where the first m rows correspond to the original constraint matrix, 
 *      the next m_t rows correspond to the multipliers on the disjunctive term inequalities
 *      the last n rows are bounds on the variables, 
 * (4) globally valid lower bounds on the disjunctive term inequalties
 *
 * With those, the strengthening requires finding values m_1,...,m_T
 *  str_coeff[k] = coeff[k] + max_t { -u^t_k + u^t_0 (D^t_0 - \ell^t) m_t }
 *
 * @return Number of coefficients strengthened
 *
 * TODO If k is nonbasic at a nonzero bound, what do we do?
 * TODO Get this to work with general disjunctive inequalities, not just bound changes (see below TODO)
 */
int strengthenCut(
    /// [out] strengthened cut coefficients
    std::vector<double>& str_coeff,
    /// [out] strengthened cut rhs
    double& str_rhs,
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    /// [in] original cut rhs
    const double rhs,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + num_disj_ineqs + n)
    const CutCertificate& v, 
    /// [in] original solver (used to get globally-valid lower bounds for the disjunctive terms)
    const OsiSolverInterface* const solver,
    /// [in] logfile for error printing
    FILE* logfile,
    /// [in] IP solution to original problem (will usually be empty unless you are debugging)
    const std::vector<double>& ip_solution) {
  // Set up original cut coeff and rhs
  str_rhs = rhs;
  str_coeff.clear();
  str_coeff.resize(solver->getNumCols(), 0.0);
  for (int i = 0; i < num_elem; i++) {
    str_coeff[ind[i]] = coeff[i];
  }

  if (!disj) { return 0; }

  // Get globally-valid lower bound for each of the disjunctive terms
  // and use this to set up term-by-term coefficients for strengthening
  // TODO get this to work with generic inequality description of disjunction, not just the bound changes version
  std::vector<std::vector<double> > disj_lb_diff(disj->num_terms); // D^t_0 - ell^t
  std::vector<double> lb_term(disj->num_terms, 0.0); // u^t_0 (D^t_0 - ell^t)
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const DisjunctiveTerm& term = disj->terms[term_ind];

    if (!term.is_feasible) {
      continue;
    }

    const int num_common = (int) disj->common_changed_var.size() + disj->common_ineqs.size();
    const int num_disj_ineqs = num_common + term.changed_var.size() + term.ineqs.size();
    // Resize current term vector
    disj_lb_diff[term_ind].resize(num_disj_ineqs, 0.0);

    for (int bound_ind = 0; bound_ind < num_disj_ineqs; bound_ind++) {
      const std::vector<int>& changed_var = (bound_ind < num_common) ? disj->common_changed_var : term.changed_var;
      const std::vector<int>& changed_bound = (bound_ind < num_common) ? disj->common_changed_bound : term.changed_bound;
      const std::vector<double>& changed_value = (bound_ind < num_common) ? disj->common_changed_value : term.changed_value;
      const int real_bound_ind = (bound_ind < num_common) ? bound_ind : bound_ind - num_common;
      const int var = changed_var[real_bound_ind];
      const double mult = (changed_bound[real_bound_ind] <= 0) ? 1.0 : -1.0;
      const double bd = (changed_bound[real_bound_ind] <= 0) ? solver->getColLower()[var] : solver->getColUpper()[var];
      const double val = changed_value[real_bound_ind];
      double lb = mult * bd;
      if (isInfinity(std::abs(lb))) {
        {
          // No bound readily available; try to optimize to find one
          OsiSolverInterface* tmpSolver = solver->clone();
          addToObjectiveFromPackedVector(tmpSolver, NULL, true);
          setConstantObjectiveFromPackedVector(tmpSolver, mult, 1, &var);
          tmpSolver->resolve();
          checkSolverOptimality(tmpSolver, false);
          if (tmpSolver->isProvenOptimal()) {
            lb = tmpSolver->getObjValue();
          }

          if (tmpSolver) { delete tmpSolver; }
        }
        if (isInfinity(std::abs(lb))) {
          warning_msg(warnstring,
              "Cannot strengthen cut using this disjunction because variable %d is missing its %s bound.\n",
              var, mult > 0 ? "lower" : "upper");
          return 0;
        }
      } // check if lb is infinite

      disj_lb_diff[term_ind][bound_ind] = mult * val - lb;
#ifdef TRACE
      assert(disj_lb_diff[term_ind][bound_ind] > -1e-3);
#endif
      if (!isZero(v[term_ind][solver->getNumRows() + bound_ind]) && !isZero(disj_lb_diff[term_ind][bound_ind])) {
        lb_term[term_ind] += v[term_ind][solver->getNumRows() + bound_ind] * disj_lb_diff[term_ind][bound_ind];
        if (isZero(lb_term[term_ind])) {
          lb_term[term_ind] = 0.;
        }
      }
    } // loop over bounds
  } // loop over terms

  // Setup monoidal cut strengthening LP
  SolverInterface* mono = NULL;

  // Now try to strengthen the coefficients on the integer-restricted variables
  int num_coeffs_changed = 0;
  for (int col = 0; col < solver->getNumCols(); col++) {
    if (!solver->isInteger(col)) continue;

    // Check if both lower and upper bound have positive multipliers across all terms
    int lt_zero_ind = -1, gt_zero_ind = -1;
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const DisjunctiveTerm& term = disj->terms[term_ind];
      const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
      const double ukt = v[term_ind][solver->getNumRows() + num_disj_ineqs + col];
      if (lessThanVal(ukt, 0)) {
        if (lt_zero_ind == -1) lt_zero_ind = term_ind;
      } else if (greaterThanVal(ukt, 0)) {
        if (gt_zero_ind == -1) gt_zero_ind = term_ind;
      }
      if (lt_zero_ind != -1 && gt_zero_ind != -1) break;
    } // loop over terms
    const double mult = (gt_zero_ind == -1) ? -1 : 1.;

    if (!mono && disj->num_terms > 2) {
      mono = new SolverInterface;
      setupMonoidalIP(mono, col, disj, lb_term, v, solver, mult);
    } else if (mono != NULL) {
      updateMonoidalIP(mono, col, disj, v, solver, mult);
    }

    if (strengthenCutCoefficient(str_coeff[col], str_rhs, col, str_coeff[col], disj, lb_term, v, solver, mono, logfile)) {
      num_coeffs_changed += 1;

      if (!ip_solution.empty()) {
        const double activity = dotProduct(str_coeff.data(), ip_solution.data(), str_coeff.size());
        if (lessThanVal(activity, str_rhs, 1e-3)) {
          error_msg(errorstring, "Cut removes optimal solution after strengthening col %d. Activity: %.10f. Rhs: %.10f.\n", col, activity, str_rhs);
          writeErrorToLog(errorstring, logfile);
          exit(1);
        }
        if (lessThanVal(activity, str_rhs)) {
          warning_msg(warnstring, "Cut removes optimal solution after strengthening col %d. Activity: %.10f. Rhs: %.10f.\n", col, activity, str_rhs);
        }
      }
    } // if coefficient strengthened, updated num_coeffs_changed and check continued correctness
  } // loop over cols

  if (mono) { delete mono; }

  return num_coeffs_changed;
} /* strengthenCut */

/// Calculate new cut coefficient (we need to minimize over monoids m)
/// \return Whether the coefficient has changed
bool strengthenCutCoefficient(
    /// [out] strengthened cut coeff
    double& str_coeff,
    /// [out] strengthened cut rhs
    double& str_rhs,
    /// [in] variable for which the coefficient is being changed
    const int var,
    /// [in] original coeff
    const double coeff,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] u^t_0 (D^t_0 - ell^t) for all t
    const std::vector<double>& lb_term,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length num_rows + num_disj_term_ineqs + num_cols)
    const CutCertificate& v, 
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in/out] monoidal cut strengthening solver (req'd for num_terms > 2),
    OsiSolverInterface* const mono,
    /// [in] logfile for error printing
    FILE* logfile) {
  const int num_terms = disj->num_terms;
  if (num_terms < 2) {
    str_coeff = coeff;
    return false;
  }

  // Check if complementing is necessary (translating may be necessary too)
  // Complementing happens if the multiplier is on the bound on var is negative, which means x_k <= bound is in effect
  // Theory assumes that x_k >= 0 is the bound, so generically strengthening a cut
  //   alpha x >= beta
  // (in which x_k <= g_k is being used)
  // is equivalent to substituting y_k = g_k - x_k and strengthening
  //   alpha_{-k} x_{-k} + alpha_k (g_k - y_k) >= beta
  // or (moving constant terms to the right-hand side)
  //  alpha_{-k} x_{-k} - alpha_k y_k >= beta - alpha_k g_k.
  // Once we have a strengthened coefficient for x_k, say gamma, we need to uncomplement to get back to the x-space:
  //  alpha_{-k} x_{-k} - gamma_k (g_k - x_k) >= beta - alpha_k g_k
  //  ==>  alpha_{-k} x_{-k} + gamma_k x_k >= beta + (gamma_k - alpha_k) g_k.
  //
  // An analogous situation happens when we have x_k >= ell_k (not equal to 0)
  // In this case, the cut the right-hand side of the strengthened cut becomes
  //   beta + (gamma_k - alpha_k) ell_k.
  int lt_zero_ind = -1, gt_zero_ind = -1;
  for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
    const DisjunctiveTerm& term = disj->terms[term_ind];
    const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
    const double ukt = v[term_ind][solver->getNumRows() + num_disj_ineqs + var];
    if (lessThanVal(ukt, 0)) {
      if (lt_zero_ind == -1) lt_zero_ind = term_ind;
    } else if (greaterThanVal(ukt, 0)) {
      if (gt_zero_ind == -1) gt_zero_ind = term_ind;
    }
    if (lt_zero_ind != -1 && gt_zero_ind != -1) break;
  } // loop over terms
  const double mult = (gt_zero_ind == -1) ? -1 : 1.;
  // Can it happen that one of the multipliers is on the lower bound, and one is on the upper bound?
  if (lt_zero_ind != -1 && gt_zero_ind != -1) {
    const int num_common = (int) disj->common_changed_var.size() + disj->common_ineqs.size();
    const double uk0 = v[lt_zero_ind][solver->getNumRows() + num_common + disj->terms[lt_zero_ind].changed_var.size() + var];
    const double uk1 = v[gt_zero_ind][solver->getNumRows() + num_common + disj->terms[gt_zero_ind].changed_var.size() + var];
    if (!isZero(uk0) && !isZero(uk1)) {
      warning_msg(warnstring,
          "CHECK: The u^t_k multipliers on variable k = %d are of different signs: u^{%d}_k = %.6e, u^{%d}_k = %.6e."
          " This is strange and may not be handled correctly in the code.\n",
          var, lt_zero_ind, uk0, gt_zero_ind, uk1);
      //str_coeff = coeff;
      //return false;
    }
  }
  str_coeff = mult * coeff;

  if (num_terms == 2 && !mono) {
    std::vector<int> m1(num_terms, 0.0);
    std::vector<int> m2(num_terms, 0.0);
    std::vector<int> m3(num_terms, 0.0);
    std::vector<int> m4(num_terms, 0.0);
    std::vector<int> m5(num_terms, 0.0);
    m2[0] = 1.0;
    m2[1] = -1.0;
    m3[0] = -1.0;
    m3[1] = 1.0;
    m4[0] = -2.0;
    m4[1] = 2.0;
    m5[0] = 2.0;
    m5[1] = -2.0;
    std::vector<std::vector<int> > m_options = { m1, m2, m3, m4, m5 };
    double min_max_term_val = std::numeric_limits<double>::max();
    for (const auto& m : m_options) {
      double max_term_val = std::numeric_limits<double>::lowest(); 
      for (int term_ind = 0; term_ind < num_terms; term_ind++) {
        const DisjunctiveTerm& term = disj->terms[term_ind];
        const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
        const double utk = v[term_ind][solver->getNumRows() + num_disj_ineqs + var];
        // TODO figure out safe floating point comparison
        const double curr_val = ((mult > 0 && utk < 0) ? 0. : -utk * mult) + lb_term[term_ind] * m[term_ind];
        if (curr_val > max_term_val) {
          max_term_val = curr_val;
        }
      }
      if (max_term_val < min_max_term_val) {
        min_max_term_val = max_term_val;
      }
    } // loop over monoid options

    if (!isZero(min_max_term_val) && !isInfinity(min_max_term_val)) {
      str_coeff += min_max_term_val;
    }
  } // if mono solver not given
  else if (mono) {
    CbcModel model(*mono);
    setIPSolverParameters(&model, mono->messageHandler()->logLevel());
    model.setModelOwnsSolver(false);
    model.branchAndBound();
    if (model.status() != 0) {
      error_msg(errorstring, "Failed to optimize monoidal strengthening solver for var %d.\n", var);
      writeErrorToLog(errorstring, logfile);
      exit(1);
    }
    str_coeff += model.getObjValue();
  } // num_terms > 2
  else {
    throw("Unable to strengthen coefficient.\n");
  }

  str_coeff *= mult; // 2024-06-06 fixed bug: previously was done after checking isVal, overcounting # str coeffs
  if (!isVal(str_coeff, coeff)) {
    // If complemented, adjust right-hand side
    if (mult < 0) {
      if (!isZero(solver->getColUpper()[var])) {
        str_rhs += (str_coeff - coeff) * solver->getColUpper()[var];
      }
    }

    // Translate if the lower bound on var is nonzero
    if (mult > 0 && !isZero(solver->getColLower()[var])) {
      str_rhs += (str_coeff - coeff) * solver->getColLower()[var];
    }
    return true;
  } // check if strengthening happened
  else {
    str_coeff = coeff;
    return false;
  }
} /* strengthenCutCoefficient */

/**
 * @brief Creates monoidal strengthening IP
 *
 * @details Let d_t := u^t_0 (D^t_0 - \ell^t)
 *
 * min_{\gamma_k, m \in \Z^T} \gamma_k + \alpha_k (we leave out \alpha_k in the formulation as it is a constant)
 *    - d_t m_t + gamma \ge -u^t_k    for all t
 *   \sum_t m_t         \ge 0
 *          m_t         \in \Z        for all t (we should add bounds to make it easier)
 */
void setupMonoidalIP(
    /// [in/out] monoidal IP solver (memory needs to already have been allocated)
    OsiSolverInterface* const mono,
    /// [in] variable for which the coefficient is being changed
    const int var,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] u^t_0 (D^t_0 - \ell^t) for all t
    const std::vector<double>& lb_term,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const CutCertificate& v, 
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] multiplier on utk (negative 1 if the upper bound is being used for all terms, so we will later complement)
    const double mult) {
  if (!disj) return;
  if (!mono) {
    fprintf(stderr,
        "*** ERROR: setupMonoidalIP: Monoidal cut strengthening LP mono needs to have been allocated.\n");
    exit(1);
  }
  const int num_terms = disj->num_terms;
  const int num_cols = num_terms + 1; // m_t for all t; gamma 
  const int num_rows = num_terms + 1; // \sum_t m_t \ge 0; gamma - d_t m_t \ge -u^t_k
  const int numElementsEstimate = 3 * num_terms;

  CoinPackedMatrix mx; // create a new col-ordered matrix
  mx.reverseOrdering(); // make it row-ordered
  mx.setDimensions(0, num_cols);
  mx.reserve(num_rows, numElementsEstimate, false);

  std::vector<int> rowStarts(num_rows);
  std::vector<int> column(numElementsEstimate);
  std::vector<double> element(numElementsEstimate);
  std::vector<double> rowLB(num_rows, 0.);

  int el_ind = 0;
  int row_ind = 0;

  // -d_t m_t + gamma \ge -u^t_k
  for (int t = 0; t < num_terms; t++) {
    rowStarts[row_ind] = el_ind;
    column[el_ind] = t;
    element[el_ind] = -1. * lb_term[t];
    el_ind++;
    column[el_ind] = num_terms;
    element[el_ind] = 1.;
    el_ind++;

    const DisjunctiveTerm& term = disj->terms[t];
    const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
    //const double utk = std::abs(v[t][solver->getNumRows() + num_disj_ineqs + var]);
    const double utk = v[t][solver->getNumRows() + num_disj_ineqs + var];
    // TODO figure out safe floating point comparison
    const double curr_val = ((mult > 0 && utk < 0) ? 0. : -utk * mult);
    rowLB[row_ind] = curr_val;

    row_ind++;
  } // loop over terms

  // \sum_t m_t \ge 0
  rowStarts[row_ind] = el_ind;
  for (int t = 0; t < num_terms; t++) {
    column[el_ind] = t;
    element[el_ind] = 1.;
    el_ind++;
  }
  rowLB[row_ind] = 0;

  // Finish setting up matrix
  //mx.appendRows(num_rows, rowStarts.data(), column.data(), element.data(), num_cols);
  for (int row = 0; row < num_rows; row++) {
    const int vecsize = (row < num_rows-1) ? 2 : num_terms;
    const int offset = rowStarts[row];
    mx.appendRow(vecsize, column.data() + offset, element.data() + offset);
  }

  // Column bounds
  std::vector<double> colLB(num_cols, -1 * solver->getInfinity());
  std::vector<double> colUB(num_cols, solver->getInfinity());

  // Objective
  std::vector<double> obj(num_cols, 0.);
  obj[num_cols - 1] = 1.;

  // Load problem
  // Defaults (in COIN-OR):
  // colLB: 0 **** We want it to be -inf
  // colUB: inf
  // obj coeff: 0 **** We set it to be minimize gamma
  // rowSense: >=
  // rhs: 0 **** We want this to be as above
  // rowRange: 0
  mono->loadProblem(mx, colLB.data(), colUB.data(), obj.data(), NULL,
      rowLB.data(), NULL);
  mono->disableFactorization();

  for (int t = 0; t < num_terms; t++) {
    mono->setInteger(t);
  }

  // Set message handling
  setLPSolverParameters(mono, solver->messageHandler()->logLevel());
} /* setupMonoidalIP */

void updateMonoidalIP(
    /// [in/out] monoidal solver (assumed to already be set up)
    OsiSolverInterface* const mono,
    /// [in] variable for which we are changing the coefficients
    const int var,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const CutCertificate& v, 
    /// [in] original solver
    const OsiSolverInterface* const solver,
    /// [in] multiplier on utk (negative 1 if the upper bound is being used for all terms, so we will later complement)
    const double mult) {
  if (!disj) { return; }
  for (int t = 0; t < disj->num_terms; t++) {
    const DisjunctiveTerm& term = disj->terms[t];
    const int num_disj_ineqs = (int) disj->common_changed_var.size() + disj->common_ineqs.size() + term.changed_var.size() + term.ineqs.size();
    //const double utk = std::abs(v[t][solver->getNumRows() + num_disj_ineqs + var]);
    const double utk = v[t][solver->getNumRows() + num_disj_ineqs + var];
    // TODO figure out safe floating point comparison
    const double curr_val = ((mult > 0 && utk < 0) ? 0. : -utk * mult);
    mono->setRowLower(t, curr_val);
  } // loop over terms
} /* updateMonoidalIP */

/**
 * @brief Calculate value of m to use in strengthening
 *
 * This algorithm assumes u^t_0 D^t_0 = 1 + u^t_0 \ell^t > 0
 * (e.g., satisfied in the nonbasic space for disjunctions involving only binary variables)
 * This leads to the cuts alpha x >= 1 (which we strengthen with the below algorithm)
 * 
 * Algorithm 2 (Balas and Jeroslow, 1979)
 * -----------
 * \gamma^t_0 = u^t_0 D^t_0
 * \gamma^t_k = u^t_0 D^t_{\cdot, k}
 * \gamma_k \gets (\sum_t \gamma^t_k) / (\sum_t \gamma^t_0) for each k \in I
 *                \max_t \gamma^t_k / \gamma^t_0 for all k \notin I
 * m^*_t \gets \gamma^t_0 \gamma_k - \gamma^t_k (implies (\gamma^t_k + m^*_t) / \gamma^t_0 = \gamma_k)
 * m_t \gets floor(m^*_t)
 * for i = 1 to -\sum_t m_t
 *   m_t \gets m_t + 1 for t = argmin_t (\gamma^t_k + m_t + 1) / \gamma^t_0
 * \gamma_k = max_t (\gamma^t_k + m_t) / \gamma^t_0
 *
 * @return Number of coefficients strengthened
 */
int BalJer79_Algorithm2(
    /// [out] strengthened cut coefficients
    std::vector<double>& str_coeff,
    /// [out] strengthened cut rhs
    double& str_rhs,
    /// [in] number of nonzero cut coefficients
    const int num_elem, 
    /// [in] indices of nonzero cut coefficients
    const int* const ind, 
    /// [in] nonzero cut coefficients
    const double* const coeff,
    /// [in] original cut rhs
    const double rhs,
    /// [in] disjunction
    const Disjunction* const disj,
    /// [in] Farkas multipliers for each of the terms of the disjunction (each is a vector of length m + n)
    const CutCertificate& v, 
    /// [in] original solver
    const OsiSolverInterface* const solver) {
  // Set up original cut coeff and rhs
  str_rhs = rhs;
  str_coeff.clear();
  str_coeff.resize(solver->getNumCols(), 0.0);
  for (int i = 0; i < num_elem; i++) {
    str_coeff[ind[i]] = coeff[i];
  }

  if (!disj) { return 0; }

  const int num_terms = disj->num_terms;

  std::vector<double> gamma_t0(num_terms, 0.);
  std::vector<double> gamma_tk(num_terms, 0.);
  double sum_gamma_t0 = 0.;
  //double sum_gamma_tk = 0.;
  for (int t = 0; t < num_terms; t++) {
    const DisjunctiveTerm& term = disj->terms[t];
    const int num_disj_ineqs = (int) disj->common_changed_var.size() + term.changed_var.size();

    // TODO get this to work with general inequalities not just bound changes
    for (int bound_ind = 0; bound_ind < num_disj_ineqs; bound_ind++) {
      const double mult = (term.changed_bound[bound_ind] <= 0) ? 1.0 : -1.0;
      gamma_t0[t] += v[t][solver->getNumRows() + bound_ind] * mult * term.changed_value[bound_ind];
      gamma_tk[t] += v[t][solver->getNumRows() + bound_ind] * mult;
      sum_gamma_t0 += gamma_t0[t];
      //sum_gamma_tk += gamma_tk[t];
    }
  } // loop over terms to set up gamma_t0 and gamma_tk

  if (isZero(sum_gamma_t0)) {
    fprintf(stderr,
        "\\sum_t \\gamma^t_0 = 0, which it should not be.\n");
    exit(1);
  }

  /*double gamma_k = 0.;
  if (solver->isInteger(var)) {
    gamma_k = sum_gamma_tk / sum_gamma_t0;
  } else {
    int max_ind = -1;
  }*/
  fprintf(stderr,
      "*** ERROR: BalJer79_Alg2 implementation unfinished.\n");
  exit(1);
} /* BalJer79_Algorithm2 */
