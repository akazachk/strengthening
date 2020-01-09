/**
 * @file CglAdvCut.cpp
 * @author A. M. Kazachkov
 * @date 2018-12-24
 */
#include "CglAdvCut.hpp"

#include <cmath> // abs, floor, ceil
#include <limits> // numeric_limits
#include <algorithm> // std::min_element, std::max_element

// COIN-OR files
#include <CbcModel.hpp>
#include <CglGMI.hpp>

// Project files
#include "SolverHelper.hpp"
#include "SolverInterface.hpp"
#include "utility.hpp"

// VPC files
#include "CglVPC.hpp"

using namespace StrengtheningParameters;

#ifdef TRACE
#include "debug.hpp"
#endif

const std::vector<std::string> CglAdvCut::ExitReasonName {
  "SUCCESS",
  "CUT_LIMIT",
  "FAIL_LIMIT",
  "TIME_LIMIT",
  "NO_CUTS_LIKELY",
  "UNKNOWN"
}; /* ExitReasonName */
const std::vector<std::string> CglAdvCut::CutTimeStatsName {
  "TOTAL_TIME",
  "INIT_SOLVE_TIME",
  "GEN_CUTS_TIME"
}; /* CutTimeStatsName */
const std::vector<std::string> CglAdvCut::CutTypeName {
  "ONE_SIDED_CUT", "OPTIMALITY_CUT", "OTHER_CUT"
}; /* CutTypeName */
const std::vector<std::string> CglAdvCut::ObjectiveTypeName {
  "DUMMY_OBJ",
  "OBJ_CUT",
  "ONE_SIDED",
  "OTHER"
}; /* ObjectiveTypeName */
const std::vector<std::string> CglAdvCut::FailureTypeName {
  "ABANDONED",
  "BAD_DYNAMISM",
  "BAD_SUPPORT",
  "BAD_VIOLATION",
  "CUT_LIMIT",
  "DUAL_INFEASIBLE",
  "DUPLICATE_SIC",
  "DUPLICATE_OTHER",
  "ITERATION_LIMIT",
  "ORTHOGONALITY_SIC",
  "ORTHOGONALITY_OTHER",
  "PRIMAL_INFEASIBLE",
  "TIME_LIMIT",
  "NUMERICAL_ISSUES_WARNING",
  "DLB_EQUALS_DUB_NO_OBJ",
  "DLB_EQUALS_LPOPT_NO_OBJ",
  "PRIMAL_INFEASIBLE_NO_OBJ",
  "NUMERICAL_ISSUES_NO_OBJ",
  "UNKNOWN"
}; /* FailureTypeName */

/**
 * @brief Universal way to check whether we reached the limit for the number of cuts for each split
 * This allows us to change between restricting number of cuts per split and total number of cuts easily
 */
int CglAdvCut::getCutLimit(const int CUTLIMIT, const int numFracVar) {
  // The cut limit is either across all cut-generating sets
  // or it is divided among the cut-generating sets (either as the limit / numFracVars, or as a fixed number per cgs)
  // If CUTLIMIT = 0 => no cuts
  // If CUTLIMIT > 0 => absolute cut limit
  // If CUTLIMIT < 0 => cut limit per cgs
  if (CUTLIMIT == 0) {
    return 0;
    //return std::numeric_limits<int>::max();
  } else if (CUTLIMIT > 0) {
    return CUTLIMIT;
  } else {
    return (numFracVar <= 0) ? 0 : std::ceil((-1. * CUTLIMIT * numFracVar));
  }
} /* getCutLimit */

/** getCutLimit */
int CglAdvCut::getCutLimit() const {
  return params.get(CUTLIMIT);
} /* getCutLimit */

/** Default constructor */
CglAdvCut::CglAdvCut() {
  initialize();
} /* default constructor */

/** Param constructor */
CglAdvCut::CglAdvCut(const Parameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy constructor */
CglAdvCut::CglAdvCut(const CglAdvCut& source) : CglCutGenerator(source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
CglAdvCut::~CglAdvCut() {
  // Delete anything owned by this class
} /* destructor */

/** Assignment operator */
CglAdvCut& CglAdvCut::operator=(const CglAdvCut& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Clone */
CglCutGenerator* CglAdvCut::clone() const {
  return new CglAdvCut(*this);
} /* clone */

/** setParams */
void CglAdvCut::setParams(const Parameters& param) {
  this->params = param;
} /* setParams */

/**
 * @brief Generate VPCs from a disjunction (e.g., arising from a partial branch-and-bound tree)
 */
void CglAdvCut::generateCuts(const OsiSolverInterface& si, OsiCuts& cuts, const CglTreeInfo info) {
  CglAdvCut::ExitReason status = CglAdvCut::ExitReason::UNKNOWN;
  if (reachedTimeLimit(CutTimeStats::TOTAL_TIME, params.get(TIMELIMIT))) {
    status = CglAdvCut::ExitReason::TIME_LIMIT_EXIT;
    finish(status);
    return;
  }

  // Reset things in preparation for round of cuts, in case a previous round was done using this generator
  setupAsNew();
  init_num_cuts = cuts.sizeCuts();
  if (init_num_cuts == 0) {
    this->cutType.resize(0);
    this->objType.resize(0);
  }
  else if (this->canReplaceGivenCuts) {
    // If we are going to be able to replace given cuts,
    // then it must be the case that the current cutType and cutHeurVec
    // should correspond to those old cuts
    const int num_old_cut_type = cutType.size();
    if (num_old_cut_type != init_num_cuts) {
      error_msg(errorstring,
          "# given cuts: %d. # old cuts: %d. "
          "Either set canReplaceGivenCuts to false, or ensure that old cuts are "
          "accounted for in the cutType and cutHeurVec members.\n",
          init_num_cuts, num_old_cut_type);
      writeErrorToLog(errorstring, params.logfile);
      exit(1);
    }
  }

  // Set cut limit
  params.set(CUTLIMIT,
      CglAdvCut::getCutLimit(params.get(CUTLIMIT),
          si.getFractionalIndices().size()));
  if (reachedCutLimit()) {
    if (si.getFractionalIndices().size() > 0) {
      status = CglAdvCut::ExitReason::CUT_LIMIT_EXIT;
    } else {
      status = CglAdvCut::ExitReason::NO_CUTS_LIKELY_EXIT;
    }
    finish(status);
    return;
  }

  // Solver can't be const because custom enableFactorization function might do initialSolve or resolve
  SolverInterface* solver;
  try {
    solver = const_cast<SolverInterface*>(dynamic_cast<const SolverInterface*>(&si));
  } catch (std::exception& e) {
    error_msg(errorstring, "Unable to cast as SolverInterface.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
  if (!solver->isProvenOptimal()) {
    error_msg(errorstring, "CglAdvCut::generateCuts: Solver not proven optimal.\n");
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Time starts here, and will end when finish is called
  timer.start_timer(CutTimeStatsName[static_cast<int>(CutTimeStats::TOTAL_TIME)]);

  // Make a copy of the solver to allow for fixing variables and changed bounds at root
  SolverInterface* mysolver = dynamic_cast<SolverInterface*>(si.clone());

  //==================================================//
  // Set up VPC generation
  VPCParametersNamespace::VPCParameters vpc_params;
  vpc_params.set(VPCParametersNamespace::DISJ_TERMS, this->params.get(StrengtheningParameters::DISJ_TERMS));
  vpc_params.set(VPCParametersNamespace::CUTLIMIT, this->params.get(StrengtheningParameters::CUTLIMIT));

  this->gen.setParams(vpc_params);

  OsiCuts tmpcuts = cuts;
  this->gen.generateCuts(*mysolver, tmpcuts);

  // Save stats
  this->timer = this->gen.timer;
  this->num_obj_tried = this->gen.num_obj_tried;
  this->num_failures = this->gen.num_failures;

  // Save cut type // TODO check this works with rounds
  for (int cut_ind = 0; cut_ind < tmpcuts.sizeCuts(); cut_ind++) {
    CglVPC::CutType vpcCutType = this->gen.cutType[cut_ind];
    CglAdvCut::CutType currCutType;
    if (vpcCutType == CglVPC::CutType::ONE_SIDED_CUT) {
      currCutType = CglAdvCut::CutType::ONE_SIDED_CUT;
    }
    else if (vpcCutType == CglVPC::CutType::OPTIMALITY_CUT) {
      currCutType = CglAdvCut::CutType::OPTIMALITY_CUT;
    }
    else {
      currCutType = CglAdvCut::CutType::OTHER_CUT;
    }

    CglVPC::ObjectiveType vpcObjType = this->gen.objType[cut_ind];
    CglAdvCut::ObjectiveType currObjType;
    if (vpcObjType == CglVPC::ObjectiveType::DUMMY_OBJ) {
      currObjType = CglAdvCut::ObjectiveType::DUMMY_OBJ;
    }
    else if (vpcObjType == CglVPC::ObjectiveType::ONE_SIDED) {
      currObjType = CglAdvCut::ObjectiveType::ONE_SIDED;
    }
    else {
      currObjType = CglAdvCut::ObjectiveType::OTHER;
    }

    // Set cutType, numCutsOfType, objType, numCutsFromHeur, num_cuts
    this->addCut(tmpcuts.rowCut(cut_ind), cuts, currCutType, currObjType);

  } // iterate over cuts

  // Set numFails (if all failure types match, then we can use same counting; otherwise we need to be more careful but for now we will just lump everything together)
  if (static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES) == static_cast<int>(CglVPC::FailureType::NUM_FAILURE_TYPES)) {
    this->numFails = this->gen.numFails;
  } else {
    this->numFails[static_cast<int>(CglAdvCut::FailureType::UNKNOWN)] = this->num_failures;
  }

  // Similarly for setting numObjFromHeur, numFailsFromHeur
  if (static_cast<int>(CglAdvCut::ObjectiveType::NUM_OBJECTIVE_TYPES) == static_cast<int>(CglVPC::ObjectiveType::NUM_OBJECTIVE_TYPES)) {
    this->numObjFromHeur = this->gen.numObjFromHeur;
    this->numFailsFromHeur = this->gen.numFailsFromHeur;
  } else {
    this->numObjFromHeur[static_cast<int>(CglAdvCut::ObjectiveType::OTHER)] = this->num_obj_tried;
    this->numFailsFromHeur[static_cast<int>(CglAdvCut::ObjectiveType::OTHER)] = this->num_failures;
  }

  // Save exit reason
  CglVPC::ExitReason vpcExitReason = this->gen.exitReason;
  if (vpcExitReason == CglVPC::ExitReason::SUCCESS_EXIT)
    status = CglAdvCut::ExitReason::SUCCESS_EXIT;
  else if (vpcExitReason == CglVPC::ExitReason::CUT_LIMIT_EXIT)
    status = CglAdvCut::ExitReason::CUT_LIMIT_EXIT;
  else if (vpcExitReason == CglVPC::ExitReason::FAIL_LIMIT_EXIT)
    status = CglAdvCut::ExitReason::FAIL_LIMIT_EXIT;
  else if (vpcExitReason == CglVPC::ExitReason::TIME_LIMIT_EXIT)
    status = CglAdvCut::ExitReason::TIME_LIMIT_EXIT;
  else if (vpcExitReason == CglVPC::ExitReason::NO_CUTS_LIKELY_EXIT)
    status = CglAdvCut::ExitReason::NO_CUTS_LIKELY_EXIT;
  else
    status = CglAdvCut::ExitReason::UNKNOWN;
  
  if (mysolver) { delete mysolver; }
  finish(status);
} /* generateCuts */

/**
 * @brief Add cut and update statistics about cut type, objective type that led to the cut, total number of cuts, ...
 */
void CglAdvCut::addCut(const OsiRowCut& cut, OsiCuts& cuts, const CutType& type,
    const ObjectiveType& cutHeur) {
  cuts.insert(cut);
  cutType.push_back(type);
  numCutsOfType[static_cast<int>(type)]++;
  objType.push_back(cutHeur);
  numCutsFromHeur[static_cast<int>(cutHeur)]++;
  num_cuts++;
} /* addCut */

/****************** PROTECTED **********************/

/**
 * @brief Reset _some_ things (those corresponding to a previous run of this generator)
 * E.g., we do not reset timing, the cutType vector, or the cutHeurVec
 * The latter two should not be changed and need to correspond to the cuts passed into generateCuts
 */
void CglAdvCut::setupAsNew() {
  this->exitReason = CglAdvCut::ExitReason::UNKNOWN;
  this->numCutsOfType.clear();
  this->numCutsOfType.resize(static_cast<int>(CutType::NUM_CUT_TYPES), 0);
  this->numCutsFromHeur.clear();
  this->numCutsFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numObjFromHeur.clear();
  this->numObjFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numFailsFromHeur.clear();
  this->numFailsFromHeur.resize(static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  this->numFails.clear();
  this->numFails.resize(static_cast<int>(FailureType::NUM_FAILURE_TYPES), 0);
  this->init_num_cuts = 0;
  this->num_cuts = 0;
  this->num_obj_tried = 0;
  this->num_failures = 0;
  this->probData.EPS = this->params.get(EPS);
  if (!this->isSetupForRepeatedUse) {
    this->canReplaceGivenCuts = false;
    this->cutType.resize(0);
    this->objType.resize(0);
  }
} /* setupAsNew */

/**
 * @brief Initialize the cut generator internal data
 */
void CglAdvCut::initialize(const CglAdvCut* const source, const Parameters* const param) {
  if (param != NULL)
    setParams(*param);
  if (source != NULL) {
    if (param == NULL)
      setParams(source->params);
    this->canReplaceGivenCuts = source->canReplaceGivenCuts;
    this->isSetupForRepeatedUse = source->isSetupForRepeatedUse;
    this->exitReason = source->exitReason;
    this->timer = source->timer;
    this->cutType = source->cutType;
    this->objType = source->objType;
    this->numCutsOfType = source->numCutsOfType;
    this->numCutsFromHeur = source->numCutsFromHeur;
    this->numObjFromHeur = source->numObjFromHeur;
    this->numFailsFromHeur = source->numFailsFromHeur;
    this->numFails = source->numFails;
    this->init_num_cuts = source->init_num_cuts;
    this->num_cuts = source->num_cuts;
    this->num_obj_tried = source->num_obj_tried;
    this->num_failures = source->num_failures;
    this->probData = source->probData;
  }
  else {
    this->isSetupForRepeatedUse = false;
    for (int t = 0; t < static_cast<int>(CutTimeStats::NUM_TIME_STATS); t++) {
      timer.register_name(CutTimeStatsName[t]);
    }
    for (int t = 0; t < static_cast<int>(ObjectiveType::NUM_OBJECTIVE_TYPES); t++) {
      timer.register_name(ObjectiveTypeName[t] + "_TIME");
    }
    setupAsNew();
  }
} /* initialize */

/**
 * @brief Get data about solver and optimal basis
 *
 * Get problem data such as min/max coeff, problem-specific epsilon,
 * nonbasic variables, row in which each variable is basic, etc.
 */
void CglAdvCut::getProblemData(
    /// [in/out] Solver being used to determine the nonbasic space; note that the basis and/or solution may change due to enableFactorization
    OsiSolverInterface* const solver,
    /// [out] Where to save the data
    ProblemData& probData, 
    /// [in] If this is a subproblem, then we may want to pass the original problem data (to save the locations of the original nonbasic variables)
    const ProblemData* const origProbData,
    /// [in] Whether to enable factorization (can change solution slightly)
    const bool enable_factorization) {
  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();

  if (enable_factorization)
    enableFactorization(solver, params.get(doubleParam::EPS)); // this may change the solution slightly

  // Set min/max reference values based on problem data
  double minReferenceValue = 1, maxReferenceValue = 1;
  const double* elements = solver->getMatrixByCol()->getElements();
  probData.minAbsCoeff = *std::min_element(elements,
      elements + solver->getMatrixByCol()->getNumElements(),
      [](double i, double j) {return std::abs(i) < std::abs(j);});
  probData.maxAbsCoeff = *std::max_element(elements,
        elements + solver->getMatrixByCol()->getNumElements(),
        [](double i, double j) {return std::abs(i) < std::abs(j);});
  probData.minAbsCoeff = std::abs(probData.minAbsCoeff);
  probData.maxAbsCoeff = std::abs(probData.maxAbsCoeff);
  minReferenceValue =
      (probData.minAbsCoeff > 0) ?
          CoinMin(minReferenceValue, probData.minAbsCoeff) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, probData.maxAbsCoeff);

  double lp_opt = solver->getObjValue();
  minReferenceValue =
      (lp_opt != 0) ?
          CoinMin(minReferenceValue, std::abs(lp_opt)) : minReferenceValue;
  maxReferenceValue = CoinMax(maxReferenceValue, std::abs(lp_opt));

  // Prepare data structures
  probData.num_cols = numCols;
  probData.lp_opt = lp_opt;
  probData.NBVarIndex.clear();
  probData.NBVarIndex.reserve(numCols);
  probData.rowOfVar.clear();
  probData.rowOfVar.resize(numCols + numRows, -1);
  probData.varBasicInRow.resize(numRows);
  solver->getBasics(&probData.varBasicInRow[0]); // get which variable is basic in each row
  if (origProbData) {
    probData.rowOfOrigNBVar.clear();
    probData.rowOfOrigNBVar.resize(numCols, -1);
  }

  // Get the basic and nonbasic original variable info
  // Note that fixed nonbasic variables might not need to be added to the nonbasic index vector
  // since they do not correspond to any rays...
  // but right now we typically assume we get n rays in some parts of the code
  for (int var = 0; var < numCols + numRows; var++) {
    // Update min/max reference values
    if (var < numCols) {
      // Count how many non-inf lower and upper bounds
      double colLB, colUB;
      colLB = solver->getColLower()[var];
      colUB = solver->getColUpper()[var];
      if (colLB > -1 * solver->getInfinity() + params.get(doubleParam::EPS)) {
        minReferenceValue =
            (colLB != 0.) ?
                CoinMin(minReferenceValue, std::abs(colLB)) : minReferenceValue;
        maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colLB));
      }

      if (colUB < solver->getInfinity() - params.get(doubleParam::EPS)) {
        minReferenceValue =
            (colUB != 0.) ?
                CoinMin(minReferenceValue, std::abs(colUB)) : minReferenceValue;
        maxReferenceValue = CoinMax(maxReferenceValue, std::abs(colUB));
      }
    } else {
      // Should only store the nonbasic slacks corresponding to inequalities since
      // the equality slack columns don't correspond to any rays
      // But, just as with fixed variables, we don't currently do this
      const int row = var - numCols;
      const double absrhs = std::abs(solver->getRightHandSide()[row]);
      minReferenceValue =
          (absrhs > 0) ? CoinMin(minReferenceValue, absrhs) : minReferenceValue;
      maxReferenceValue = CoinMax(maxReferenceValue, absrhs);

      // Assign var basic in this row to rowOfVar
      const int basic_var = probData.varBasicInRow[row];
      probData.rowOfVar[basic_var] = row;

      // If this basic_var was non-basic in the basis at v,
      // we need to add this row in the right spot in rowOfOrigNBVar
      if (origProbData) {
        const int NBIndexOfBasicVar = origProbData->getVarNBIndex(basic_var);
        if (NBIndexOfBasicVar >= 0) {
          probData.rowOfOrigNBVar[NBIndexOfBasicVar] = row;
        }
      }
    }

    if (!isBasicVar(solver, var)) {
      // Recall that rowOfVar stores -1 - nb index for nb variables
      // The -1 is to prevent the conflict of the 0th nb var and var basic in row 0
      probData.rowOfVar[var] -= probData.NBVarIndex.size();
      probData.NBVarIndex.push_back(var);
    } else {
      // Quick check that basic slack vars are basic in their row
      const int row = var - numCols;
      if (row >= 0) {
        if (probData.varBasicInRow[row] != var) {
          // Elsewhere we are using that each slack is basic in its own row,
          // so if this is not true, we will have to adjust
          // We use this, for example, in PCut, to calculate the right order
          // for the packed NB rays in genCornerNB
          error_msg(errstr,
              "Basic variable in row %d is variable %d, but it should be the slack on this row (%d).\n",
              row, probData.varBasicInRow[row], var);
          writeErrorToLog(errstr, params.logfile);
          exit(1);
        }
      }
    }
  } // loop over variables

  // May also need to save rays of C1, where the coefficients of inv(B) * A
  // are sometimes negated because we have
  // \bar x = inv(B) * b - inv(B) * A * x_N
  const int numNB = probData.NBVarIndex.size();
  int tempIndex = 0;
  probData.NBReducedCost.resize(numNB);
  for (int j = 0; j < numNB; j++) {
    const int NBVar = probData.NBVarIndex[j];

    if (NBVar < numCols) {
      if (isNonBasicUBCol(solver, NBVar))
        probData.NBReducedCost[j] = -1.0 * solver->getReducedCost()[NBVar];
      else
        probData.NBReducedCost[j] = solver->getReducedCost()[NBVar];
    } else {
      tempIndex = NBVar - numCols;
      if (isNonBasicUBSlack(solver, tempIndex))
        probData.NBReducedCost[j] = solver->getRowPrice()[tempIndex];
      else
        probData.NBReducedCost[j] = -1.0 * solver->getRowPrice()[tempIndex];
    }
    if (lessThanVal(probData.NBReducedCost[j], 0.)) {
      if (lessThanVal(probData.NBReducedCost[j], 0., params.get(doubleConst::DIFFEPS))) {
        error_msg(errorstring,
            "Nonbasic reduced cost should be >= 0 in the complemented nonbasic space. "
            "However, for nb var %d (real var %d), it is %e.\n",
            j, NBVar, probData.NBReducedCost[j]);
        writeErrorToLog(errorstring, params.logfile);
        exit(1);
      } else {
        warning_msg(warnstring,
            "Nonbasic reduced cost should be >= 0 in the complemented nonbasic space. "
            "However, for nb var %d (real var %d), it is %e. Small enough error that we only send warning.\n",
            j, NBVar, probData.NBReducedCost[j]);
        this->numFails[static_cast<int>(FailureType::NUMERICAL_ISSUES_WARNING)]++;
      }
    }

//    if (isZero(nonBasicReducedCost[j], EPS)) {
//      numDualDegeneratePivots++;
//    }
  } // loop over nonbasic vars to collect reduced costs

  // Set data-specific epsilon
  probData.EPS = CoinMin(params.get(doubleParam::EPS), minReferenceValue / maxReferenceValue);

  if (enable_factorization)
    solver->disableFactorization();
} /* getProblemData */

/**
 * @brief Checks whether too many unsuccessful objective attempts have been made
 *
 * There are four types of checks. The first three have non-decreasing success requirements.
 * 1. Few cuts have been generated:
 *      default is FEW_CUTS = 1 cut, and the success threshold is at least 1 cut every 20 obj (fail ratio = .95).
 * 2. Many cuts have been generated:
 *      default is MANY_CUTS = .25 * CUT_LIMIT, and the success threshold is at least 1 cut every 10 obj (fail ratio = .90).
 * 3. Many obj have been tried, i.e., we have been trying for a long time so we better be successful super often:
 *      MANY_OBJ = max(FEW_CUTS / (1-few_cuts_fail_threshold), MANY_CUTS / (1-many_cuts_fail_threshold));
 *      default = max(20, 2.5 * CUT_LIMIT) and the success threshold is at least 1 cut every 5 obj (fail ratio = .80).
 * 4. Time is too long and we are not too successful:
 *       \# obj tried >= MANY_OBJ && time >= 10 && average time / obj >= CUTSOLVER_TIMELIMIT + 1
 *       the success threshold is at least 1 cut every 3 obj
 *
 * Examples:
 * If the failure ratio is .85, then this will return true only if # obj >= MANY_OBJ.
 * This means that, as long as we have not tried too many times unsuccessfully,
 * then even if we have MANY_CUTS cuts, we would like more, and we feel that we are doing pretty well generating them.
 * (Of course, we have not hit the cut limit yet.)
 *
 * If your current failure ratio (# fails / # obj) is .93, then this will return true if # cuts >= MANY_CUTS, as .93 > .90.
 * Note that if # cuts < MANY_CUTS, but # obj >= MANY_OBJ,
 * then we would have failed earlier based on the more stringent limit for the MANY_OBJ case.
 *
 * If the failure ratio is .99, then we got here with # cuts < MANY_CUTS and # obj < MANY_OBJ.
 * Otherwise, if # cuts >= MANY_CUTS, then we would have hit the failure limit earlier (the first time it went above .90).
 * Similarly, if # obj >= MANY_OBJ, then we would have hit the limit even earlier.
 *
 * If the failure ratio is 1 (all failures), then we will reject if # obj > FEW_CUTS / (1.-few_cuts_fail_threshold).
 */
bool CglAdvCut::reachedFailureLimit(const int num_cuts, const int num_fails, //const double time,
    const double few_cuts_fail_threshold, const double many_cuts_fail_threshold,
    const double many_obj_fail_threshold, const double time_fail_threshold) const {
  const int num_obj_tried = num_cuts + num_fails;
  if (num_obj_tried == 0) {
    return false;
  }
  const int CUT_LIMIT = getCutLimit();
  const int FEW_CUTS = 1;
  const int NO_CUTS_OBJ_LIMIT = std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold));
  const int MANY_CUTS = std::ceil(.25 * CUT_LIMIT);
  int MANY_OBJ = CoinMax(
      std::ceil(FEW_CUTS / (1. - few_cuts_fail_threshold)),
      std::ceil(MANY_CUTS / (1. - many_cuts_fail_threshold)));
  if (MANY_OBJ < 0) // in case of overflow
    MANY_OBJ = std::numeric_limits<int>::max();
  const double fail_ratio = (double) num_fails / num_obj_tried;
  bool reached_limit = false;
  if (num_obj_tried >= MANY_OBJ && greaterThanVal(fail_ratio, many_obj_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts >= MANY_CUTS && greaterThanVal(fail_ratio, many_cuts_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts >= FEW_CUTS && greaterThanVal(fail_ratio, few_cuts_fail_threshold)) {
    reached_limit = true;
  } else if (num_cuts < FEW_CUTS && num_obj_tried >= NO_CUTS_OBJ_LIMIT) {
    reached_limit = true;
  }
  const double time = 0; //TODO Add here if desired, e.g., timer.get_total_time(CutTimeStatsName[static_cast<int>(CutTimeStats::TOTAL_TIME)]);
  const double max_avg_time = 1000; // e.g., before it was params.get(PRLP_TIMELIMIT)
  if (!reached_limit && num_obj_tried >= MANY_OBJ && time > 10.
      && time / num_obj_tried >= max_avg_time) { // checks if average time is too high
    reached_limit = true;
  }
//  if (reached_limit) {
//    this->exitReason = CglAdvCut::ExitReason::FAIL_LIMIT_EXIT;
//  }
  return reached_limit;
} /* reachedFailureLimit */
