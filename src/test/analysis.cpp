/**
 * @file analysis.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-24
 */
#include "analysis.hpp"

// COIN-OR
#include <OsiSolverInterface.hpp>
#include <OsiCuts.hpp>
#include <CglGMI.hpp>

// Project files
#include "CglAdvCut.hpp"
#include "CutHelper.hpp" // isRowDifferent
#include "BBHelper.hpp"
#include "Disjunction.hpp"
#include "SolverHelper.hpp"
#include "Parameters.hpp"
using namespace StrengtheningParameters;
#include "utility.hpp" // isInfinity, stringValue

template <typename T>
std::vector<double> computeStats(const std::vector<T>& vec) {
  std::vector<double> stats(static_cast<int>(Stat::num_stats), 0.);
  if (vec.empty()) {
    return stats;
  }
  stats[static_cast<int>(Stat::total)] = 0.;
  stats[static_cast<int>(Stat::min)] = std::numeric_limits<double>::max();
  stats[static_cast<int>(Stat::max)] = std::numeric_limits<double>::lowest();
  for (const auto& val : vec) {
    stats[static_cast<int>(Stat::total)] += val;
    stats[static_cast<int>(Stat::min)] = std::min(stats[static_cast<int>(Stat::min)], static_cast<double>(val));
    stats[static_cast<int>(Stat::max)] = std::max(stats[static_cast<int>(Stat::max)], static_cast<double>(val));
  }
  stats[static_cast<int>(Stat::avg)] = stats[static_cast<int>(Stat::total)] / static_cast<double>(vec.size());
  for (const auto& val : vec) {
    stats[static_cast<int>(Stat::stddev)] += (val - stats[static_cast<int>(Stat::avg)]) * (val - stats[static_cast<int>(Stat::avg)]);
  }
  stats[static_cast<int>(Stat::stddev)] = sqrt(stats[static_cast<int>(Stat::stddev)] / static_cast<double>(vec.size()));
  return stats;
} /* computeStats */

// Below values are used to doublecheck that the columns in the header are correctly counted
const int countBoundInfoEntries = 20;
const int countGapInfoEntries = 9;
const int countSummaryBBInfoEntries = 4 * 2;
const int countFullBBInfoEntries = static_cast<int>(BB_INFO_CONTENTS.size()) * 4 * 2;
const int countOrigProbEntries = 13;
const int countPostCutProbEntries = 10;
const int countDisjInfoEntries = 12;
const int countCutInfoEntries = 13;
const int countObjInfoEntries = 1;
const int countFailInfoEntries = 1 + static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES);
const std::vector<std::string> STR_INFO_CONTENTS = { "NUM STR AFFECTED CUTS", "NUM COEFFS STR AVG", "NUM COEFFS STR STDDEV", "NUM COEFFS STR MIN", "NUM COEFFS STR MAX" };
const int countStrInfoEntries = STR_INFO_CONTENTS.size() * 2;
const std::vector<std::string> CERT_INFO_CONTENTS = {
  "NUM UNMATCHED BOUNDS",
  "AVG NUM CGS FACETS",
  "NUM IRREG LESS",
  "NUM REGULAR",
  "NUM IRREG MORE",
  "NUM UNCONVERGED",
  "RCVMIP ITER AVG",
  "RCVMIP ITER STDDEV",
  "RCVMIP ITER MIN",
  "RCVMIP ITER MAX",
  "RCVMIP ITER NUM LIMIT",
  "RCVMIP TIME AVG",
  "RCVMIP TIME STDDEV",
  "RCVMIP TIME MIN",
  "RCVMIP TIME MAX",
}; /* CERT_INFO_CONTENTS */
const int countCertInfoEntries = CERT_INFO_CONTENTS.size() * 2;
const int countParamInfoEntries = intParam::NUM_INT_PARAMS + doubleParam::NUM_DOUBLE_PARAMS;
int countTimeInfoEntries = 0; // set in printHeader
const int countVersionInfoEntries = 6;
const int countExtraInfoEntries = 4;

void printHeader(const StrengtheningParameters::Parameters& params,
    const std::vector<std::string>& time_name,
    const char SEP) {
  FILE* logfile = params.logfile;
  if (logfile == NULL)
    return;

  countTimeInfoEntries = time_name.size();

  // First line of the header details the categories of information displayed
  std::string tmpstring = "";
  fprintf(logfile, "%c", SEP);
  fprintf(logfile, "%s", "PARAM INFO");
  tmpstring.assign(countParamInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BOUND INFO");
  tmpstring.assign(countBoundInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "GAP INFO");
  tmpstring.assign(countGapInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "BB INFO");
  tmpstring.assign(countSummaryBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "ORIG PROB");
  tmpstring.assign(countOrigProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "POST-CUT PROB");
  tmpstring.assign(countPostCutProbEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "DISJ INFO");
  tmpstring.assign(countDisjInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "STRENGTHENING INFO");
  tmpstring.assign(countStrInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "CERTIFICATE INFO");
  tmpstring.assign(countCertInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "CUT INFO");
  tmpstring.assign(countCutInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "OBJ INFO");
  tmpstring.assign(countObjInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "FAIL INFO");
  tmpstring.assign(countFailInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "FULL BB INFO");
  tmpstring.assign(countFullBBInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "TIME INFO");
  tmpstring.assign(countTimeInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "VERSION INFO");
  tmpstring.assign(countVersionInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "WRAP UP INFO");
  tmpstring.assign(countExtraInfoEntries, SEP);
  fprintf(logfile, "%s", tmpstring.c_str());
  fprintf(logfile, "%s", "END");
  fprintf(logfile, "\n");

  fprintf(logfile, "%s%c", "INSTANCE", SEP);
  { // PARAM INFO
    printParams(params, logfile, 1); // only int/double param names
  } // PARAM INFO
  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", "LP OBJ", SEP); count++;  // 1
    fprintf(logfile, "%s%c", "BEST DISJ OBJ", SEP); count++; // 2
    fprintf(logfile, "%s%c", "WORST DISJ OBJ", SEP); count++; // 3
    fprintf(logfile, "%s%c", "IP OBJ", SEP); count++; // 4
    fprintf(logfile, "%s%c", "NUM CHANGED ROOT BOUNDS", SEP); count++; // 5
    fprintf(logfile, "%s%c", "ROOT OBJ", SEP); count++; // 6
    fprintf(logfile, "%s%c", "NUM GMIC", SEP); count++; // 7
    fprintf(logfile, "%s%c", "GMIC OBJ", SEP); count++; // 8
    fprintf(logfile, "%s%c", "NUM L&PC", SEP); count++; // 9
    fprintf(logfile, "%s%c", "L&PC OBJ", SEP); count++; // 10
    fprintf(logfile, "%s%c", "NUM MYCUTS", SEP); count++; // 11
    fprintf(logfile, "%s%c", "MYCUTS OBJ", SEP); count++; // 12
    fprintf(logfile, "%s%c", "MYCUTS+GMIC BOUND", SEP); count++; // 13
    fprintf(logfile, "%s%c", "NUM AFFECTED STR MYCUTS", SEP); count++; // 14
    fprintf(logfile, "%s%c", "UNSTR MYCUTS OBJ", SEP); count++; // 15
    fprintf(logfile, "%s%c", "UNSTR MYCUTS+GMIC BOUND", SEP); count++; // 16
    fprintf(logfile, "%s%c", "NUM AFFECTED RCVMIP-STR MYCUTS", SEP); count++; // 17
    fprintf(logfile, "%s%c", "RCVMIP MYCUTS OBJ", SEP); count++; // 18
    fprintf(logfile, "%s%c", "RCVMIP MYCUTS+GMIC BOUND", SEP); count++; // 19
    fprintf(logfile, "%s%c", "RCVMIP ALL BOUND", SEP); count++; // 20
    assert(count == countBoundInfoEntries);
  } // BOUND INFO
  { // GAP INFO
    int count = 0;
    fprintf(logfile, "%s%c", "GMIC % GAP CLOSED", SEP); count++; // 1
    fprintf(logfile, "%s%c", "L&PC % GAP CLOSED", SEP); count++; // 2
    fprintf(logfile, "%s%c", "MYCUTS % GAP CLOSED", SEP); count++; // 3
    fprintf(logfile, "%s%c", "GMIC+MYCUTS % GAP CLOSED", SEP); count++; // 4
    fprintf(logfile, "%s%c", "UNSTR MYCUTS % GAP CLOSED", SEP); count++; // 5
    fprintf(logfile, "%s%c", "UNSTR GMIC+MYCUTS % GAP CLOSED", SEP); count++; // 6
    fprintf(logfile, "%s%c", "RCVMIP MYCUTS % GAP CLOSED", SEP); count++; // 7
    fprintf(logfile, "%s%c", "RCVMIP GMIC+MYCUTS % GAP CLOSED", SEP); count++; // 8
    fprintf(logfile, "%s%c", "RCVMIP ALL % GAP CLOSED", SEP); count++; // 9
    assert(count == countGapInfoEntries);
  } // GAP INFO
  { // BB INFO
    int count = 0;
    std::vector<std::string> nameVec = {"NODES", "TIME"};
    for (auto name : nameVec) {
      fprintf(logfile, "%s%c", ("FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("FIRST GUR+CUTS " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR+CUTS " + name).c_str(), SEP); count++;
    }
    assert(count == countSummaryBBInfoEntries);
  } // BB INFO
  { // ORIG PROB
    int count = 0;
    fprintf(logfile, "%s%c", "ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "COLS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "EQ ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "BOUND ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "ASSIGN ROWS", SEP); count++;
    fprintf(logfile, "%s%c", "FIXED COLS", SEP); count++;
    fprintf(logfile, "%s%c", "GEN INT", SEP); count++;
    fprintf(logfile, "%s%c", "BINARY", SEP); count++;
    fprintf(logfile, "%s%c", "CONTINUOUS", SEP); count++;
    fprintf(logfile, "%s%c", "A-DENSITY", SEP); count++;
    assert(count == countOrigProbEntries);
  } // ORIG PROB
  { // POST-CUT PROB
    int count = 0;
    fprintf(logfile, "%s%c", "NEW NUM FRAC", SEP); count++;
    fprintf(logfile, "%s%c", "NEW MIN FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "NEW MAX FRACTIONALITY", SEP); count++;
    fprintf(logfile, "%s%c", "NEW A-DENSITY", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC (gmics)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE MYCUTS (gmics)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC (my cuts)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE MYCUTS (my cuts)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE GMIC (all cuts)", SEP); count++;
    fprintf(logfile, "%s%c", "ACTIVE MYCUTS (all cuts)", SEP); count++;
    assert(count == countPostCutProbEntries);
  } // POST-CUT PROB
  { // DISJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM DISJ TERMS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM INTEGER SOL", SEP); count++;
    fprintf(logfile, "%s%c", "NUM DISJ", SEP); count++;
    fprintf(logfile, "%s%c", "AVG DENSITY PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG ROWS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG COLS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG POINTS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG RAYS PRLP", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB EXPLORED NODES", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB PRUNED NODES", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB MIN DEPTH", SEP); count++;
    fprintf(logfile, "%s%c", "AVG PARTIAL BB MAX DEPTH", SEP); count++;
    assert(count == countDisjInfoEntries);
  } // DISJ INFO
  
  // Certificate types will be used for STR and CERT info
  std::vector<std::string> CERTIFICATE_TYPES = { "ORIG ", "RCVMIP " };
  { // STR INFO
    int count = 0;
    for (const std::string& type : CERTIFICATE_TYPES) {
      for (const std::string& name : STR_INFO_CONTENTS) {
        fprintf(logfile, "%s%c", (type + name).c_str(), SEP); count++;
      }
    }
    // fprintf(logfile, "%s%c", "ORIG NUM STR AFFECTED CUTS", SEP); count++;
    // fprintf(logfile, "%s%c", "ORIG NUM COEFFS STR AVG", SEP); count++;
    // fprintf(logfile, "%s%c", "ORIG NUM COEFFS STR STDDEV", SEP); count++;
    // fprintf(logfile, "%s%c", "ORIG NUM COEFFS STR MIN", SEP); count++;
    // fprintf(logfile, "%s%c", "ORIG NUM COEFFS STR MAX", SEP); count++;
    // fprintf(logfile, "%s%c", "RCVMIP NUM STR AFFECTED CUTS", SEP); count++;
    // fprintf(logfile, "%s%c", "RCVMIP NUM COEFFS STR AVG", SEP); count++;
    // fprintf(logfile, "%s%c", "RCVMIP NUM COEFFS STR STDDEV", SEP); count++;
    // fprintf(logfile, "%s%c", "RCVMIP NUM COEFFS STR MIN", SEP); count++;
    // fprintf(logfile, "%s%c", "RCVMIP NUM COEFFS STR MAX", SEP); count++;
    assert(count == countStrInfoEntries);
  } // STR INFO

  { // CERTIFICATE INFO
    int count = 0;
    for (const std::string& type : CERTIFICATE_TYPES) {
      for (const std::string& name : CERT_INFO_CONTENTS) {
        fprintf(logfile, "%s%c", (type + name).c_str(), SEP); count++;
      }
    }
    assert(count == countCertInfoEntries);
  } // CERTIFICATE INFO

  { // CUT INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM ROUNDS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM CUTS", SEP); count++; // repeat, but it's ok
    fprintf(logfile, "%s%c", "NUM ONE SIDED CUTS", SEP); count++;
    fprintf(logfile, "%s%c", "NUM OPTIMALITY CUTS", SEP); count++;
    fprintf(logfile, "%s%c", "GOM MIN SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "GOM MAX SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "GOM AVG SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "MYCUTS MIN SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "MYCUTS MAX SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "MYCUTS AVG SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "UNSTR MYCUTS MIN SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "UNSTR MYCUTS MAX SUPPORT", SEP); count++;
    fprintf(logfile, "%s%c", "UNSTR MYCUTS AVG SUPPORT", SEP); count++;
    assert(count == countCutInfoEntries);
  } // CUT INFO
  { // OBJ INFO
    // For each objective: num obj, num fails, num active
    int count = 0;
    fprintf(logfile, "%s%c", "NUM OBJ", SEP); count++;
    /*for (int obj_ind = 0; obj_ind < static_cast<int>(CglAdvCut::ObjectiveType::NUM_OBJECTIVE_TYPES); obj_ind++) {
      fprintf(logfile, "NUM OBJ %s%c", CglAdvCut::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM CUTS %s%c", CglAdvCut::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM FAILS %s%c", CglAdvCut::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
      fprintf(logfile, "NUM ACTIVE %s%c", CglAdvCut::ObjectiveTypeName[obj_ind].c_str(), SEP); count++;
    }*/
    assert(count == countObjInfoEntries);
  } // OBJ INFO
  { // FAIL INFO
    int count = 0;
    fprintf(logfile, "%s%c", "NUM FAILS", SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES); fail_ind++) {
      fprintf(logfile, "%s%c", CglAdvCut::FailureTypeName[fail_ind].c_str(), SEP); count++;
    }
    assert(count == countFailInfoEntries);
  } // FAIL INFO
  { // FULL BB INFO
    int count = 0;
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("FIRST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("FIRST GUR+CUTS " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("BEST GUR+CUTS " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG GUR " + name).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", ("AVG GUR+CUTS " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL GUR " + name).c_str(), SEP); count++;
    }
    for (std::string name : BB_INFO_CONTENTS) {
      fprintf(logfile, "%s%c", ("ALL GUR+CUTS " + name).c_str(), SEP); count++;
    }
    assert(count == countFullBBInfoEntries);
  } // FULL BB INFO
  { // TIME INFO
    int count = 0;
    for (int t = 0; t < (int) time_name.size(); t++) {
      fprintf(logfile, "%s%c", time_name[t].c_str(), SEP); count++;
    }
    assert(count == countTimeInfoEntries);
  } // TIME INFO
  { // VERSION INFO
    fprintf(logfile, "%s%c", "code_version", SEP);
    fprintf(logfile, "%s%c", "vpc_version", SEP);
    fprintf(logfile, "%s%c", "cbc_version", SEP);
    fprintf(logfile, "%s%c", "clp_version", SEP);
    fprintf(logfile, "%s%c", "gurobi_version", SEP);
    fprintf(logfile, "%s%c", "cplex_version", SEP);
  } // VERSION INFO
  { // WRAP UP INFO
    fprintf(logfile, "%s%c", "ExitReason", SEP);
    fprintf(logfile, "%s%c", "end_time_string", SEP);
    fprintf(logfile, "%s%c", "time elapsed", SEP);
    fprintf(logfile, "%s%c", "instname", SEP);
  } // WRAP UP INFO

  fprintf(logfile, "\n");
  fflush(logfile);
} /* printHeader */

void printBoundAndGapInfo(const SummaryBoundInfo& boundInfo, FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  std::vector<double> rcvmip_bounds = { boundInfo.rcvmip_mycut_obj, boundInfo.rcvmip_gmic_mycut_obj, boundInfo.rcvmip_all_cuts_obj };

  { // BOUND INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(boundInfo.lp_obj, "%2.20f").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(boundInfo.best_disj_obj, "%2.20f").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(boundInfo.worst_disj_obj, "%2.20f").c_str(), SEP); count++;
    if (!isInfinity(std::abs(boundInfo.ip_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.ip_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_root_bounds_changed).c_str(), SEP); count++;
    if (!isInfinity(std::abs(boundInfo.root_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.root_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_gmic).c_str(), SEP); count++;
    if (boundInfo.num_gmic > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.gmic_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_lpc).c_str(), SEP); count++;
    if (boundInfo.num_lpc > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.lpc_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_mycut).c_str(), SEP); count++;
    if (boundInfo.num_mycut > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.mycut_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    if (!isInfinity(std::abs(boundInfo.gmic_mycut_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.gmic_mycut_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_str_affected_cuts).c_str(), SEP); count++;
    if (boundInfo.num_mycut > 0) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.unstr_mycut_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    if (!isInfinity(std::abs(boundInfo.unstr_gmic_mycut_obj))) {
      fprintf(logfile, "%s%c", stringValue(boundInfo.unstr_gmic_mycut_obj, "%2.20f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%c", SEP); count++;
    }
    fprintf(logfile, "%s%c", stringValue(boundInfo.num_rcvmip_str_affected_cuts).c_str(), SEP); count++;
    for (const auto& val : rcvmip_bounds) {
      if (!isInfinity(std::abs(val))) {
        fprintf(logfile, "%s%c", stringValue(val, "%2.20f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++;
      }
    }
    assert(count == countBoundInfoEntries);
  } // BOUND INFO

  { // GAP INFO
    int count = 0;
    if (!isInfinity(std::abs(boundInfo.ip_obj))) {
      if (!isInfinity(std::abs(boundInfo.gmic_obj))) {
        double val = 100. * (boundInfo.gmic_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // gmic
      }
      if (!isInfinity(std::abs(boundInfo.lpc_obj))) {
        double val = 100. * (boundInfo.lpc_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // lpc
      }
      if (!isInfinity(std::abs(boundInfo.mycut_obj))) {
        double val = 100. * (boundInfo.mycut_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // mycuts
      }
      if (!isInfinity(std::abs(boundInfo.gmic_mycut_obj))) {
        double val = 100. * (boundInfo.gmic_mycut_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // gmic_mycuts
      }
      if (!isInfinity(std::abs(boundInfo.unstr_mycut_obj))) {
        double val = 100. * (boundInfo.unstr_mycut_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // unstr_mycuts
      }
      if (!isInfinity(std::abs(boundInfo.unstr_gmic_mycut_obj))) {
        double val = 100. * (boundInfo.unstr_gmic_mycut_obj - boundInfo.lp_obj)
            / (boundInfo.ip_obj - boundInfo.lp_obj);
        fprintf(logfile, "%s%c", stringValue(val, "%2.6f").c_str(), SEP); count++;
      } else {
        fprintf(logfile, "%c", SEP); count++; // unstr_gmic_mycuts
      }
      for (const auto& bound : rcvmip_bounds) {
        if (!isInfinity(std::abs(bound))) {
          double val = 100. * (bound - boundInfo.lp_obj) / (boundInfo.ip_obj - boundInfo.lp_obj);
          fprintf(logfile, "%s%c", stringValue(val, "%2.20f").c_str(), SEP); count++;
        } else {
          fprintf(logfile, "%c", SEP); count++;
        }
      }
    } else {
      fprintf(logfile, "%c", SEP); count++; // gmic
      fprintf(logfile, "%c", SEP); count++; // lpc
      fprintf(logfile, "%c", SEP); count++; // mycuts
      fprintf(logfile, "%c", SEP); count++; // gmic_mycuts
      fprintf(logfile, "%c", SEP); count++; // unstr_mycuts
      fprintf(logfile, "%c", SEP); count++; // unstr_gmic_mycuts
      fprintf(logfile, "%c", SEP); count++; // rcvmip_mycuts
      fprintf(logfile, "%c", SEP); count++; // rcvmip_gmic_mycuts
      fprintf(logfile, "%c", SEP); count++; // rcvmip_all_cuts
    }
    assert(count == countGapInfoEntries);
  }
  fflush(logfile);
} /* printBoundAndGapInfo */

void printSummaryBBInfo(const std::vector<SummaryBBInfo>& info_vec, FILE* logfile,
    const bool print_blanks, const char SEP) {
  if (!logfile)
      return;

  int count = 0;
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%ld%c", info.first_bb_info.nodes, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%ld%c", info.first_bb_info.nodes, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%2.3f%c", info.first_bb_info.time, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  for (auto info : info_vec) {
    if (!print_blanks)
      fprintf(logfile, "%2.3f%c", info.first_bb_info.time, SEP);
    else
      fprintf(logfile, "%c", SEP);
    count++;
  }
  fflush(logfile);
  assert(count == countSummaryBBInfoEntries);
} /* printSummaryBBInfo */

void printFullBBInfo(const std::vector<SummaryBBInfo>& info_vec, FILE* logfile,
    const bool print_blanks, const char SEP) {
  if (!logfile)
    return;

//  const std::vector<bool> did_branch(info_vec.size());
//  for (unsigned i = 0; i < info_vec.size(); i++) {
//    did_branch[i] = info_vec[i].vec_bb_info.size() > 0;
//  }

  int count = 0;
  if (!print_blanks) {
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.first_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.best_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto& info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.avg_bb_info.obj, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.first_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.best_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%s%c", stringValue(info.avg_bb_info.bound, "%2.20f").c_str(), SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.iters, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.nodes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.first_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.best_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%ld%c", info.avg_bb_info.root_passes, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.first_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.best_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.avg_bb_info.first_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.first_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.best_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.20f%c", info.avg_bb_info.last_cut_pass, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.root_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.last_sol_time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.first_bb_info.time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.best_bb_info.time, SEP);
      count++;
    }
    for (auto info : info_vec) {
      fprintf(logfile, "%2.3f%c", info.avg_bb_info.time, SEP);
      count++;
    }

    // Finally, all
    for (auto info : info_vec) {
      std::vector<std::string> vec_str;
      createStringFromBBInfoVec(info.vec_bb_info, vec_str);
      for (unsigned i = 0; i < vec_str.size(); i++) {
        fprintf(logfile, "%s%c", vec_str[i].c_str(), SEP);
        count++;
      }
    }
  } else {
    // *4 for first, best, avg, all
    for (unsigned i = 0; i < BB_INFO_CONTENTS.size() * info_vec.size() * 4; i++) {
      fprintf(logfile, "%c", SEP); count++;
    }
  }
  fflush(logfile);
  assert(count == countFullBBInfoEntries);
} /* printFullBBInfo */

void printOrigProbInfo(const OsiSolverInterface* const solver, FILE* logfile,
    const char SEP) {
  if (!logfile)
    return;

  const int num_rows = solver->getNumRows();
  const int num_cols = solver->getNumCols();
  // Get row stats
  int num_eq_rows = 0, num_bound_rows = 0, num_assign_rows = 0;
  const CoinPackedMatrix* mat = solver->getMatrixByRow();
  for (int row = 0; row < num_rows; row++) {
    const double row_lb = solver->getRowLower()[row];
    const double row_ub = solver->getRowUpper()[row];
    if (isVal(row_lb, row_ub))
      num_eq_rows++;
    if (mat->getVectorSize(row) == 1) {
      if (isVal(row_lb, row_ub))
        num_assign_rows++;
      else
        num_bound_rows++;
    }
  }
  // Calculate fractionality
  int num_frac = 0;
  int num_fixed = 0, num_gen_int = 0, num_bin = 0, num_cont = 0;
  double min_frac = 1., max_frac = 0.;
  for (int col = 0; col < num_cols; col++) {
    const double col_lb = solver->getColLower()[col];
    const double col_ub = solver->getColUpper()[col];
    if (isVal(col_lb, col_ub))
      num_fixed++;
    if (!solver->isInteger(col)) {
      num_cont++;
      continue;
    }
    if (solver->isBinary(col))
      num_bin++;
    else
      num_gen_int++;
    const double val = solver->getColSolution()[col];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);
    const double frac = CoinMin(val - floorxk, ceilxk - val);
    if (!isVal(frac, 0., 1e-5)) {
      num_frac++;
      if (frac < min_frac)
        min_frac = frac;
      if (frac > max_frac)
        max_frac = frac;
    }
  }

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(num_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_cols).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_frac).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(min_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(max_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_eq_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_bound_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_assign_rows).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_fixed).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_gen_int).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_bin).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(num_cont).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue((double) mat->getNumElements() / (num_rows * num_cols)).c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countOrigProbEntries);
} /* printOrigProbInfo */

/**
 * Assumed that solver is already with cuts added
 */
void printPostCutProbInfo(const OsiSolverInterface* const solver,
    const SummaryCutInfo& cutInfoGMICs, const SummaryCutInfo& cutInfo,
    FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  const int num_rows = solver->getNumRows();
  const int num_cols = solver->getNumCols();

  // Calculate fractionality
  int num_frac = 0;
  double min_frac = 1., max_frac = 0.;
  for (int col = 0; col < num_cols; col++) {
    if (!solver->isInteger(col)) {
      continue;
    }
    const double val = solver->getColSolution()[col];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);
    const double frac = CoinMin(val - floorxk, ceilxk - val);
    if (!isVal(frac, 0., 1e-5)) {
      num_frac++;
      if (frac < min_frac)
        min_frac = frac;
      if (frac > max_frac)
        max_frac = frac;
    }
  }

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(num_frac).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(min_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(max_frac, "%.5f").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue((double) solver->getMatrixByCol()->getNumElements() / (num_rows * num_cols)).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_gmic).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfo.num_active_gmic).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_mycut).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfo.num_active_mycut).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.num_active_all).c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(cutInfo.num_active_all).c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countPostCutProbEntries);
} /* printPostCutProbInfo */

void printDisjInfo(const SummaryDisjunctionInfo& disjInfo, FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  int count = 0;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_terms, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.num_integer_sol, "%d").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.num_disj, "%d").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_density_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rows_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_cols_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_points_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_num_rays_prlp, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_explored_nodes, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_pruned_nodes, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_min_depth, "%g").c_str(), SEP); count++;
  fprintf(logfile, "%s%c", stringValue(disjInfo.avg_max_depth, "%g").c_str(), SEP); count++;
  fflush(logfile);
  assert(count == countDisjInfoEntries);
} /* printDisjInfo */

void printStrInfo(const SummaryStrengtheningInfo& orig_info, const SummaryStrengtheningInfo& rcvmip_info, FILE* const logfile, const char SEP) {
  if (!logfile)
    return;

  std::vector<SummaryStrengtheningInfo> info_vec = {orig_info, rcvmip_info};

  int count = 0;
  for (const auto& info : info_vec) {
    fprintf(logfile, "%s%c", stringValue(info.num_str_affected_cuts, "%d").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_coeffs_strengthened[(int) Stat::avg], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_coeffs_strengthened[(int) Stat::stddev], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_coeffs_strengthened[(int) Stat::min], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_coeffs_strengthened[(int) Stat::max], "%g").c_str(), SEP); count++;
  }
  fflush(logfile);
  assert(count == countStrInfoEntries);
} /* printStrInfo */

void printCertificateInfo(
    const SummaryCertificateInfo& orig_info,
    const SummaryCertificateInfo& rcvmip_info,
    const int RCVMIP_ITER_LIMIT,
    FILE* const logfile,
    const char SEP) {
  if (!logfile)
    return;
  
  std::vector<SummaryCertificateInfo> info_vec = {orig_info, rcvmip_info};

  int count = 0;
  for (const auto& info : info_vec) {
    fprintf(logfile, "%s%c", stringValue(info.num_unmatched_bounds, "%d").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.avg_num_cgs_facet, "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_irreg_less, "%d").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_reg, "%d").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_irreg_more, "%d").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(info.num_unconverged, "%d").c_str(), SEP); count++;
    
    std::vector<double> iter_stats = computeStats(info.num_iterations);
    fprintf(logfile, "%s%c", stringValue(iter_stats[(int) Stat::avg], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(iter_stats[(int) Stat::stddev], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(iter_stats[(int) Stat::min], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(iter_stats[(int) Stat::max], "%g").c_str(), SEP); count++;

    int num_limit = 0;
    for (const int curr_iter : info.num_iterations) {
      if (curr_iter >= RCVMIP_ITER_LIMIT)
        num_limit++;
    }
    fprintf(logfile, "%s%c", stringValue(num_limit, "%d").c_str(), SEP); count++;

    std::vector<double> time_stats = computeStats(info.rcvmip_time);
    fprintf(logfile, "%s%c", stringValue(time_stats[(int) Stat::avg], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(time_stats[(int) Stat::stddev], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(time_stats[(int) Stat::min], "%g").c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(time_stats[(int) Stat::max], "%g").c_str(), SEP); count++;
  }
  fflush(logfile);
  assert(count == countCertInfoEntries);
} /* printCertificateInfo */

void printCutInfo(const SummaryCutInfo& cutInfoGMICs,
    const SummaryCutInfo& cutInfo, const SummaryCutInfo& cutInfoUnstr,
    FILE* logfile, const char SEP) {
  if (!logfile)
    return;

  { // CUT INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_rounds).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_cuts).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.numCutsOfType[static_cast<int>(CglAdvCut::CutType::ONE_SIDED_CUT)]).c_str(), SEP); count++;
    fprintf(logfile, "%s%c", stringValue(cutInfo.numCutsOfType[static_cast<int>(CglAdvCut::CutType::OPTIMALITY_CUT)]).c_str(), SEP); count++;
    if (cutInfoGMICs.num_cuts > 0) {
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.min_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.max_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoGMICs.avg_support, "%.3f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
    }
    if (cutInfo.num_cuts > 0) {
      fprintf(logfile, "%s%c", stringValue(cutInfo.min_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.max_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfo.avg_support, "%.3f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
    }
    if (cutInfoUnstr.num_cuts > 0) {
      fprintf(logfile, "%s%c", stringValue(cutInfoUnstr.min_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoUnstr.max_support).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(cutInfoUnstr.avg_support, "%.3f").c_str(), SEP); count++;
    } else {
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
      fprintf(logfile, "%s%c", stringValue(0).c_str(), SEP); count++;
    }
    assert(count == countCutInfoEntries);
  } // CUT INFO
  { // OBJ INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_obj_tried).c_str(), SEP); count++;
    assert(count == countObjInfoEntries);
  } // OBJ INFO
  { // FAIL INFO
    int count = 0;
    fprintf(logfile, "%s%c", stringValue(cutInfo.num_failures).c_str(), SEP); count++;
    for (int fail_ind = 0; fail_ind < static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES); fail_ind++) {
      fprintf(logfile, "%s%c", stringValue(cutInfo.numFails[fail_ind]).c_str(), SEP); count++;
    }
    assert(count == countFailInfoEntries);
  } // FAIL INFO
  fflush(logfile);
} /* printCutInfo */

/// @details Gets cut support size and updates min/max component of \p cutInfo
int checkCutDensity(
    /// [in,out] Where to save min and max support
    SummaryCutInfo& cutInfo,
    /// [in] Row that we want to check
    const OsiRowCut* const cut,
    /// [in] What counts as a zero coefficient
    const double EPS) {
  int num_elem = cut->row().getNumElements();
  const double* el = cut->row().getElements();
  for (int i = 0; i < cut->row().getNumElements(); i++) {
    if (isZero(el[i], EPS)) {
      num_elem--;
    }
  }
  if (num_elem < cutInfo.min_support)
    cutInfo.min_support = num_elem;
  if (num_elem > cutInfo.max_support)
    cutInfo.max_support = num_elem;
  return num_elem;
} // checkCutDensity

bool checkCutActivity(
  const OsiSolverInterface* const solver,
  const OsiRowCut* const cut) {
  if (solver && solver->isProvenOptimal()) {
    const double activity = dotProduct(cut->row(), solver->getColSolution());
    return isVal(activity, cut->rhs());
  } else {
    return false;
  }
} /* checkCutActivity */

/// @details Returns how many cuts are violated by a given integer-feasible solution
int checkCutsAgainstFeasibleSolution(
    const OsiCuts& currCuts, ///< [in] Cuts to check
    const std::vector<double> ip_solution ///< [in] Feasible solution
) {
  int num_violated = 0;
  for (int cut_ind = 0; cut_ind < currCuts.sizeCuts(); cut_ind++) {
    OsiRowCut currCut = currCuts.rowCut(cut_ind);
    const double rhs = currCut.rhs();
    const int num_el = currCut.row().getNumElements();
    const int* ind = currCut.row().getIndices();
    const double* el = currCut.row().getElements();
    const double activity = dotProduct(num_el, ind, el, ip_solution.data());

    if (lessThanVal(activity, rhs)) {
      num_violated++;
      warning_msg(warnstring, "Unstrengthened cut %d removes optimal solution. Activity: %.10f. Rhs: %.10f.\n", cut_ind, activity, rhs);
    }
  } // loop over cuts
  return num_violated;
} /* checkCutsAgainstFeasibleSolution */

/**
 * @details The cut properties we want to look at are:
 * 1. Gap closed
 * 2. Activity (after adding cuts)
 * 3. Density
 */
void analyzeStrength(
    const StrengtheningParameters::Parameters& params, 
    const OsiSolverInterface* const solver_gmic,
    const OsiSolverInterface* const solver_mycut,
    const OsiSolverInterface* const solver_all,
    SummaryCutInfo& cutInfoGMICs, SummaryCutInfo& cutInfo,
    const OsiCuts* const gmics, const OsiCuts* const mycuts,
    const SummaryBoundInfo& boundInfo, std::string& output) {
  cutInfoGMICs.num_active_gmic = 0;
  cutInfoGMICs.num_active_mycut = 0;
  cutInfoGMICs.num_active_all = 0;
  cutInfo.num_active_gmic = 0;
  cutInfo.num_active_mycut = 0;
  cutInfo.num_active_all = 0;
  cutInfo.numActiveFromHeur.resize(static_cast<int>(CglAdvCut::ObjectiveType::NUM_OBJECTIVE_TYPES), 0);
  if (mycuts) {
    const int num_mycuts = mycuts->sizeCuts();
    int total_support = 0;
    for (int cut_ind = 0; cut_ind < num_mycuts; cut_ind++) {
      const OsiRowCut* const cut = mycuts->rowCutPtr(cut_ind);
      if (checkCutActivity(solver_gmic, cut)) {
        cutInfo.num_active_gmic++;
      }
      if (checkCutActivity(solver_mycut, cut)) {
        cutInfo.num_active_mycut++;
        if (cut_ind < static_cast<int>(cutInfo.objType.size())) {
          cutInfo.numActiveFromHeur[static_cast<int>(cutInfo.objType[cut_ind])]++;
        }
      }
      if (checkCutActivity(solver_all, cut)) {
        cutInfo.num_active_all++;
      }
      total_support += checkCutDensity(cutInfo, cut, params.get(EPS) / 2.);
    }
    cutInfo.avg_support = (double) total_support / num_mycuts;
  }
  if (gmics) {
    const int num_gmics = gmics->sizeCuts();
    cutInfoGMICs.num_cuts = num_gmics;
    int total_support = 0;
    for (int cut_ind = 0; cut_ind < num_gmics; cut_ind++) {
      const OsiRowCut* const cut = gmics->rowCutPtr(cut_ind);
      if (checkCutActivity(solver_gmic, cut)) {
        cutInfoGMICs.num_active_gmic++;
      }
      if (checkCutActivity(solver_mycut, cut)) {
        cutInfoGMICs.num_active_mycut++;
      }
      if (checkCutActivity(solver_all, cut)) {
        cutInfoGMICs.num_active_all++;
      }
      total_support += checkCutDensity(cutInfoGMICs, cut, params.get(EPS) / 2.);
    }
    cutInfoGMICs.avg_support = (double) total_support / num_gmics;
  }

  // Print results from adding cuts
  int NAME_WIDTH = 25;
  int NUM_DIGITS_BEFORE_DEC = 7;
  int NUM_DIGITS_AFTER_DEC = 7;
  const double INF = std::numeric_limits<double>::max();
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
      "\n## Results from adding cuts ##\n");
  output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
      NAME_WIDTH, NAME_WIDTH, "LP: ",
      stringValue(boundInfo.lp_obj, "% -*.*g",
        INF,
        NUM_DIGITS_BEFORE_DEC,
        NUM_DIGITS_AFTER_DEC).c_str());
  output += tmpstring;
  if (!isInfinity(std::abs(boundInfo.root_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d bounds changed)\n", NAME_WIDTH, NAME_WIDTH, "Root: ",
        stringValue(boundInfo.root_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_root_bounds_changed);
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.gmic_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "GMICs: ",
        stringValue(boundInfo.gmic_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_gmic);
    output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active GMICs", cutInfoGMICs.num_active_gmic);
    output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active MYCUTs", cutInfo.num_active_gmic);
    output += tmpstring;
    output += ")\n";
  }
  if (boundInfo.num_str_affected_cuts > 0 && !isInfinity(std::abs(boundInfo.unstr_mycut_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "unstr MYCUTs: ",
        stringValue(boundInfo.unstr_mycut_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_mycut);
    output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        " -> %d strengthened", boundInfo.num_str_affected_cuts);
    output += tmpstring;
    output += ")\n";
  }
  if (!isInfinity(std::abs(boundInfo.mycut_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "MYCUTs: ",
        stringValue(boundInfo.mycut_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_mycut);
    output += tmpstring;
    if (gmics && gmics->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active GMICs", cutInfoGMICs.num_active_mycut);
      output += tmpstring;
    }
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        ", %d active MYCUTs", cutInfo.num_active_mycut);
    output += tmpstring;
    output += ")\n";
  }
  if (boundInfo.num_str_affected_cuts > 0 && !isInfinity(std::abs(boundInfo.unstr_mycut_obj))) {
    if (boundInfo.num_gmic + boundInfo.num_lpc > 0) { // Even if there are no MYCUTs, but not if there are *only* MYCUTs
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          "%-*.*s%s (%d cuts)\n", NAME_WIDTH, NAME_WIDTH, "unstr All: ",
          stringValue(boundInfo.unstr_all_cuts_obj, "% -*.*g",
            INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(),
          boundInfo.num_gmic + boundInfo.num_lpc + boundInfo.num_mycut);
      output += tmpstring;
    }
  }
  if (boundInfo.num_gmic + boundInfo.num_lpc > 0) { // Even if there are no MYCUTs, but not if there are *only* MYCUTs
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts", NAME_WIDTH, NAME_WIDTH, "All: ",
        stringValue(boundInfo.all_cuts_obj, "% -*.*g", INF, NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str(), 
        boundInfo.num_gmic + boundInfo.num_lpc + boundInfo.num_mycut);
    output += tmpstring;
    if (gmics && gmics->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active GMICs", cutInfoGMICs.num_active_all);
      output += tmpstring;
    }
    if (mycuts && mycuts->sizeCuts() > 0) {
      snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
          ", %d active MYCUTs", cutInfo.num_active_all);
      output += tmpstring;
    }
    output += ")\n";
  }

  // Finally, print effect when RCVMIP-strengthened cuts are used, if applicable
  if (boundInfo.num_rcvmip_str_affected_cuts > 0) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char),
        "%-*.*s%s (%d cuts affected by RCVMIP", NAME_WIDTH, NAME_WIDTH, "RCVMIP: ",
        stringValue(boundInfo.rcvmip_all_cuts_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str(),
        boundInfo.num_rcvmip_str_affected_cuts);
    output += tmpstring;
    output += ")\n";
  }

  if (!isInfinity(std::abs(boundInfo.best_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive lb: ",
        stringValue(boundInfo.best_disj_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.worst_disj_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "Disjunctive ub: ",
        stringValue(boundInfo.worst_disj_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
  if (!isInfinity(std::abs(boundInfo.ip_obj))) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*.*s%s\n",
        NAME_WIDTH, NAME_WIDTH, "IP: ",
        stringValue(boundInfo.ip_obj, "% -*.*g",
          INF,
          NUM_DIGITS_BEFORE_DEC,
          NUM_DIGITS_AFTER_DEC).c_str());
    output += tmpstring;
  }
} /* analyzeStrength */

void analyzeBB(const StrengtheningParameters::Parameters& params, SummaryBBInfo& info_nocuts,
    SummaryBBInfo& info_mycuts, SummaryBBInfo& info_allcuts, std::string& output) {
  if (params.get(BB_RUNS) == 0) {
    return;
  }
  // B&B mode: ones bit = no_cuts, tens bit = w/mycuts, hundreds bit = w/gmics
  const int mode_param = params.get(intParam::BB_MODE);
  const int mode_ones = mode_param % 10;
  const int mode_tens = (mode_param % 100 - (mode_param % 10)) / 10;
  const int mode_hundreds = (mode_param % 1000 - (mode_param % 100)) / 100;
  const bool branch_with_no_cuts = (mode_ones > 0);
  const bool branch_with_mycuts = (mode_tens > 0) && (info_mycuts.num_cuts > 0);
  const bool branch_with_gmics = (mode_hundreds > 0) && (info_allcuts.num_cuts > 0);
  if (branch_with_no_cuts + branch_with_mycuts + branch_with_gmics == 0) {
    return;
  }

  // Save results to string and also print to the logfile
  int NAME_WIDTH = 10; //25
  int NUM_DIGITS_BEFORE_DEC = 15; //10
  int NUM_DIGITS_AFTER_DEC = 2; //2
  char tmpstring[300];

  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n## Branch-and-bound results ##\n"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, ""); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Obj"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Bound"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Iters"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Nodes"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Root passes"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "First cut pass"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Last cut pass"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Root time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Last sol time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NUM_DIGITS_BEFORE_DEC, "Time"); output += tmpstring;
  snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  if (branch_with_no_cuts) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_nocuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_nocuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // no cuts
  if (branch_with_mycuts) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur+MyCuts"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_mycuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_mycuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // mycuts
  if (branch_with_gmics) {
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*s", NAME_WIDTH, "Gur+MyCuts+G"); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.obj, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.bound, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.iters); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.nodes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%-*ld", NUM_DIGITS_BEFORE_DEC, info_allcuts.avg_bb_info.root_passes); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.first_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_cut_pass, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.root_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.last_sol_time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "%s", stringValue(info_allcuts.avg_bb_info.time, "%-*.*f", NUM_DIGITS_BEFORE_DEC, NUM_DIGITS_AFTER_DEC).c_str()); output += tmpstring;
    snprintf(tmpstring, sizeof(tmpstring) / sizeof(char), "\n"); output += tmpstring;
  } // gmics
} /* analyzeBB */

double getNumGomoryRounds(const StrengtheningParameters::Parameters& params,
    const OsiSolverInterface* const origSolver,
    const OsiSolverInterface* const postCutSolver) {
  // Get number rounds of SICs needed to meet bound from mycuts+SICs
#ifdef TRACE
  printf("\nGetting number rounds of Gomory cuts req'd to get bound.\n");
#endif
  const int num_cuts = postCutSolver->getNumRows() - origSolver->getNumRows();
  const double post_cut_opt = postCutSolver->getObjValue();
  const int min_sic_rounds = (params.get(STRENGTHEN) == 2) ? 2 : 0;
  int max_rounds = 1000;

  int total_num_sics = 0;
  int num_sic_rounds = 0;
  double curr_sic_opt = 0.;
  std::vector<int> numCutsByRoundSIC;
  std::vector<double> boundByRoundSIC;
  OsiSolverInterface* copySolver = origSolver->clone();
  if (!copySolver->isProvenOptimal()) {
    copySolver->initialSolve();
    checkSolverOptimality(copySolver, false);
  }
  while (num_sic_rounds < min_sic_rounds
      || (lessThanVal(curr_sic_opt, post_cut_opt) && total_num_sics < num_cuts)) {
    OsiCuts GMICs;
    CglGMI gen;
    gen.generateCuts(*copySolver, GMICs);
    const int curr_num_cuts = GMICs.sizeCuts();
    if (curr_num_cuts == 0)
      break;

    num_sic_rounds++;
    total_num_sics += curr_num_cuts;
    numCutsByRoundSIC.push_back(curr_num_cuts);
    curr_sic_opt = applyCutsCustom(copySolver, GMICs);
    boundByRoundSIC.push_back(curr_sic_opt);

    // Other stopping conditions:
    // Bound does not improve at all after one round
    if (num_sic_rounds >= 2
        && !greaterThanVal(curr_sic_opt, boundByRoundSIC[num_sic_rounds - 2])) {
      break;
    }
    // Bound does not significantly improve after five rounds
    if (num_sic_rounds > 4) {
      const double delta = curr_sic_opt - boundByRoundSIC[num_sic_rounds - 4];
      if (!greaterThanVal(delta, 1e-3)) {
        break;
      }
    }
  } // do rounds of Gomory cuts
  if (max_rounds < num_sic_rounds) {
    max_rounds = boundByRoundSIC.size();
  }
  const double final_sic_bound = copySolver->getObjValue();
  return final_sic_bound;
} /* getNumGomoryRounds */

/**
 * @brief Use this to add to cutInfo
 *
 * @details Use this to add to cutInfo (but within one round,
 * because the cutType and objType vectors are cleared in gen in each round
 * (so tracking that based on isSetupForRepeatedUse does not work,
 * and the old cutType and objType stored in cutInfo would be overwritten)
 */
void updateCutInfo(SummaryCutInfo& cutInfo, const CglAdvCut* const gen) {
  cutInfo.num_cuts += gen->num_cuts;
  cutInfo.num_obj_tried += gen->num_obj_tried;
  cutInfo.num_failures += gen->num_failures;

  // For cutType and objType, what we do depends on whether the generator is setup for repeated use or not
  if (gen->isSetupForRepeatedUse) {
    cutInfo.cutType = gen->cutType;
    cutInfo.objType = gen->objType;
  } else {
    cutInfo.cutType.insert(cutInfo.cutType.end(), gen->cutType.begin(), gen->cutType.end());
    cutInfo.objType.insert(cutInfo.objType.end(), gen->objType.begin(), gen->objType.end());
  }

  if (cutInfo.numCutsOfType.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglAdvCut::CutType::NUM_CUT_TYPES); i++) {
      cutInfo.numCutsOfType[i] += gen->numCutsOfType[i];
    }
  } else {
    cutInfo.numCutsOfType = gen->numCutsOfType;
  }

  if (cutInfo.numCutsFromHeur.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglAdvCut::ObjectiveType::NUM_OBJECTIVE_TYPES); i++) {
      cutInfo.numCutsFromHeur[i] += gen->numCutsFromHeur[i];
      cutInfo.numObjFromHeur[i] += gen->numObjFromHeur[i];
      cutInfo.numFailsFromHeur[i] += gen->numFailsFromHeur[i];
    }
  } else {
    cutInfo.numCutsFromHeur = gen->numCutsFromHeur;
    cutInfo.numObjFromHeur = gen->numObjFromHeur;
    cutInfo.numFailsFromHeur = gen->numFailsFromHeur;
  }

  if (cutInfo.numFails.size() > 0) {
    for (int i = 0; i < static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES); i++) {
      cutInfo.numFails[i] += gen->numFails[i];
    }
  } else {
    cutInfo.numFails = gen->numFails;
  }
} /* updateCutInfo (within one round) */

/**
 * @brief Use this to merge cut info from multiple rounds
 */
void setCutInfo(SummaryCutInfo& cutInfo, const int num_rounds,
    const SummaryCutInfo* const oldCutInfos) {
  const int numCutTypes = static_cast<int>(CglAdvCut::CutType::NUM_CUT_TYPES);
  const int numObjTypes = static_cast<int>(CglAdvCut::ObjectiveType::NUM_OBJECTIVE_TYPES);
  const int numFailTypes = static_cast<int>(CglAdvCut::FailureType::NUM_FAILURE_TYPES);

  cutInfo.num_cuts = 0;
  cutInfo.num_active_gmic = 0;
  cutInfo.num_active_mycut = 0;
  cutInfo.num_active_all = 0;
  cutInfo.num_obj_tried = 0;
  cutInfo.num_failures = 0;
  cutInfo.num_rounds = num_rounds;
  cutInfo.cutType.resize(0);
  cutInfo.objType.resize(0);
  cutInfo.numCutsOfType.clear();
  cutInfo.numCutsOfType.resize(numCutTypes, 0);
  cutInfo.numCutsFromHeur.clear();
  cutInfo.numCutsFromHeur.resize(numObjTypes, 0);
  cutInfo.numObjFromHeur.clear();
  cutInfo.numObjFromHeur.resize(numObjTypes, 0);
  cutInfo.numFailsFromHeur.clear();
  cutInfo.numFailsFromHeur.resize(numObjTypes, 0);
  cutInfo.numActiveFromHeur.clear();
  cutInfo.numActiveFromHeur.resize(numObjTypes, 0);
  cutInfo.numFails.clear();
  cutInfo.numFails.resize(numFailTypes, 0);

  for (int round = 0; round < num_rounds; round++) {
    cutInfo.num_cuts += oldCutInfos[round].num_cuts;
    cutInfo.num_active_gmic += oldCutInfos[round].num_active_gmic;
    cutInfo.num_active_mycut += oldCutInfos[round].num_active_mycut;
    cutInfo.num_active_all += oldCutInfos[round].num_active_all;
    cutInfo.num_obj_tried += oldCutInfos[round].num_obj_tried;
    cutInfo.num_failures += oldCutInfos[round].num_failures;

    for (int i = 0; i < numCutTypes; i++) {
      cutInfo.numCutsOfType[i] += oldCutInfos[round].numCutsOfType[i];
    }
    for (int i = 0; i < numObjTypes; i++) {
      cutInfo.numCutsFromHeur[i] += oldCutInfos[round].numCutsFromHeur[i];
      cutInfo.numObjFromHeur[i] += oldCutInfos[round].numObjFromHeur[i];
      cutInfo.numFailsFromHeur[i] += oldCutInfos[round].numFailsFromHeur[i];
    }
    if (oldCutInfos[round].numActiveFromHeur.size() > 0) {
      for (int i = 0; i < numObjTypes; i++) {
        cutInfo.numActiveFromHeur[i] += oldCutInfos[round].numActiveFromHeur[i];
      }
    }
    for (int i = 0; i < numFailTypes; i++) {
      cutInfo.numFails[i] += oldCutInfos[round].numFails[i];
    }
  }

  cutInfo.cutType.resize(cutInfo.num_cuts);
  cutInfo.objType.resize(cutInfo.num_cuts);
  int cut_ind = 0;
  for (int round = 0; round < num_rounds; round++) {
    for (int i = 0; i < oldCutInfos[round].num_cuts; i++) {
      if (cut_ind < static_cast<int>(oldCutInfos[round].cutType.size())) {
        cutInfo.cutType[cut_ind] = oldCutInfos[round].cutType[i];
      } else {
        cutInfo.cutType[cut_ind] = CglAdvCut::CutType::OTHER_CUT;
      }
      if (cut_ind < static_cast<int>(oldCutInfos[round].objType.size())) {
        cutInfo.objType[cut_ind] = oldCutInfos[round].objType[i];
      } else {
        cutInfo.objType[cut_ind] = CglAdvCut::ObjectiveType::OTHER;
      }
      cut_ind++;
    }
  }
} /* setCutInfo (merge from multiple rounds) */

/// @details For a collection of cuts, count number of nonzero multipliers,
/// and also find number of distinct facets of the convex cut-generating set
/// defined by these multipliers (applied to \p disj)
///
/// TODO only works with bounds defining disjunctive terms for now
void setCertificateInfo(
    SummaryCertificateInfo& info,
    const Disjunction* const disj,
    /// [in] Index is [term][multiplier] where the first m multipliers are on original rows, then there are m0 on the disjunctive term ineqs, and finally n on columns
    const std::vector<CutCertificate>& v_vec,
    const int num_rows,
    const int num_cols,
    /// [in] Indices of strengthened cuts
    const std::vector<int>& str_cut_ind,
    const double EPS) {
  if (!disj) { return; }

  assert( disj->common_ineqs.size() == 0 ); // TODO handle this case
  const int num_common_rows = disj->common_changed_bound.size() + disj->common_ineqs.size();
  const int num_str_cuts = (int) str_cut_ind.size();
  
  for (const int cut_ind : str_cut_ind) {
    const CutCertificate& v = v_vec[cut_ind];

    // K counts the number of nonzero multipliers for each constraint
    // If the constraints form a basis, then the cut is "regular",
    // assuming the corresponding point is not feasible for the disjunction
    std::vector<int> K(num_rows + num_common_rows + num_cols, 0);
    // int num_nonzero_coeff = 0;

    // Compute the convex set S using the multipliers on the disjunctive term ineqs
    // These will all be stored as inequalities in >= form
    std::vector<CoinPackedVector> facetLHS;
    std::vector<double> facetRHS;
    facetLHS.reserve(disj->num_terms);
    facetRHS.reserve(disj->num_terms);
    for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
      const DisjunctiveTerm& term = disj->terms[term_ind];
      const int num_disj_ineqs = (int) term.changed_var.size();

      // For each original inequality, check whether it has a nonzero multiplier
      for (int row = 0; row < num_rows + num_common_rows; row++) {
        const double ukt = v[term_ind][row];
        if (!isZero(ukt, EPS)) {
          // if (K[row] == 0) num_nonzero_coeff++;
          K[row]++;
        }
      } // loop over original rows
      
      // For each disj term inequality, do the aggregation using u^t_0
      // TODO generalize for general inequalities (need to account for coefficients)
      std::vector<int> indices;
      indices.reserve(num_disj_ineqs);
      std::vector<double> elements;
      elements.reserve(num_disj_ineqs);
      double rhs = 0.;
      for (int bound_ind = 0; bound_ind < num_disj_ineqs; bound_ind++) {
        const double uk0 = v[term_ind][num_rows + num_common_rows + bound_ind];
        if (isZero(uk0, EPS)) continue;

        const int var = term.changed_var[bound_ind];
        const int bd = term.changed_bound[bound_ind]; // <= 0: lower bound, 1: upper bound
        const double val = term.changed_value[bound_ind];
        const double mult = (bd <= 0) ? 1. : -1.; // if bd == 1, then -x_k >= -val is the term
        const double coeff = uk0 * mult;
        const double curr_rhs = uk0 * val * mult;
        
        // Check if this variable appeared in a previous inequality for this term
        // This could happen if we branch (x_k <= 1) v (x_k >= 2) and later (x_k <= 0) v (x_k >= 1)
        // But in that case, only one of those should have a nonzero multiplier
        // On the other hand, for disjunctions not restricted to variable branching,
        // the situation can arise that a nonzero multiplier exists on two constraints containing the same var
        int prev_ind = -1;
        for (int i = 0; i < (int) indices.size(); i++) {
          if (indices[i] == var) {
            prev_ind = i;
            break;
          }
        }

        rhs += curr_rhs;
        if (prev_ind == -1) {
          indices.push_back(var);
          elements.push_back(coeff);
        } else {
          elements[prev_ind] += coeff;
        }
      } // loop over disj term ineqs
      if ((int) indices.size() > 0) {
        facetLHS.push_back(CoinPackedVector((int) indices.size(), indices.data(), elements.data()));
        facetRHS.push_back(rhs);
      }
    } // loop over terms

  #ifdef TRACE
    // Print facets generated
    std::string cgsName = "";
    for (int facet_ind = 0; facet_ind < (int) facetLHS.size(); facet_ind++) {
      CoinPackedVector& vec = facetLHS[facet_ind];
      const double rhs = facetRHS[facet_ind];
      const int num_elem = vec.getNumElements();
      const int* indices = vec.getIndices();
      const double* elements = vec.getElements();
      Disjunction::setCgsName(cgsName, num_elem, indices, elements, rhs, false);
    }
    printf("setCertificateInfo: After aggregating disjunctive terms, convex cut-generating set has following name:\n");
    printf("%s\n", cgsName.c_str());
  #endif

    // Loop over the new facets and check which are distinct
    int num_facets = 0;
    std::vector<int> sameAsFacet(facetLHS.size(), -1);
    for (int facet_ind = 0; facet_ind < (int) facetLHS.size(); facet_ind++) {
      CoinPackedVector& vec1 = facetLHS[facet_ind];
      const double rhs1 = facetRHS[facet_ind];
      int f = 0;
      for (f = 0; f < facet_ind; f++) {
        if (sameAsFacet[f] != -1) continue;
        CoinPackedVector& vec2 = facetLHS[f];
        const double rhs2 = facetRHS[f];
        const int howDifferent = isRowDifferent(&vec1, rhs1, &vec2, rhs2, EPS);
        if (howDifferent == 0) {
          sameAsFacet[facet_ind] = f;
          break;
        }
      } // loop over previous facets
      if (sameAsFacet[facet_ind] == -1) {
        num_facets++;
      }
    } // loop over cgs facets
    info.avg_num_cgs_facet += (double) num_facets / num_str_cuts;

    // Update number of unmatched bounds (times both a lower and upper bound on a variable have
    // nonzero multipliers, which can be detected by the sign on the corresponding bound in v)
    for (int var = 0; var < num_cols; var++) {
      int lt_zero_ind = -1, gt_zero_ind = -1;
      for (int term_ind = 0; term_ind < disj->num_terms; term_ind++) {
        const int num_disj_ineqs = (int) disj->terms[term_ind].changed_var.size();
        const double ukt = v[term_ind][num_rows + num_common_rows + num_disj_ineqs + var];
        if (lessThanVal(ukt, 0, EPS)) {
          if (lt_zero_ind == -1) lt_zero_ind = term_ind;
        } else if (greaterThanVal(ukt, 0, EPS)) {
          if (gt_zero_ind == -1) gt_zero_ind = term_ind;
        }
        if (lt_zero_ind != -1 && gt_zero_ind != -1) break;

        if (!isZero(ukt, EPS)) {
          // if (K[num_rows + num_common_rows + var] == 0) num_nonzero_coeff++;
          K[num_rows + num_common_rows + var]++;
        }
      } // loop over terms
      info.num_unmatched_bounds += (lt_zero_ind != -1 && gt_zero_ind != -1);
    } // loop over vars
  } // loop over each strengthened cut
} /* setCertificateInfo */
