/**
 * @file lift.cpp
 * @author A. M. Kazachkov
 * @date 2019-11-18
 */
#include "lift.hpp"

/*****************************************/
/** Lifting in the (u,v)-variable space **/
/*****************************************/

/**
 * Set up integrated lifting LP in the u, v space.
 */
void setupLandPUVSpace(SolverInterface* Dk,
    SolverInterface* subprob, OsiCut &cutinfo, 
    const int currvar) {
#ifdef TRACE
  printf("\n## Setting up Dk LP. ##\n");
#endif

  genLandPUVSpace(Dk, cutinfo, subprob);
#ifdef TRACE
  printf("D%d: # rows: %d. # cols: %d.\n", currvar, Dk.getNumRows(),
      Dk.getNumCols());
  printf("\n");
#endif
} /* setupLandPUVSPace */

/**
 * Update integrated lifting LP in the u, v space
 * with info for the new split.
 */
void updateLandPUVSpace(SolverInterface* Dk, OsiCut &cutinfo, const int currvar,
    const int oldSplitIndex, const double floorxk, const double ceilxk, 
    const int numSubRows, const int numSubCols, const Parameters& params) {
  // First the right sides for the baralpha rows
  for (int col = 0; col < (int) cutinfo.num_coeff; col++) {
    double alphaj = cutinfo[col];
    Dk.setRowBounds(col, alphaj, alphaj); // x_j^0
    Dk.setRowBounds(numSubCols + 1 + col, alphaj, alphaj); // x_j^1
  }

  // Now change the column u0 and v0
  Dk.modifyCoefficient(numSubCols, numSubRows, floorxk);
  Dk.modifyCoefficient(2 * (numSubCols + 1) - 1, 2 * (numSubRows + 1) - 1,
      -1 * ceilxk);

  // We also switch so that column u_0 has the 1 in the right row
  // and v_0 has a -1 in the correct row
  Dk.modifyCoefficient(oldSplitIndex, numSubRows, 0.0);
  Dk.modifyCoefficient(numSubCols + 1 + oldSplitIndex, 2 * numSubRows + 1,
      0.0);
  Dk.modifyCoefficient(currvar, numSubRows, -1.0);
  Dk.modifyCoefficient(numSubCols + 1 + currvar, 2 * numSubRows + 1, 1.0);
} /* updateLandPUVSpace */

/**
 * Generates lifting in the L&P style.
 */
void liftingLandPUVSpace(OsiCuts &liftedCut, OsiCuts &orderedLiftedCut,
    OsiCut &cutinfo, SolverInterface* subprob,
    SolverInterface* solver, 
    /*
    const std::vector<int>& deletedCols,
    const std::vector<int>& oldRowIndices, 
    const std::vector<int>& oldColIndices,
    */
    const int max_pivots, int &num_pivots, 
    const std::string &cutname, 
    SolverInterface* Dk_init,
    const int subspace_option, 
    FILE *inst_info_out) {

  SolverInterface Dk(Dk_init);

  // Save Dk for posterity
  std::string dir = "", instname = "", in_file_ext = "", filename = "";
  parseFilename(dir, instname, in_file_ext, params);
  filename = dir + "/" + instname;
  char Dk_name[256];
  snprintf(Dk_name, sizeof(Dk_name) / sizeof(char), "%s-D%s",
      filename.c_str(), cutname.c_str());
#ifdef WRITE_CGLP
  if (WRITE_LP_FILE) {
    //Dk.writeLp(Dk_name);
    Dk.writeMps(Dk_name, "mps");
  }
#endif

  // Perform lifting using first solution
#ifdef TRACE
  printf("\n## Obtaining first lifting. ##\n");
#endif
  Dk.enableFactorization();
  Dk.initialSolve();

  if (Dk.isProvenOptimal()) {
#ifdef TRACE
    printf("Dk opt value: %f\n", Dk.getObjValue());
#endif
  } else {
    error_msg(errorstring, "%s not proven optimal!\n", Dk_name);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }
//  printf("\n## Saving alternate copy. ##\n");
//  SolverInterface DkCopy(Dk.getModelPtr());
//  DkCopy.enableFactorization();
//  DkCopy.initialSolve();

#ifdef TRACE
  printf("Saving dual LP solution information.\n");
#endif
  SolutionInfo solnInfoDual(Dk, false);
  if (SHOULD_WRITE_BRIEF_SOLN)
    solnInfoDual.printBriefSolutionInfo(Dk, Dk_name);

  if (SHOULD_WRITE_SOLN) {
    char soln_name[300];
    snprintf(soln_name, sizeof(soln_name) / sizeof(char), "%sInit",
        Dk_name);
    writeSoln(&Dk, soln_name);
    char OptTabName[300];
    snprintf(OptTabName, sizeof(OptTabName) / sizeof(char),
        "%s-Opt_Tableau.csv", soln_name);
    FILE* OptTabOut = fopen(OptTabName, "w");
    if (OptTabOut != NULL)
      printSimplexTableauWithNames(OptTabOut, Dk);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write optimal simplex tableau.\n",
          OptTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }

    char NBTabName[300];
    snprintf(NBTabName, sizeof(NBTabName) / sizeof(char),
        "%s-NB_Opt_Tableau.csv", soln_name);
    FILE* NBTabOut = fopen(NBTabName, "w");
    if (NBTabOut != NULL)
      printNBSimplexTableauWithNames(NBTabOut, Dk);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write NB optimal simplex tableau.\n",
          NBTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }
  }

  liftUsingCoeffsUVSpace(liftedCut, orderedLiftedCut, cutinfo, subprob,
      solver, Dk, Dk, deletedCols, oldRowIndices, oldColIndices,
      subspace_option, inst_info_out);

  /***********************************************************************************
   * Lift using any neighboring alternative optimal lifting coefficients
   ***********************************************************************************/
  if (max_pivots == 0) // We do not do the following section if we do not want to do any pivots
    return;

  vector<int> colIn, colInNBIndex, inStatus, colOut, outStatus;
  vector<double> minRatio;

#ifdef TRACE
  printf("\n## Getting pivots. ##\n");
#endif
  getPivots(Dk, solnInfoDual, colIn, colInNBIndex, inStatus, colOut,
      outStatus, minRatio, num_pivots, inst_info_out);

//  SolverInterface tmpdualsolver(dualsolver); // In case of singular factorization
//  tmpdualsolver.initialSolve();
//  CoinWarmStart* initWarmStart = dualsolver.getWarmStart();
//  SolverInterface tmpdualsolver;
//  dualsolver.disableFactorization();
//  dualsolver.markHotStart();

#ifdef TRACE
  printf("\n## Attempting alternate liftings. ##\n");
#endif
  int numpivots = 0;
  for (int i = 0; i < (int) colIn.size(); i++) {
#ifdef TRACE
    printf(
        "Pivot %d\n\tcolIn: %d\n\tcolInNBIndex: %d\n\tcolOut: %d\n\toutStatus: %d\n\tminRatio: %f.\n",
        i, colIn[i], colInNBIndex[i], colOut[i], outStatus[i],
        minRatio[i]);
#endif
//    if (i != 13) {
//      continue;
//    }
    if (colOut[i] != -1 && numpivots < max_pivots && std::abs(minRatio[i]) > PIVOTEPS) {
#ifdef TRACE
      printf("\n## Generating lifting on pivot %d. ##\n", i);
#endif
      SolverInterface tmpdualsolver(Dk_init);

      if (Dk.getObjSense() != tmpdualsolver.getObjSense()) {
        cerr << "Objective senses unequal!" << endl;
        writeErrorToII("Objective senses unequal!\n", inst_info_out);
        exit(1);
      }

//      cout << "Dk obj coeff: " << Dk.getObjCoefficients()[]

//tmpdualsolver.copyParameters(Dk);

//      tmpdualsolver.setObjSense(Dk.getObjSense());

//      tmpdualsolver.enableFactorization();
//      bool accepted = tmpdualsolver.setWarmStart(initWarmStart);
//      if (!accepted) {
//        printf("*** ERROR: Warm start not accepted for pivot %d.\n", i);
//      }
//      tmpdualsolver.initialSolve();
//      tmpdualsolver.setBasisStatus(&solnInfoDual.cstat[0],
//          &solnInfoDual.rstat[0]);
//      tmpdualsolver.setColSolution(Dk.getColSolution());
//      tmpdualsolver.setRowPrice(Dk.getRowPrice());
      tmpdualsolver.enableFactorization();
      tmpdualsolver.initialSolve();

      if (tmpdualsolver.isProvenOptimal()) {
#ifdef TRACE
        printf("\n## Proven optimal after initial solve! ##\n");
#endif

        if (SHOULD_WRITE_SOLN) {
          char soln_name[300];
          snprintf(soln_name, sizeof(soln_name) / sizeof(char),
              "%s-%d", Dk_name, i);
          char OptTabName[300];
          snprintf(OptTabName, sizeof(OptTabName) / sizeof(char),
              "%s-Opt_Tableau.csv", soln_name);
          FILE* OptTabOut = fopen(OptTabName, "w");
          if (OptTabOut != NULL)
            printSimplexTableauWithNames(OptTabOut, Dk);
          else {
            char errorstring[300];
            snprintf(errorstring,
                sizeof(errorstring) / sizeof(char),
                "*** ERROR: Could not open file %s to write optimal simplex tableau.\n",
                OptTabName);
            cerr << errorstring << endl;
            writeErrorToII(errorstring, inst_info_out);
            exit(1);
          }
        }

      } else {
        char errorstring[300];
        snprintf(errorstring, sizeof(errorstring) / sizeof(char),
            "*** ERROR: NOT proven optimal after initial solve for %s for pivot %d!\n",
            Dk_name, i);
        cerr << errorstring << endl;
        char tmpname[300];
        snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
            Dk_name, i);
        writeSoln(&tmpdualsolver, tmpname);

        tmpdualsolver.writeLp(tmpname);
        tmpdualsolver.writeMps(tmpname);

        writeErrorToII(errorstring, inst_info_out);

        exit(1);
      }

      char tmpname[300];
      snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Copy-%d",
          Dk_name, i);
      writeSoln(&tmpdualsolver, tmpname);
//      tmpdualsolver.resolve();

#ifdef TRACE
      printf("Ensuring starting solution is the same for pivot %d.\n", i);
#endif
      bool succ_copy = solversEqual(Dk, tmpdualsolver, inst_info_out);
      if (!succ_copy) {
        // Let's now try forcing the basis status of everything
#ifdef TRACE
        printf(
            "First attempt at making solutions be the same failed. Trying to set cstat and rstat explicitly.\n");
#endif
        tmpdualsolver.setBasisStatus(&solnInfoDual.cstat[0],
            &solnInfoDual.rstat[0]);
        tmpdualsolver.setColSolution(Dk.getColSolution());
        tmpdualsolver.setRowPrice(Dk.getRowPrice());

        tmpdualsolver.initialSolve();

        succ_copy = solversEqual(Dk, tmpdualsolver, inst_info_out);
        if (!succ_copy) {
          char errorstring[300];
          snprintf(errorstring, sizeof(errorstring) / sizeof(char),
              "*** ERROR: Pivot %d for %s does not have initial solution same as before.\n",
              i, Dk_name);
          cerr << errorstring << endl;

          char tmpname[300];
          snprintf(tmpname, sizeof(tmpname) / sizeof(char),
              "%s-Error-%d", Dk_name, i);
          writeSoln(&tmpdualsolver, tmpname);

          tmpdualsolver.writeLp(tmpname);

          writeErrorToII(errorstring, inst_info_out);

          exit(1);
        }
      }

#ifdef TRACE
      printf("\n## Successful copy created! ##\n");
#endif

      // This is an actual pivot we can perform
//      dualsolver.disableFactorization();
      tmpdualsolver.enableSimplexInterface(true);
      if (tmpdualsolver.isProvenOptimal()) {
#ifdef TRACE
        printf(
            "\n## Proven optimal after enable simplex interface! ##\n");
#endif
      } else {
        char errorstring[300];
        snprintf(errorstring, sizeof(errorstring) / sizeof(char),
            "*** ERROR: Not proven optimal after enable simplex interface for pivot %d for %s.\n",
            i, Dk_name);
        cerr << errorstring << endl;

        char tmpname[300];
        snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
            Dk_name, i);
        writeSoln(&tmpdualsolver, tmpname);

        tmpdualsolver.writeLp(tmpname);

        writeErrorToII(errorstring, inst_info_out);

        exit(1);
      }

      // Pivot the colIn column in, colOut column out, where colOut will be at ub
      // if outStatus is 1 and at lb if outStatus is -1.
      int pivotStatus = tmpdualsolver.pivot(colIn[i], colOut[i],
          outStatus[i]);

      tmpdualsolver.disableSimplexInterface();
//      dualsolver.enableFactorization();

      if (pivotStatus == -1) {
        char errorstring[300];
        snprintf(errorstring, sizeof(errorstring) / sizeof(char),
            "*** ERROR: Pivot status, in %s, after pivoting in colIn %d and out colOut %d returned %d, which indicates singular factorization.\n",
            Dk_name, colIn[i], colOut[i], pivotStatus);
        cerr << errorstring << endl;

        // Reset basis
//        dualsolver.reset();
//
////        dualsolver.setWarmStart(initWarmStart);
////        dualsolver = tmpdualsolver;
////        dualsolver.initialSolve();
//        if (dual_ext.compare("lp") == 0)
//          dualsolver.readLp(dual_f_name.c_str());
//        else if (dual_ext.compare("mps") == 0)
//          dualsolver.readMps(dual_f_name.c_str());
//        else {
//          char errorstring[256];
//          snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//              "*** ERROR: Unrecognized extension when reading dual file: %s.\n",
//              dual_ext.c_str());
//          cerr << errorstring << endl;
//          writeErrorToII(errorstring, inst_info_out);
//          fclose(inst_info_out);
//          exit(1);
//        }
//
//#ifndef TRACE
//        dualsolver.messageHandler()->setLogLevel(0);
//        dualsolver.getModelPtr()->messageHandler()->setLogLevel(0);
//#endif
//
//        dualsolver.disableSimplexInterface();
//        dualsolver.enableFactorization();
//
//        dualsolver.setWarmStart(initWarmStart);
//
//        dualsolver.initialSolve();
        // Set basis from initial solve
//        dualsolver.setBasisStatus(solnInfoDual.cstat.data(),
//            solnInfoDual.rstat.data());
//////        dualsolver.resolve();
////
//
//        dualsolver.solveFromHotStart();

//        vector<int> tmpcstat(dualsolver.getNumCols()), tmprstat(
//            dualsolver.getNumRows());
//        dualsolver.getBasisStatus(&tmpcstat[0], &tmprstat[0]);
//
//        printf("\nCstat of variable 21 with name %s is %d. Tmpcstat: %d.\n",
//            dualsolver.getColName(21).c_str(), solnInfoDual.cstat[21],
//            tmpcstat[21]);

//        printf("\nReset to initial state after finding bad pivot %d.\n", i);
//        dualsolver.reset();
//
//        if (dual_ext.compare("lp") == 0)
//          dualsolver.readLp(dual_f_name.c_str());
//        else if (dual_ext.compare("mps") == 0)
//          dualsolver.readMps(dual_f_name.c_str());
//        else {
//          char errorstring[256];
//          snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//              "*** ERROR: Unrecognized extension when reading dual file: %s.\n",
//              dual_ext.c_str());
//          cerr << errorstring << endl;
//          writeErrorToII(errorstring, inst_info_out);
//          fclose(inst_info_out);
//          exit(1);
//        }
//
//#ifndef TRACE
//        dualsolver.messageHandler()->setLogLevel(0);
//        dualsolver.getModelPtr()->messageHandler()->setLogLevel(0);
//#endif
//
//        dualsolver.setWarmStart(initWarmStart);
//
//        dualsolver.initialSolve();
//        // Set basis from initial solve
//        dualsolver.setBasisStatus(solnInfoDual.cstat.data(),
//            solnInfoDual.rstat.data());
//        dualsolver.resolve();

//        printf("Ensuring solution is the same.\n");
//        vector<int> tmpcstat(dualsolver.getNumCols()), tmprstat(
//            dualsolver.getNumRows());
//        dualsolver.getBasisStatus(&tmpcstat[0], &tmprstat[0]);
//        for (int col = 0; col < dualsolver.getNumCols(); col++) {
//          double oldval = solnInfoDual.primalSoln[col];
//          double newval = dualsolver.getColSolution()[col];
//          if (std::abs(oldval - newval) > params.get(EPS)
//              || tmpcstat[col] != solnInfoDual.cstat[col]) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset primal values after finding singular factorization (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\nNew cstat: %d. Old cstat: %d.\n",
//                i, col, oldval, newval, oldval - newval, tmpcstat[col],
//                solnInfoDual.cstat[col]);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//        }
//        for (int row = 0; row < dualsolver.getNumRows(); row++) {
//          double oldval = solnInfoDual.slackSoln[row];
//          double newval = std::abs(
//              dualsolver.getRightHandSide()[row]
//                  - dualsolver.getRowActivity()[row]);
//          if (std::abs(oldval - newval) > params.get(EPS)
//              || tmprstat[row] != solnInfoDual.rstat[row]) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset slack values afterfinding singular factorization (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\nNew rstat: %d. Old rstat: %d.\n",
//                i, row, oldval, newval, oldval - newval, tmprstat[row],
//                solnInfoDual.rstat[row]);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//
//        }
//        for (int row = 0; row < dualsolver.getNumRows(); row++) {
//          double oldval = solnInfoDual.dualSoln[row];
//          double newval = dualsolver.getRowPrice()[row];
//          if (std::abs(oldval - newval) > params.get(EPS)) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset dual values after finding singular factorization (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\n",
//                i, row, oldval, newval, oldval - newval);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//        }

        char tmpname[300];
        snprintf(tmpname, sizeof(tmpname) / sizeof(char),
            "%s-PostSingFact-%d", Dk_name, i);
        writeSoln(&tmpdualsolver, tmpname);

        writeErrorToII(errorstring, inst_info_out);
        exit(1);
      } else {
        numpivots++;

        if (SHOULD_WRITE_SOLN) {
          char tmpname[300];
          snprintf(tmpname, sizeof(tmpname) / sizeof(char),
              "%s-PostPivot-%d", Dk_name, i);
          writeSoln(&tmpdualsolver, tmpname);
        }

        // Use this solution to lift
        liftUsingCoeffsUVSpace(liftedCut, orderedLiftedCut, cutinfo,
            subprob, solver, Dk, tmpdualsolver, deletedCols,
            oldRowIndices, oldColIndices, subspace_option,
            inst_info_out);

        // Return everything back to how it was
//        dualsolver.enableSimplexInterface(true);
//        pivotStatus = dualsolver.pivot(colOut[i], colIn[i], inStatus[i]);
//        dualsolver.disableSimplexInterface();
//        if (pivotStatus == -1) {
//          char errorstring[300];
//          snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//              "*** ERROR: Pivot status after pivoting in colOut %d and out colIn %d returned %d, which indicates singular factorization.\n",
//              colOut[i], colIn[i], pivotStatus);
//          cerr << errorstring << endl;
//          writeErrorToII(errorstring, inst_info_out);
//          exit(1);
//        }

//        printf("Returning system back to initial state.\n");
//        // Set basis from initial solve
//        dualsolver.solveFromHotStart();
//        dualsolver.setBasisStatus(solnInfoDual.cstat.data(),
//            solnInfoDual.rstat.data());
//        dualsolver.setWarmStart(initWarmStart);
//        dualsolver.resolve();
//        dualsolver.initialSolve();

//        dualsolver.reset();
//
//        if (dual_ext.compare("lp") == 0)
//          dualsolver.readLp(dual_f_name.c_str());
//        else if (dual_ext.compare("mps") == 0)
//          dualsolver.readMps(dual_f_name.c_str());
//        else {
//          char errorstring[256];
//          snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//              "*** ERROR: Unrecognized extension when reading dual file: %s.\n",
//              dual_ext.c_str());
//          cerr << errorstring << endl;
//          writeErrorToII(errorstring, inst_info_out);
//          fclose(inst_info_out);
//          exit(1);
//        }
//
//#ifndef TRACE
//        dualsolver.messageHandler()->setLogLevel(0);
//        dualsolver.getModelPtr()->messageHandler()->setLogLevel(0);
//#endif

//        dualsolver.setWarmStart(initWarmStart);
//
////        dualsolver.initialSolve();
//        // Set basis from initial solve
//        dualsolver.setBasisStatus(solnInfoDual.cstat.data(),
//            solnInfoDual.rstat.data());
//        dualsolver.resolve();

//        printf("Ensuring solution is the same.\n");
//        vector<int> tmpcstat(dualsolver.getNumCols()), tmprstat(
//            dualsolver.getNumRows());
//        dualsolver.getBasisStatus(&tmpcstat[0], &tmprstat[0]);
//        for (int col = 0; col < dualsolver.getNumCols(); col++) {
//          double oldval = solnInfoDual.primalSoln[col];
//          double newval = dualsolver.getColSolution()[col];
//          if (std::abs(oldval - newval) > params.get(EPS)
//              || tmpcstat[col] != solnInfoDual.cstat[col]) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset primal values after doing lifting using alt soln (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\nNew cstat: %d. Old cstat: %d.\n",
//                i, col, oldval, newval, oldval - newval, tmpcstat[col],
//                solnInfoDual.cstat[col]);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//        }
//        for (int row = 0; row < dualsolver.getNumRows(); row++) {
//          double oldval = solnInfoDual.slackSoln[row];
//          double newval = std::abs(
//              dualsolver.getRightHandSide()[row]
//                  - dualsolver.getRowActivity()[row]);
//          if (std::abs(oldval - newval) > params.get(EPS)
//              || tmprstat[row] != solnInfoDual.rstat[row]) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset slack values after doing lifting using alt soln (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\nNew rstat: %d. Old rstat: %d.\n",
//                i, row, oldval, newval, oldval - newval, tmprstat[row],
//                solnInfoDual.rstat[row]);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//
//        }
//        for (int row = 0; row < dualsolver.getNumRows(); row++) {
//          double oldval = solnInfoDual.dualSoln[row];
//          double newval = dualsolver.getRowPrice()[row];
//          if (std::abs(oldval - newval) > params.get(EPS)) {
//            char errorstring[300];
//            snprintf(errorstring, sizeof(errorstring) / sizeof(char),
//                "*** ERROR: Reset dual values after doing lifting using alt soln (pivot %d) not the same as old values. Variable %d has oldval: %f and newval: %f. Diff: %e.\n",
//                i, row, oldval, newval, oldval - newval);
//            cerr << errorstring << endl;
//            writeErrorToII(errorstring, inst_info_out);
//
//            char tmpname[300];
//            snprintf(tmpname, sizeof(tmpname) / sizeof(char), "%s-Error-%d",
//                dual_f_name_stub.c_str(), i);
//            writeSoln(&dualsolver, tmpname);
//
//            exit(1);
//          }
//        }
      }
    }
  }
//  printf("\n Here?\n");
//  dualsolver.unmarkHotStart();
//  dualsolver.enableFactorization();

  // Now explore if there are any alternate optimal solutions to the dual.
  // What we need to find is a row i with rstat[i] == 1 (so that the
  // corresponding dual variable u_i is non-basic, with reduced cost = 0).
  // Rc is simply b - a_i x, so rhs - row activity, or the slack s_i in that row.
  // Then we want to pivot out s_i from the primal basis, which corresponds to
  // pivoting in u_i in the dual basis. What do we pivot in for the primal?

  // Suppose u_i is nb at lower bound (can we tell whether nb at ub?).
  //

} /* liftingLandPUVSpace */

/***********************************************************************/
/**
 * Generates the entire L&P program (in the extended formulation space)
 * Num rows:
 *   m + 1 rows for the first side of the disjunction,
 *  + m + 1 rows for the second side of the disjuction,
 *  + n rows for the x = x^1 + x^2 constraints
 *  + 1 row for the lambda1 + lambda2 = 1 constraint
 *  = 2m + n + 3
 * Num cols:
 *    x, x^1, x^2, lambda1, lambda2
 *  = 3n + 2
 */
void genLandPUVSpace(SolverInterface* liftingSolver, OsiCut &cutinfo,
    SolverInterface* subsolver) {
#ifdef TRACE
  printf(
      "\n############### Starting genLandPUVSpace routine to generate lift-and-project solver. ###############\n");
#endif

  int numrows = subsolver.getNumRows(), numcols = subsolver.getNumCols();
  int splitIndex = cutinfo.splitVarIndex;
  double xk = subsolver.getColSolution()[splitIndex];
  int floorxk = floor(xk);
  int ceilxk = ceil(xk);
#ifdef TRACE
  printf("xk: %f. floorxk: %d. ceilxk: %d.\n", xk, floorxk, ceilxk);
#endif

  const double *rhs = subsolver.getRightHandSide();

  vector<double> negrhs(numrows);
  for (int row = 0; row < numrows; row++)
    negrhs[row] = -1 * rhs[row];

  vector<int> indices(numrows);
  vector<double> vec_zeroes(numrows);
  vector<double> vals;

  for (int i = 0; i < numrows; i++)
    indices[i] = i;

  CoinPackedMatrix mx;
  CoinPackedMatrix mx_zeroes;
  CoinPackedMatrix liftmx, tmpmx;

  mx.copyOf(*(subsolver.getMatrixByCol())); // Copy A into mx
  mx.transpose(); // We now have a row-ordered transposed matrix [n x m]

  // Add -1 * rhs as a row in mx
  mx.appendRow(numrows, indices.data(), negrhs.data());

  // Matrix of [n+1 x m] zeroes
  mx_zeroes.setDimensions(numcols + 1, numrows);
  //for (int col = 0; col <= numcols; col++)
  //  mx_zeroes.appendRow(numrows, indices.data(), vec_zeroes.data());

  mx_zeroes.reverseOrdering(); // Make it row-ordered

  indices.clear();
  vec_zeroes.clear();

  // First side of the disjunction
  liftmx = mx;
  liftmx.majorAppendSameOrdered(mx_zeroes); // mx_zeroes is also row ordered

  indices.push_back(splitIndex); // This is for the column -e_k u_0 and pi_0 u_0
  indices.push_back(numcols);
  vals.push_back(-1);
  vals.push_back(floorxk);
  liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
  indices.clear();
  vals.clear();

  // Second side of the disjunction
  tmpmx = mx_zeroes; // tmpmx is mx_zeroes, an n+1 x m matrix of zeroes, row-ordered
  tmpmx.majorAppendSameOrdered(mx); // mx is row-ordered too
  liftmx.minorAppendSameOrdered(tmpmx); // liftmx is row-ordered too

  indices.push_back(numcols + 1 + splitIndex); // This is for the column e_k v_0 and -(pi_0 + 1) v_0
  indices.push_back(numcols + 1 + numcols);
  vals.push_back(1);
  vals.push_back(-1 * ceilxk);
  liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
  indices.clear();
  vals.clear();

  // Add col corresponding to beta
  indices.push_back(numcols);
  indices.push_back(2 * numcols + 1);
  vals.resize(2, 1);
  liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
  indices.clear();
  vals.clear();

  // Add cols corresponding to bounds on the variables. We need these explicitly
  // because x^i <= ub * lambda_i, and x^i >= lb * lambda_i.
  for (int col = 0; col < numcols; col++) {
    double currlb = subsolver.getColLower()[col];
    double currub = subsolver.getColUpper()[col];
    vals.resize(2, 0);

    // Lower-bound (x^0 >= lb * lambda_0)
    indices.push_back(col);
    indices.push_back(numcols);
    vals[0] = 1;
    vals[1] = -1 * currlb;
    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    // Upper-bound (will be stored as -x >= -ub * lambda1)
    vals[0] = -1;
    vals[1] = currub;
    if (currub < subsolver.getInfinity() - params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    // Lower-bound (x^1 >= lb * lambda_1)
    indices[0] = numcols + 1 + col;
    indices[1] = numcols + 1 + numcols;
    vals[0] = 1;
    vals[1] = -1 * currlb;
    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    // Upper-bound (will be stored as -x >= -ub * lambda1)
    vals[0] = -1;
    vals[1] = currub;
    if (currub < subsolver.getInfinity() - params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
    indices.clear();
    vals.clear();
  }

  // Add liftmx to the liftingSolver and set right hand sides, bounds, etc.
  // How many rows and cols in the extended formulation? (Same as written above.)
//  int EFnumrows = liftmx.getMajorDim(); // Because we are row-ordered
//  int EFnumcols = liftmx.getMinorDim();

  /*
   // Bound vector
   vector<double> colLB(EFnumcols), colUB(EFnumcols);
   for (int i = 0; i < 2*(numrows + 1); i++) {
   colLB[i] = -DBL_MAX;
   colUB[i] = 10.0;
   }

   liftingSolver.loadProblem(liftmx, colLB.data(), colUB.data(), NULL, NULL, NULL, NULL);
   */
  liftingSolver.loadProblem(liftmx, NULL, NULL, NULL, NULL, NULL, NULL);

  // Set col bounds (decide which are >=, which are <=)
  for (int row = 0; row < numrows; row++) {
    char rowsense = subsolver.getRowSense()[row];
    if (rowsense == 'E') {
      liftingSolver.setColBounds(row, -1.0 * liftingSolver.getInfinity(),
          liftingSolver.getInfinity());
      liftingSolver.setColBounds(numrows + 1 + row,
          -1 * liftingSolver.getInfinity(),
          liftingSolver.getInfinity());
    }
    if (rowsense == 'L') {
      // Dual variable would be negative, so let us negate all the coefficients of the row
//      for (int col = 0; col < EFnumcols; col++) {
//        CoinPackedMatrix *mutable_mx = liftingSolver.getMutableMatrixByCol();
//        mutable_mx->modifyCoefficient(row, col, -1*mutable_mx->getCoefficient(row, col));
//        mutable_mx->modifyCoefficient(numrows + 1 + row, col, -1*mutable_mx->getCoefficient(row, col));
//      }
      liftingSolver.setColBounds(row, -1.0 * liftingSolver.getInfinity(),
          0.0);
      liftingSolver.setColBounds(numrows + 1 + row,
          -1 * liftingSolver.getInfinity(), 0.0);
    }
  }

  // // Set bounds for x = x^0 + x^1 rows. One such row per variable.
  // Also set bounds for x variables (make them free)
  for (int col = 0; col < numcols; col++) {

    // Make $x$ variables free
    liftingSolver.setRowBounds(col, cutinfo[col], cutinfo[col]); // x_i^0
    liftingSolver.setRowBounds(numcols + 1 + col, cutinfo[col],
        cutinfo[col]); // x_i^1
    //liftingSolver.setColBounds(currCol, -1 * liftingSolver.getInfinity(),
    //    liftingSolver.getInfinity()); // x_i
  }
  //liftingSolver.setRowBounds(2 * (numrows + 1) + numcols, 1.0, 1.0); // Set bound for lambda1 + lambda2 row

  // Also we need that the beta rows are <= 0
  liftingSolver.setRowBounds(numcols, -1.0 * liftingSolver.getInfinity(),
      0.0);
  liftingSolver.setRowBounds(2 * numcols + 1,
      -1 * liftingSolver.getInfinity(), 0.0);

//  // beta should be equal to the previous rhs (had constraint for this, but we removed this)
//  indices.push_back(2 * (numrows + 1));
//  vals.push_back(1);
//  liftingSolver.addRow(1, indices.data(), vals.data(),
//      -1 * liftingSolver.getInfinity(), cutinfo.RHS);
//  indices.clear();
//  vals.clear();

  // Set beta to be urs
  liftingSolver.setColBounds(2 * (numrows + 1),
      -1 * liftingSolver.getInfinity(), liftingSolver.getInfinity());

  // Set objective
  liftingSolver.setObjCoeff(2 * (numrows + 1), -1.0);
  liftingSolver.setObjSense(1.0);
}

/***********************************************************************/
/**
 * Lift using the current optimal solution to dualsolver
 */
void liftUsingCoeffsUVSpace(OsiCuts &liftedCut, OsiCuts &orderedLiftedCut,
    OsiCut &cutinfo, SolverInterface* subprob,
    SolverInterface* solver, SolverInterface* Dk,
    SolverInterface* dualsolver, vector<int> &deletedCols,
    vector<int> oldRowIndices, std::vector<int> &oldColIndices,
    int subspace_option, FILE* inst_info_out, const Parameters& params) {
#ifdef TRACE
  printf("\n## Starting liftUsingCoeffsUVSpace routine. ##\n");
#endif

//  double* dualsoln = dualsolver.getColSolution();

  std::vector<double> DkDual(dualsolver.getColSolution(),
      dualsolver.getColSolution() + dualsolver.getNumCols());
//  double DkOpt;

  // Use this solution to lift
//  orderDualVars(PkDual, subprob, dualsolver, Dk);
//  PkDual.assign(dualsolver.getColSolution(),
//      dualsolver.getColSolution() + dualsolver.getNumCols());

  // It should be true that the objective we get
  // is at least the rhs
  // Recall that we are solving max beta
  // but that Clp changes that to min -beta
  // So dualsolver.getObjValue() is negative beta.
  // Thus, beta^* should be >= cutinfo.RHS,
  // so -1*dualsolver.getObjValue - cutinfo.RHS >= -params.get(EPS)
  // i.e., cutinfo.RHS + dualsolver.getObjValue <= params.get(EPS)
  if (cutinfo.rhs() + dualsolver.getObjValue() > params.get(EPS)) {
    error_msg(errorstring,
        "When checking OsiCut on var %d,\n"
        "\t(1) solver objective: %f\n"
        "\t(2) rhs of the OsiCut we are lifting: %f\n(2) should be <= (1).\n",
        cutinfo.splitVarIndex, dualsolver.getObjValue(), cutinfo.RHS);
    writeErrorToLog(errorstring, params.logfile);
    exit(1);
  }

  // Lift coefficients
  OsiCut tmpDualLiftedCut, tmpDualOrderedLiftedCut;
  liftLandP(tmpDualLiftedCut, solver, deletedCols, oldRowIndices, DkDual,
      cutinfo, subspace_option, inst_info_out);

  // Change rhs
  tmpDualLiftedCut.setLb(cutinfo.rhs()); // Because Clp switches it from max to min

  // For lifted OsiCut in order of original variables
  orderLiftedCut(tmpDualLiftedCut, &solver, &subprob, deletedCols,
      tmpDualOrderedLiftedCut);

  // Print this lifted OsiCut
#ifdef TRACE
  printf("\n## Ordered lifted OsiCut information: ##\n");
  printf("Split var index: %d\n", tmpDualOrderedLiftedCut.splitVarIndex);
  printf("Split var name: %s\n",
      solver.getColName(tmpDualOrderedLiftedCut.splitVarIndex).c_str());
  printf("Rhs: %f\n", tmpDualOrderedLiftedCut.RHS);
  printf("%s,%s\n", "Orig var", "Coeff");
  for (int k = 0; k < (int) tmpDualOrderedLiftedCut.coeff.size(); k++)
  printf("%d,%f\n", k, tmpDualOrderedLiftedCut.coeff[k]);
  printf("\n");
#endif

  // Correct the coefficients, because some may be incorrect due to complemented vars
  complementCoeffs(tmpDualOrderedLiftedCut, deletedCols, solver);

  // Re-normalize
  double absrhs = std::abs(tmpDualOrderedLiftedCut.rhs());
  if (absrhs > params.get(EPS)) {
    for (int col = 0; col < (int) tmpDualOrderedLiftedCut.coeff.size();
        col++) {
      tmpDualOrderedLiftedCut.coeff[col] =
          tmpDualOrderedLiftedCut.coeff[col] / absrhs;
    }
    if (tmpDualOrderedLiftedCut.RHS > params.get(EPS)) {
      tmpDualOrderedLiftedCut.RHS = 1.0;
    } else if (tmpDualOrderedLiftedCut.RHS < params.get(EPS)) {
      tmpDualOrderedLiftedCut.RHS = -1.0;
    }
  } else {
    tmpDualOrderedLiftedCut.RHS = 0.0;
  }

  // Print this ordered lifted OsiCut after the complementing
#ifdef TRACE
  printf("\n## After complementing, ordered lifted OsiCut information: ##\n");
  printf("Split var index: %d\n", tmpDualOrderedLiftedCut.splitVarIndex);
  printf("Split var name: %s\n",
      solver.getColName(tmpDualOrderedLiftedCut.splitVarIndex).c_str());
  printf("Rhs: %f\n", tmpDualOrderedLiftedCut.RHS);
  printf("%s,%s\n", "Orig var", "Coeff");
  for (int k = 0; k < (int) tmpDualOrderedLiftedCut.coeff.size(); k++)
  printf("%d,%f\n", k, tmpDualOrderedLiftedCut.coeff[k]);
  printf("\n");
#endif

  // Check if this OsiCut exists in the previously generated cuts,
  // and if it doesn't, store the OsiCut
  bool newcut = true;
  for (int i = 0; i < (int) orderedLiftedCut.size(); i++) {
    if (!cutsDifferent(tmpDualOrderedLiftedCut, orderedLiftedCut[i])) {
#ifdef TRACE
      printf(
          "This OsiCut is old, so we are not storing it. Previous OsiCut %d for this split is the same.\n",
          i);
#endif
      newcut = false;
      break;
    }
  }

  if (newcut) {
#ifdef TRACE
    printf("This OsiCut is new, so we are storing it.\n");
#endif
    liftedCut.push_back(tmpDualLiftedCut);
    orderedLiftedCut.push_back(tmpDualOrderedLiftedCut);
  }
}

/*******************************************************/
/** Lifting in the (u,v)-variable space, separate LPs **/
/*******************************************************/

/***********************************************************************/
/**
 * Set up split lifting LPs in the u, v space.
 */
void setupLandPSeparate(SolverInterface* Dku, SolverInterface* Dkv,
    SolverInterface* subprob, OsiCut &cutinfo, int currvar) {
#ifdef TRACE
  printf("\n## Setting up Dku and Dkv LPs. ##\n");
#endif

  genLandPU(Dku, cutinfo, subprob);
  genLandPV(Dkv, cutinfo, subprob);

#ifdef TRACE
  printf("D%du: # rows: %d. # cols: %d.\n", currvar, Dku.getNumRows(),
      Dku.getNumCols());
  printf("D%dv: # rows: %d. # cols: %d.\n", currvar, Dkv.getNumRows(),
      Dkv.getNumCols());
  printf("\n");
#endif
}

/**
 * Update split lifting LPs in the u, v space
 * with info for the new split.
 */
void updateLandPSeparate(SolverInterface* Dku, SolverInterface* Dkv,
    OsiCut &cutinfo, const int currvar, const int oldSplitIndex, 
    const double floorxk, const double ceilxk, 
    const int numRows, const int numCols) {
  // First the right sides
  for (int col = 0; col < (int) cutinfo.coeff.size(); col++) {
    double alphaj = cutinfo.coeff[col];
    Dku.setRowBounds(col, alphaj, alphaj); // x_j
    Dkv.setRowBounds(col, alphaj, alphaj); // x_j
  }

  // Now change the objective coefficient for column u0 and v0
  Dku.setObjCoeff(numSubRows, floorxk);
  Dkv.setObjCoeff(numSubRows, -1 * ceilxk);

  // We also switch so that column u_0 has the 1 in the right row
  // and v_0 has a -1 in the correct row
  Dku.modifyCoefficient(oldSplitIndex, numSubRows, 0.0);
  Dkv.modifyCoefficient(oldSplitIndex, numSubRows, 0.0);
  Dku.modifyCoefficient(currvar, numSubRows, -1.0);
  Dkv.modifyCoefficient(currvar, numSubRows, 1.0);
} /* updateLandPSeparate */

/***********************************************************************/
/**
 * Generates lifting in the L&P style.
 */
void liftingLandPSeparate(OsiCuts &liftedCut, OsiCuts &orderedLiftedCut,
    OsiCut &cutinfo, SolverInterface* subprob,
    SolverInterface* solver, std::vector<int> &deletedCols,
    std::vector<int> &oldRowIndices, std::vector<int> &oldColIndices,
    int max_pivots, const std::string &out_f_name_stub,
    const std::string &cutname, SolverInterface* Dku,
    SolverInterface* Dkv, int subspace_option, FILE *inst_info_out) {

//  int currsplit = cutinfo.splitVarIndex;

  /***********************************************************************************
   * Generate LP
   ***********************************************************************************/

  // Save Dku and Dkv for posterity
  char Dku_name[256];
  snprintf(Dku_name, sizeof(Dku_name) / sizeof(char), "%s-D%su",
      out_f_name_stub.c_str(), cutname.c_str());
  Dku.setObjName(Dku_name);
  if (WRITE_LP_FILE) {
    Dku.writeLp(Dku_name);
    Dku.writeMps(Dku_name, "mps");
  }

  char Dkv_name[256];
  snprintf(Dkv_name, sizeof(Dkv_name) / sizeof(char), "%s-D%sv",
      out_f_name_stub.c_str(), cutname.c_str());
  Dkv.setObjName(Dkv_name);
  if (WRITE_LP_FILE) {
    Dkv.writeLp(Dkv_name);
    Dkv.writeMps(Dkv_name, "mps");
  }

  /***********************************************************************************
   * Perform lifting using first solution
   ***********************************************************************************/
#ifdef TRACE
  printf("\n## Obtaining first lifting. ##\n");
#endif

#ifdef TRACE
  printf("\n## Solving Dku. ##\n");
#endif
  Dku.enableFactorization();
  Dku.initialSolve();

  if (Dku.isProvenOptimal()) {
#ifdef TRACE
    printf("Dku opt value: %f\n", Dku.getObjValue());
#endif
  } else {
    char errorstring[300];
    snprintf(errorstring, sizeof(errorstring) / sizeof(char),
        "%s not proven optimal!\n", Dku_name);
    cerr << errorstring << endl;
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }
#ifdef TRACE
  printf("Saving dual LP solution information.\n");
#endif
  SolutionInfo solnInfoDualU(Dku, false);
  if (SHOULD_WRITE_BRIEF_SOLN)
    solnInfoDualU.printBriefSolutionInfo(Dku, Dku_name);

  if (SHOULD_WRITE_SOLN) {
    char soln_name[300];
    snprintf(soln_name, sizeof(soln_name) / sizeof(char), "%sInit",
        Dku_name);
    writeSoln(&Dku, soln_name);
    char OptTabName[300];
    snprintf(OptTabName, sizeof(OptTabName) / sizeof(char),
        "%s-Opt_Tableau.csv", soln_name);
    FILE* OptTabOut = fopen(OptTabName, "w");
    if (OptTabOut != NULL)
      printSimplexTableauWithNames(OptTabOut, Dku);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write optimal simplex tableau.\n",
          OptTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }

    char NBTabName[300];
    snprintf(NBTabName, sizeof(NBTabName) / sizeof(char),
        "%s-NB_Opt_Tableau.csv", soln_name);
    FILE* NBTabOut = fopen(NBTabName, "w");
    if (NBTabOut != NULL)
      printNBSimplexTableauWithNames(NBTabOut, Dku);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write NB optimal simplex tableau.\n",
          NBTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }
  }

  // Solve Dkv
#ifdef TRACE
  printf("\n## Solving Dkv. ##\n");
#endif
  Dkv.enableFactorization();
  Dkv.initialSolve();

  if (Dkv.isProvenOptimal()) {
#ifdef TRACE
    printf("Dkv opt value: %f\n", Dkv.getObjValue());
#endif
  } else {
    char errorstring[300];
    snprintf(errorstring, sizeof(errorstring) / sizeof(char),
        "%s not proven optimal!\n", Dkv_name);
    cerr << errorstring << endl;
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }

#ifdef TRACE
  printf("Saving dual LP solution information.\n");
#endif
  SolutionInfo solnInfoDualV(Dkv, false);
  solnInfoDualV.printBriefSolutionInfo(Dkv, Dkv_name);

  if (SHOULD_WRITE_SOLN) {
    char soln_name[300];
    snprintf(soln_name, sizeof(soln_name) / sizeof(char), "%sInit",
        Dkv_name);
    writeSoln(&Dkv, soln_name);
    char OptTabName[300];
    snprintf(OptTabName, sizeof(OptTabName) / sizeof(char),
        "%s-Opt_Tableau.csv", soln_name);
    FILE* OptTabOut = fopen(OptTabName, "w");
    if (OptTabOut != NULL)
      printSimplexTableauWithNames(OptTabOut, Dkv);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write optimal simplex tableau.\n",
          OptTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }

    char NBTabName[300];
    snprintf(NBTabName, sizeof(NBTabName) / sizeof(char),
        "%s-NB_Opt_Tableau.csv", soln_name);
    FILE* NBTabOut = fopen(NBTabName, "w");
    if (NBTabOut != NULL)
      printNBSimplexTableauWithNames(NBTabOut, Dkv);
    else {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "*** ERROR: Could not open file %s to write NB optimal simplex tableau.\n",
          NBTabName);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }
  }

  liftUsingCoeffsSplit(liftedCut, orderedLiftedCut, cutinfo, subprob, solver,
      Dku, Dkv, deletedCols, oldRowIndices, oldColIndices,
      subspace_option, inst_info_out);
}

/***********************************************************************/
/**
 * Generates the L&P program (in the extended formulation space)
 * for the u variables only
 *
 * We have the problem
 *
 * min \bar\alpha^T x
 * s.t.
 *     Ax \ge b         (u)
 *     -x_k \ge -\pi_0  (u_0)
 *      x \ge \ell      (\lambda)
 *     -x \ge -g        (\mu)
 *
 * Dual is
 *
 * max b^T u - pi_0 u_0 + \ell^T \lambda - g^T \mu
 *     (same as)
 * min -b^T u + pi_0 u_0 - \ell^T \lambda + g^T \mu
 * s.t.
 *     A^T u - e_k u_0 + \lambda - \mu = \bar\alpha
 *     u, u_0 \ge 0 (depends actually on rows)
 *
 */
void genLandPU(SolverInterface* liftingSolver, OsiCut &cutinfo,
    SolverInterface* subsolver) {
#ifdef TRACE
  printf(
      "\n############### Starting genLandPU routine to generate lift-and-project solver. ###############\n");
#endif

  int numrows = subsolver.getNumRows(), numcols = subsolver.getNumCols();
  int splitIndex = cutinfo.splitVarIndex;
  double xk = subsolver.getColSolution()[splitIndex];
  int floorxk = floor(xk);
  int ceilxk = ceil(xk);
#ifdef TRACE
  printf("xk: %f. floorxk: %d. ceilxk: %d.\n", xk, floorxk, ceilxk);
#endif

  const double *rhs = subsolver.getRightHandSide();

  vector<int> indices;
  vector<double> vals;

  CoinPackedMatrix mx;
  CoinPackedMatrix liftmx;

  mx.copyOf(*(subsolver.getMatrixByCol())); // Copy A into mx
  mx.transpose(); // We now have a row-ordered transposed matrix [n x m]

  indices.clear();

  // A^T u
  liftmx = mx;
//  liftmx.majorAppendSameOrdered(mx_zeroes); // mx_zeroes is also row ordered

  // -e_k u_0
  indices.push_back(splitIndex);
  vals.push_back(-1);
  liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
  indices.clear();
  vals.clear();

  // Add cols corresponding to bounds on the variables.
  // We need these explicitly
  // because x^i <= ub (i.e., -x^i >= -ub),
  // and x^i >= lb.
  for (int col = 0; col < numcols; col++) {
    double currlb = subsolver.getColLower()[col];
    double currub = subsolver.getColUpper()[col];
    vals.resize(1, 0);

    // Lower-bound (+ \ell)
    indices.push_back(col);
    vals[0] = 1;
    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    // Upper-bound (- g)
    vals[0] = -1;
    if (currub < subsolver.getInfinity() - params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    indices.clear();
    vals.clear();
  }

  // Add liftmx to the liftingSolver and set right hand sides, bounds, etc.
  // How many rows and cols in the extended formulation? (Same as written above.)
//  int EFnumrows = liftmx.getMajorDim(); // Because we are row-ordered
//  int EFnumcols = liftmx.getMinorDim();

  liftingSolver.loadProblem(liftmx, NULL, NULL, NULL, NULL, NULL, NULL);

  // Set col bounds (decide which are >=, which are <=)
  for (int row = 0; row < numrows; row++) {
    char rowsense = subsolver.getRowSense()[row];
    if (rowsense == 'E') {
      liftingSolver.setColBounds(row, -1 * liftingSolver.getInfinity(),
          liftingSolver.getInfinity());
//      liftingSolver.setColBounds(numrows + 1 + row,
//          -1 * liftingSolver.getInfinity(), liftingSolver.getInfinity());
    }
    if (rowsense == 'L') {
      liftingSolver.setColBounds(row, -1 * liftingSolver.getInfinity(),
          0.0);
//      liftingSolver.setColBounds(numrows + 1 + row,
//          -1 * liftingSolver.getInfinity(), 0.0);
    }
  }

  // Also set bounds for x variables (make them free)
  for (int col = 0; col < numcols; col++) {
    liftingSolver.setRowBounds(col, cutinfo.coeff[col], cutinfo.coeff[col]); // x_i^0
//    liftingSolver.setRowBounds(numcols + 1 + col, cutinfo.coeff[col],
//        cutinfo.coeff[col]); // x_i^1
  }

  // Set objective
  // -b^T u + pi_0 u_0 - \ell^T \lambda + g^T \mu
  for (int row = 0; row < numrows; row++) {
    liftingSolver.setObjCoeff(row, -1 * rhs[row]);
  }
  liftingSolver.setObjCoeff(numrows, floorxk);
  int i = 1;
  for (int col = 0; col < numcols; col++) {
    double currlb = subsolver.getColLower()[col];
    double currub = subsolver.getColUpper()[col];

    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS)) {
      liftingSolver.setObjCoeff(numrows + i, -1 * currlb);
      i++;
    }

    if (currub < subsolver.getInfinity() - params.get(EPS)) {
      liftingSolver.setObjCoeff(numrows + i, currub);
      i++;
    }
  }
  liftingSolver.setObjSense(1.0);

  // Add constraint corresponding to u0 <= M
  // First calculate M, which we can say is 1000 * \sum_{j \in S} \bar\alpha_j \bar y_j
  if (MMULT != 0) {
    double M = 0.0;
    for (int col = 0; col < numcols; col++) {
      M += cutinfo.coeff[col] * subsolver.getColSolution()[col];
    }

    M *= MMULT;

    indices.push_back(numrows);
    vals.push_back(1.0);
    liftingSolver.addRow((int) indices.size(), indices.data(), vals.data(),
        -1 * liftingSolver.getInfinity(), M);
    indices.clear();
    vals.clear();
  }
}

/***********************************************************************/
/**
 * Generates the L&P program (in the extended formulation space)
 * for the v variables only
 *
 * We have the problem
 *
 * min \bar\alpha^T x
 * s.t.
 *     Ax \ge b         (v)
 *    x_k \ge ceilxk    (v_0)
 *      x \ge \ell      (\lambda)
 *     -x \ge -g        (\mu)
 *
 * Dual is
 *
 * max b^T v + ceilxk v_0 + \ell^T \lambda - g^T \mu
 *     (same as)
 * min -b^T v - ceilxk v_0 - \ell^T \lambda + g^T \mu
 * s.t.
 *     A^T v + e_k v_0 + \lambda - \mu = \bar\alpha
 *     v, v_0 \ge 0 (depends actually on rows)
 *
 */
void genLandPV(SolverInterface* liftingSolver, OsiCut &cutinfo,
    SolverInterface* subsolver) {
#ifdef TRACE
  printf(
      "\n############### Starting genLandPV routine to generate lift-and-project solver. ###############\n");
#endif

  int numrows = subsolver.getNumRows(), numcols = subsolver.getNumCols();
  int splitIndex = cutinfo.splitVarIndex;
  double xk = subsolver.getColSolution()[splitIndex];
  int floorxk = floor(xk);
  int ceilxk = ceil(xk);
#ifdef TRACE
  printf("xk: %f. floorxk: %d. ceilxk: %d.\n", xk, floorxk, ceilxk);
#endif

  const double *rhs = subsolver.getRightHandSide();

  vector<int> indices;
  vector<double> vals;

  CoinPackedMatrix mx;
  CoinPackedMatrix liftmx;

  mx.copyOf(*(subsolver.getMatrixByCol())); // Copy A into mx
  mx.transpose(); // We now have a row-ordered transposed matrix [n x m]

  indices.clear();

  // A^T v
  liftmx = mx;
//  liftmx.majorAppendSameOrdered(mx_zeroes); // mx_zeroes is also row ordered

  // + e_k v_0
  indices.push_back(splitIndex);
  vals.push_back(1);
  liftmx.appendCol((int) indices.size(), indices.data(), vals.data());
  indices.clear();
  vals.clear();

  // Add cols corresponding to bounds on the variables.
  // We need these explicitly
  // because x^i <= ub (i.e., -x^i >= -ub),
  // and x^i >= lb.
  for (int col = 0; col < numcols; col++) {
    double currlb = subsolver.getColLower()[col];
    double currub = subsolver.getColUpper()[col];
    vals.resize(1, 0);

    // Lower-bound (+ \ell)
    indices.push_back(col);
    vals[0] = 1;
    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    // Upper-bound (- g)
    vals[0] = -1;
    if (currub < subsolver.getInfinity() - params.get(EPS))
      liftmx.appendCol((int) indices.size(), indices.data(), vals.data());

    indices.clear();
    vals.clear();
  }

  // Add liftmx to the liftingSolver and set right hand sides, bounds, etc.
  // How many rows and cols in the extended formulation? (Same as written above.)
//  int EFnumrows = liftmx.getMajorDim(); // Because we are row-ordered
//  int EFnumcols = liftmx.getMinorDim();

  liftingSolver.loadProblem(liftmx, NULL, NULL, NULL, NULL, NULL, NULL);

  // Set col bounds (decide which are >=, which are <=)
  for (int row = 0; row < numrows; row++) {
    char rowsense = subsolver.getRowSense()[row];
    if (rowsense == 'E') {
      liftingSolver.setColBounds(row, -1 * liftingSolver.getInfinity(),
          liftingSolver.getInfinity());
//      liftingSolver.setColBounds(numrows + 1 + row,
//          -1 * liftingSolver.getInfinity(), liftingSolver.getInfinity());
    }
    if (rowsense == 'L') {
      liftingSolver.setColBounds(row, -1 * liftingSolver.getInfinity(),
          0.0);
//      liftingSolver.setColBounds(numrows + 1 + row,
//          -1 * liftingSolver.getInfinity(), 0.0);
    }
  }

  // Also set bounds for x variables (make them free)
  for (int col = 0; col < numcols; col++) {
    liftingSolver.setRowBounds(col, cutinfo.coeff[col], cutinfo.coeff[col]); // x_i^0
//    liftingSolver.setRowBounds(numcols + 1 + col, cutinfo.coeff[col],
//        cutinfo.coeff[col]); // x_i^1
  }

  // Set objective
  // -b^T v - ceilxk v_0 - \ell^T \lambda + g^T \mu
  for (int row = 0; row < numrows; row++) {
    liftingSolver.setObjCoeff(row, -1 * rhs[row]);
  }
  liftingSolver.setObjCoeff(numrows, -1 * ceilxk);
  int i = 1;
  for (int col = 0; col < numcols; col++) {
    double currlb = subsolver.getColLower()[col];
    double currub = subsolver.getColUpper()[col];

    if (currlb > -1 * subsolver.getInfinity() + params.get(EPS)) {
      liftingSolver.setObjCoeff(numrows + i, -1 * currlb);
      i++;
    }

    if (currub < subsolver.getInfinity() - params.get(EPS)) {
      liftingSolver.setObjCoeff(numrows + i, currub);
      i++;
    }
  }
  liftingSolver.setObjSense(1.0);

  // Add constraint corresponding to u0 <= M
  // First calculate M, which we can say is 1000 * \sum_{j \in S} \bar\alpha_j \bar y_j
  if (MMULT != 0) {
    double M = 0.0;
    for (int col = 0; col < numcols; col++) {
      M += cutinfo.coeff[col] * subsolver.getColSolution()[col];
    }

    M *= MMULT;

    indices.push_back(numrows);
    vals.push_back(1.0);
    liftingSolver.addRow((int) indices.size(), indices.data(), vals.data(),
        -1 * liftingSolver.getInfinity(), M);
    indices.clear();
    vals.clear();
  }
}

/***********************************************************************/
/**
 * Lift using the current optimal solution to dualsolvers Dku and Dkv
 */
void liftUsingCoeffsSplit(std::OsiCuts &liftedCut,
    std::OsiCuts &orderedLiftedCut, OsiCut &cutinfo,
    SolverInterface* subprob, SolverInterface* solver,
    SolverInterface* Dku, SolverInterface* Dkv,
    std::vector<int> &deletedCols, std::vector<int> oldRowIndices,
    std::vector<int> &oldColIndices, int subspace_option,
    FILE* inst_info_out) {
#ifdef TRACE
  printf("\n## Starting liftUsingCoeffsUVSpace routine. ##\n");
#endif

//  double* dualsoln = dualsolver.getColSolution();

  int numSubRows = subprob.getNumRows();
  std::vector<double> DkDual(2 * (numSubRows + 1));

  // We first set the u and u0 solution variables:
  for (int col = 0; col <= numSubRows; col++) {
    DkDual[col] = Dku.getColSolution()[col];
  }

  // Next we set the v and v0 solution variables:
  for (int col = 0; col <= numSubRows; col++) {
    DkDual[numSubRows + 1 + col] = Dkv.getColSolution()[col];
  }

  // Use this solution to lift

  // It should be true that the objective we get
  // is at least the rhs
  // Recall that we are solving min -beta (essentially)
  // So Dk_.getObjValue() is negative beta.
  // Thus, beta^* should be >= cutinfo.RHS,
  // so -1*DkOpt - cutinfo.RHS >= -params.get(EPS)
  // i.e. DkOpt + cutinfo.RHS <= params.get(EPS)

  double DkOpt;

  if (Dku.getObjValue() < Dkv.getObjValue()) {
    DkOpt = Dkv.getObjValue();
  } else {
    DkOpt = Dku.getObjValue();
  }

  if (cutinfo.RHS + DkOpt > params.get(EPS)) {
    char errorstring[300];
    snprintf(errorstring, sizeof(errorstring) / sizeof(char),
        "*** ERROR: \n\t(1) solver objective: %f\n\t(2) rhs of the OsiCut we are lifting: %f\n(2) should be <= -(1).\n",
        DkOpt, cutinfo.RHS);
    cerr << errorstring << endl;
    writeErrorToII(errorstring, inst_info_out);
    exit(1);
  }

  // Lift coefficients
  OsiCut tmpDualLiftedCut, tmpDualOrderedLiftedCut;
  liftLandP(tmpDualLiftedCut, solver, deletedCols, oldRowIndices, DkDual,
      cutinfo, subspace_option, inst_info_out);

  // Change rhs
  tmpDualLiftedCut.RHS = cutinfo.RHS; // Because Clp switches it from max to min

  // For lifted OsiCut in order of original variables
  orderLiftedCut(tmpDualLiftedCut, &solver, &subprob, deletedCols,
      tmpDualOrderedLiftedCut);

  // Print this lifted OsiCut
#ifdef TRACE
  printf("\n## Ordered lifted OsiCut information: ##\n");
  printf("Split var index: %d\n", tmpDualOrderedLiftedCut.splitVarIndex);
  printf("Split var name: %s\n",
      solver.getColName(tmpDualOrderedLiftedCut.splitVarIndex).c_str());
  printf("Rhs: %f\n", tmpDualOrderedLiftedCut.RHS);
  printf("%s,%s\n", "Orig var", "Coeff");
  for (int k = 0; k < (int) tmpDualOrderedLiftedCut.coeff.size(); k++)
  printf("%d,%f\n", k, tmpDualOrderedLiftedCut.coeff[k]);
  printf("\n");
#endif

  // Correct the coefficients, because some may be incorrect due to complemented vars
  complementCoeffs(tmpDualOrderedLiftedCut, deletedCols, solver);

  // Re-normalize
  double absrhs = std::abs(tmpDualOrderedLiftedCut.RHS);
  if (absrhs > params.get(EPS)) {
    for (int col = 0; col < (int) tmpDualOrderedLiftedCut.coeff.size();
        col++) {
      tmpDualOrderedLiftedCut.coeff[col] =
          tmpDualOrderedLiftedCut.coeff[col] / absrhs;
    }
    if (tmpDualOrderedLiftedCut.RHS > params.get(EPS)) {
      tmpDualOrderedLiftedCut.RHS = 1.0;
    } else if (tmpDualOrderedLiftedCut.RHS < params.get(EPS)) {
      tmpDualOrderedLiftedCut.RHS = -1.0;
    }
  } else {
    tmpDualOrderedLiftedCut.RHS = 0.0;
  }

  // Print this ordered lifted OsiCut after the complementing
#ifdef TRACE
  printf("\n## After complementing, ordered lifted OsiCut information: ##\n");
  printf("Split var index: %d\n", tmpDualOrderedLiftedCut.splitVarIndex);
  printf("Split var name: %s\n",
      solver.getColName(tmpDualOrderedLiftedCut.splitVarIndex).c_str());
  printf("Rhs: %f\n", tmpDualOrderedLiftedCut.RHS);
  printf("%s,%s\n", "Orig var", "Coeff");
  for (int k = 0; k < (int) tmpDualOrderedLiftedCut.coeff.size(); k++)
  printf("%d,%f\n", k, tmpDualOrderedLiftedCut.coeff[k]);
  printf("\n");
#endif

  // Check if this OsiCut exists in the previously generated cuts,
  // and if it doesn't, store the OsiCut
  bool newcut = true;
  for (int i = 0; i < (int) orderedLiftedCut.size(); i++) {
    if (!cutsDifferent(tmpDualOrderedLiftedCut, orderedLiftedCut[i])) {
#ifdef TRACE
      printf(
          "This OsiCut is old, so we are not storing it. Previous OsiCut %d for this split is the same.\n",
          i);
#endif
      newcut = false;
      break;
    }
  }

  if (newcut) {
#ifdef TRACE
    printf("This OsiCut is new, so we are storing it.\n");
#endif
    liftedCut.push_back(tmpDualLiftedCut);
    orderedLiftedCut.push_back(tmpDualOrderedLiftedCut);
  }
}

/*************************************/
/** Lifting methods used generally  **/
/*************************************/

/***********************************************************************/
/**
 * @brief Lift coefficients from subspace to full space.
 * @param solver        ::  Original problem. Call this (LP).
 * @param deletedCols   ::  Indices of columns deleted from (LP).
 * @param oldRowIndices ::  Row i in subprob is oldRowIndices[i] of (LP).
 * @param DDual         ::  Dual solution to D.
 * @param cutinfo       ::  Info about the OsiCut generated in the subspace.
 * @param liftedCut     ::  Coefficients to be transmitted about lifted OsiCut.
 */
void liftLandP(OsiCut &liftedCut, SolverInterface* solver,
    std::vector<int> &deletedCols, std::vector<int> &oldRowIndices,
    std::vector<double> &DDual, OsiCut &cutinfo, int subspace_option,
    FILE* inst_info_out) {
#ifdef TRACE
  std::cout
  << "\n############### Starting liftLandP routine to generate lifted coefficients. ###############"
  << std::endl;
#endif

  // Get pointer to elements from original problem
  const CoinPackedMatrix *A = solver.getMatrixByCol();

  // # original columns, subspace columns,
  // rows of A that were met at equality for Au and Av
  int numOrigCols = solver.getNumCols();
//  int numOrigRows = solver.getNumRows();
  int numSubRows = (int) oldRowIndices.size();
//  int numSubCols = numOrigCols - (int) deletedCols.size();

  // Reserve space in liftedCut std::vector
  liftedCut.coeff.reserve(numOrigCols);

  // Coefficients for subspace variables remain the same
  // Note here we are also include the split index and rhs values
  liftedCut.splitVarIndex = cutinfo.splitVarIndex;
  liftedCut.splitIndex = cutinfo.splitIndex;
  liftedCut.RHS = cutinfo.RHS;
  liftedCut.coeff.insert(liftedCut.coeff.begin(), cutinfo.coeff.begin(),
      cutinfo.coeff.end());

  // For the remaining coefficients, choose max of ucoeff and vcoeff
  // Recall that some of these columns were complemented to form the subspace
  // (thus generating a different rhs from which cuts were found).
  // So for columns that were complemented, use -1 * the column to generate a valid lifting
  // into the complemented full space.
  // We will later make this valid for the original, non-complemented space.
  // Also, remember to use the dual multipliers for the bound constraints
  for (int var = 0; var < (int) deletedCols.size(); var++) {
    double ucoeff, vcoeff;
    int currvar = deletedCols[var];

#ifdef TRACE
    printf("Choosing coefficient for original variable %d (%s).\n", currvar,
        solver.getColName(currvar).c_str());
#endif
//    int currub = solver.getColUpper()[currvar];
//    int currlb = solver.getColLower()[currvar];
    // We do not use the bounds because when extending \bar u to the full space,
    // we put a multiplier of zero on those rows that did not exist in the subspace.

    int mult = varToDelete(currvar, &solver, subspace_option); // -1 if complemented, 1 if not, 0 if not deleted (shouldn't happen).
    if (std::abs(mult) < params.get(EPS)) {
      char errorstring[300];
      snprintf(errorstring, sizeof(errorstring) / sizeof(char),
          "Deleted var index: %d. The variable deleted: %d. Mult: %d.\n",
          var, currvar, mult);
      cerr << errorstring << endl;
      writeErrorToII(errorstring, inst_info_out);
      exit(1);
    }

    // Compute ucoeff: \bar u \cdot (a_j)^u
    ucoeff = 0;
    for (int i = 0; i < numSubRows; i++) {
//      printf("u corresponding to row %d (%s) is %f (mult %d). A coeff: %f.\n",
//          oldRowIndices[i], solver.getRowName(oldRowIndices[i]).c_str(),
//          DDual[i], mult, A->getCoefficient(oldRowIndices[i], currvar));
      ucoeff += DDual[i] * A->getCoefficient(oldRowIndices[i], currvar)
          * mult;
    }

    // Compute vcoeff: \bar v \cdot (a_j)^v
    vcoeff = 0;
    for (int i = 0; i < numSubRows; i++) {
//      printf("v corresponding to row %d (%s) is %f (mult %d). A coeff: %f.\n",
//          oldRowIndices[i], solver.getRowName(oldRowIndices[i]).c_str(),
//          DDual[numSubRows + 1 + i], mult,
//          A->getCoefficient(oldRowIndices[i], currvar));
      vcoeff += DDual[numSubRows + 1 + i]
          * A->getCoefficient(oldRowIndices[i], currvar) * mult;
    }

    // Choose max of the two
#ifdef TRACE
    printf("var: %d\tucoeff: %f\tvcoeff: %f\n", currvar, ucoeff, vcoeff);
#endif
    if (ucoeff > vcoeff) {
#ifdef TRACE
      printf("Choosing ucoeff with value %f for variable %d.\n", ucoeff,
          currvar);
#endif
      liftedCut.coeff.push_back(ucoeff);
    } else {
#ifdef TRACE
      printf("Choosing vcoeff with value %f for variable %d.\n", vcoeff,
          currvar);
#endif
      liftedCut.coeff.push_back(vcoeff);
    }
  }
}

/***********************************************************************/
/**
 * Reorders the lifted OsiCut so that the order of the columns is the same
 * as the order of the columns in the original problem. (This order has
 * changed as a result of creating the subspace and deleting variables.)
 */
void orderLiftedCut(OsiCut &liftedCut, SolverInterface *solver,
    SolverInterface *subprob, vector<int> &deletedCols,
    OsiCut &orderedLiftedCut) {
#ifdef TRACE
  std::cout
  << "\n############### Starting orderLiftedCut routine to generate ordered lifted coefficients. ###############"
  << std::endl;
#endif
  orderedLiftedCut.coeff.reserve(solver->getNumCols());
  orderedLiftedCut.splitVarIndex = liftedCut.splitVarIndex;
  orderedLiftedCut.splitIndex = liftedCut.splitIndex;
  orderedLiftedCut.RHS = liftedCut.RHS;

  int k = 0;
  int currdel = -1;
  int numSScols = 0;
  int totalSScols = subprob->getNumCols();
  for (int i = 0; i < solver->getNumCols(); i++) {
    // What's the next deleted variable?
    if ((int) deletedCols.size() >= 1 && currdel < i
        && k < (int) deletedCols.size()) {
      currdel = deletedCols[k];
      k++;
    }
    // If we have reached a deleted column, pull it from latter part of liftedCut
    if (currdel == i) {
      int index = i - numSScols + totalSScols; // index in liftedCut of currdel
      orderedLiftedCut.coeff.push_back(liftedCut.coeff[index]);
    } else { // Else this is a subspace variable
      orderedLiftedCut.coeff.push_back(liftedCut.coeff[numSScols]);
      numSScols++;
    }
  }
#ifdef TRACE
  std::cout
  << "\n## Lifted OsiCut has been ordered according to original columns. ##"
  << std::endl;
#endif
}

/**
 * @param inStatus  ::  1 if currently at ub, -1 if currently at lb
 */
void getPivots(SolverInterface* solver, SolutionInfo &solnInfo,
    vector<int> &colIn, vector<int> &colInNBIndex, vector<int> &inStatus,
    vector<int> &colOut, vector<int> &outStatus, vector<double> &minRatio,
    int &num_pivots, FILE* inst_info_out) {
  // Find a variable to pivot in (nb var with rc = 0)
#ifdef TRACE
  printf("Find variable to pivot in.\n");
#endif
  for (int i = 0; i < (int) solnInfo.nonBasicOrigVarIndex.size(); i++) {
    int col = solnInfo.nonBasicOrigVarIndex[i];
    double rc = solnInfo.reducedCost[col];
    if (rc == 0) {
      colIn.push_back(col);
      colInNBIndex.push_back(i);

      if (solnInfo.primalSoln[col] == solver.getColLower()[col]) {
        inStatus.push_back(-1);
      } else {
        inStatus.push_back(1);
      }
#ifdef TRACE
      printf("Found column to pivot in: %d.\n", col);
      printf("Now looking for variable to pivot out.\n");
#endif
      double minR = solver.getInfinity();
      int outcol = -1;
      int outstat = 0;
      //      printf("The variable to pivot in is at its lower-bound, so rays are in their normal directions.\n");
      vector<double> ray = solnInfo.raysOfC1[i];

      for (int row = 0; row < solnInfo.numRows; row++) {
        int var = solnInfo.varBasicInRow[row];
        double rayDirn = ray[row];
        double currRatio = solver.getInfinity();
        double val, lb, ub;
        int tmpstat = 0;
        if (var < solnInfo.numCols) {
          val = solnInfo.primalSoln[var];
          lb = solver.getColLower()[var];
          ub = solver.getColUpper()[var];
        } else {
          int currRow = var - solnInfo.numCols;
          val = solnInfo.slackSoln[currRow];
          if (solver.getObjSense() == 1) {
            if (solver.getRowSense()[currRow] == 'G') {
              lb = 0;
              ub = solver.getInfinity();
            } else if (solver.getRowSense()[currRow] == 'L') {
              lb = -1 * solver.getInfinity();
              ub = 0;
            } else {
              lb = -1 * solver.getInfinity();
              ub = solver.getInfinity();
            }
          } else {
            if (solver.getRowSense()[currRow] == 'L') {
              lb = 0;
              ub = solver.getInfinity();
            } else if (solver.getRowSense()[currRow] == 'G') {
              lb = -1 * solver.getInfinity();
              ub = 0;
            } else {
              lb = -1 * solver.getInfinity();
              ub = solver.getInfinity();
            }
          }
        }

//        printf(
//            "Row %d\n\tVar: %d\n\tvalue %f\n\tlb: %f\n\tub: %f\n\trayDirn: %f.\n",
//            row, var, val, lb, ub, rayDirn);

        // If rayDirn > 0, then we are increasing the variable in that row
        if (rayDirn > params.get(EPS)) {
          // Only provides any bound if there is an upper bound on the var
          // basic in this row.
          if (ub < solver.getInfinity() - params.get(EPS)) {
            currRatio = (ub - val) / rayDirn;
            tmpstat = 1;
//            printf("Ray is positive, ub exists, and currRatio is %f.\n",
//                currRatio);
          }
        } else if (rayDirn < -1 * params.get(EPS)) {
          // Only provides any bound if there is a lower bound on the var
          // basic in this row.
          if (lb > -1 * solver.getInfinity() + params.get(EPS)) {
            currRatio = (lb - val) / rayDirn;
            tmpstat = -1;
//            printf("Ray is negative, lb exists, and currRatio is %f.\n",
//                currRatio);
          }
        }

//        printf("CurrRatio: %f.\n", currRatio);

        if (currRatio < minR) {
          minR = currRatio;
          outcol = var;
          outstat = tmpstat;
        }

//        printf("\n");
      }

      colOut.push_back(outcol);
      outStatus.push_back(outstat);
      minRatio.push_back(minR);

      if (outcol == -1) {
#ifdef TRACE
        printf(
            "Could not find column to pivot out for alternate optimal solution. Problem is unbounded in this direction.\n"
            "This does not mean the problem is unbounded. Since rc on this variable is zero, we are not improving the objective anyway.\n"
            "We can simply take 1 or 1000 units in the direction of this ray, and that would give us an alternate optimal solution.\n"
            "We will do this at some point, maybe.\n\n");
#endif
        // Actually an extreme ray here may correspond to a really weak lifted inequality, because if the lifted coefficients were
        // \gamma_j before, now they are \gamma_j + \theta (r^j)^\T a^j
        // But we can keep track of how often this happens, and do a lifting with some large theta, to compare.
      } else {
#ifdef TRACE
        printf(
            "Found column to pivot out: %d. Amount to increase var %d by: %f. Outstatus: %d.\n\n",
            outcol, col, minR, outstat);
#endif
        if (std::abs(minR) > PIVOTEPS) { // TODO: Change this to whatever is reasonable for a minimum ratio
          num_pivots++;
        }
      }
    }
  }
}

/**
 * @brief Puts DDual ordered according to rows of the original subprob correctly.
 *
 * Suppose that the subproblem has m rows and n columns
 *  Constraints 0 through m-1 correspond to Ax^0 \ge b \lambda_0
 *  Constraint m is -x_k^0 \ge -\pi_0 \lambda_0
 *  Constraints m+1 through 2m correspond to Ax^1 \ge b \lambda_1
 *  Constraint 2m+1 is x_k^1 \ge (\pi_0 + 1) \lambda_1
 *  Constraint 2m+2 is \lambda_0 + \lambda_1 = 1
 * Then we go into lower-bounds and upper-bounds, but maybe not keep those?
 */
void orderDualVarsFromCplexOutput(vector<double> &DDual,
    SolverInterface* subprob, SolverInterface* dualsolver,
    SolverInterface* Pk) {
  const double* dualSoln = dualsolver.getColSolution();

  DDual.clear();
//  DDual.resize(2*(subprob.getNumRows()+1), -1);
  DDual.resize(dualsolver.getNumCols(), -1.0);
  for (int col = 0; col < dualsolver.getNumCols(); col++) {
    string colname = dualsolver.getColName(col);
    // Extract from colname the number of the row this corresponds to
    int row = atoi(colname.substr(4, colname.length() - 4).c_str());
//    if (row < 2*subprob.getNumRows()+2) {
//      dualVarIndex[row] = col;
//    }
//    printf("\nDual variable %d corresponds to row %d.\n", col, row);

    DDual[row] = dualSoln[col];

    // If this row was a <= row in the subproblem,
    // it would normally correspond to a negative slack
    // but CPLEX wants non-negative slacks, so instead
    // it negates coefficients. To get the true value
    // we should use for lifting for such rows, we need to
    // negate the value. This is all assuming a minimization problem,
    // which Pk should be (min \bar \alpha^\T x).
    if (Pk.getRowSense()[row] == 'L')
      DDual[row] *= -1;
  }
}

/**
 * @brief Puts DDual ordered according to rows of the original subprob correctly.
 *
 * Suppose that the subproblem has m rows and n columns
 *  Constraints 0 through m-1 correspond to Ax^0 \ge b \lambda_0
 *  Constraint m is -x_k^0 \ge -\pi_0 \lambda_0
 *  Constraints m+1 through 2m correspond to Ax^1 \ge b \lambda_1
 *  Constraint 2m+1 is x_k^1 \ge (\pi_0 + 1) \lambda_1
 *  Constraint 2m+2 is \lambda_0 + \lambda_1 = 1
 * Then we go into lower-bounds and upper-bounds, but maybe not keep those?
 */
void orderDualVars(vector<double> &DDual, SolverInterface* subprob,
    SolverInterface* dualsolver, SolverInterface* Pk) {
  const double* dualSoln = dualsolver.getColSolution();

  DDual.clear();
//  DDual.resize(2*(subprob.getNumRows()+1), -1);
  DDual.resize(dualsolver.getNumCols(), -1.0);
  for (int col = 0; col < dualsolver.getNumCols(); col++) {
    string colname = dualsolver.getColName(col);
    // Extract from colname the number of the row this corresponds to
    int row = atoi(colname.substr(4, colname.length() - 4).c_str());
//    if (row < 2*subprob.getNumRows()+2) {
//      dualVarIndex[row] = col;
//    }
//    printf("\nDual variable %d corresponds to row %d.\n", col, row);

    DDual[row] = dualSoln[col];

    // If this row was a <= row in the subproblem,
    // it would normally correspond to a negative slack
    // but CPLEX wants non-negative slacks, so instead
    // it negates coefficients. To get the true value
    // we should use for lifting for such rows, we need to
    // negate the value. This is all assuming a minimization problem,
    // which Pk should be (min \bar \alpha^\T x).
    if (Pk.getRowSense()[row] == 'L')
      DDual[row] *= -1;
  }
} /* orderDualVars */
