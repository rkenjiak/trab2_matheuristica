/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*      This file is based on other part of the program and library          */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  It is an example of a branch-and-bound code using SCIP as library 
 *
 * @author Edna Hoshino (based on template codified by Timo Berthold and Stefan Heinz)
 *
 * This is an example for solving the knapsack problem. 
 * The goal of this problem is finding a subset of items from a input set,
 * to maximize the sum of the values of selected items subject to their weights do not exceed a given capacity C.
 * We also use a naive primal heuristic to generate a feasible solution quickly.
 * 
 **/ 

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include<stdio.h>
#include<time.h>
#include<string.h>

#include "scip/scip.h"
#include "problem.h"
#include "utils.h"

int main(int argc, char **argv)
{
  SCIP* scip;
  instanceT* in;
  clock_t start, end;
  char outputname[SCIP_MAXSTRLEN];
  parametersT param;

  // set default+user parameters
  if(!setParameters(argc, argv, &param))
     return 0;

  // load instance file
  if(!loadInstance(argv[1], &in)){
    printf("\nProblem to read instance file %s\n", argv[1]);
    return 1;
  }
#ifdef DEBUG
  printInstance(in);
#endif
  // create scip and set scip configurations
  configScip(&scip, &param);
  // load problem into scip
  if(!loadProblem(scip,argv[1],in,0,NULL,&param)){
    printf("\nProblem to load instance problem\n");
    return 1;
  }
  // print problem
  SCIP_CALL( SCIPwriteOrigProblem(scip, "knapsack.lp", "lp", FALSE) );
  srand(time(NULL));
  // solve scip problem
  start=clock();
  SCIP_CALL( SCIPsolve(scip) );
  end = clock();
  // config output filename
  configOutputName(outputname, argv[1], argv[0], &param);
  // print statistics and print resume in output file
  printStatistic(scip,((double) (end-start))/CLOCKS_PER_SEC, outputname);
  // write the best solution in a file
  printSol(scip, outputname);
  SCIP_CALL( SCIPfree(&scip) ); 
  freeInstance(in);
  BMScheckEmptyMemory();
  return 0;
}
