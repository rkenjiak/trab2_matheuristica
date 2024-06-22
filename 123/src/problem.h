#ifndef __PROBLEM__
#define __PROBLEM__
#include<stdio.h>
#include "scip/scip.h"
#include "parameters_mochila.h"

/** structure for each item */
typedef struct{
  int label;
  int value;
  int weight;
  int nsets; /**< total of forfeit set where item belongs to */
  int *set;  /**< forfeit sets that contain the item */
}itemType;

/** structure for each forfeit set */
typedef struct{
   int j; /**< label for the forfeit set */
   int h; /**< maximum of items without to pay the forfeit cost */
   int d; /**< forfeit cost */
   int n; /**< total of items in the forfeit set */
   int *items; /**< list of items in the forfeit set */
}forfeitType;

typedef struct{
  int n;  /**< total of items */
  int nS; /**< total of forfeit sets */
  int C;  /**< knapsack capacity */
  int k;  /**< forfeit limits */
  forfeitType *S; /**< data for each forfeit set */
  itemType *item; /**< data for each item in 0..n-1 */
} instanceT;

void freeInstance(instanceT* I);
void createInstance(instanceT** I, int n, int nS, int C);
void printInstance(instanceT* I);
// load instance from a file
int loadInstance(char* filename, instanceT** I);
// load instance problem into SCIP
int loadProblem(SCIP* scip, char* probname, instanceT* in, int relaxed, int* fixed, parametersT* param);
#endif
