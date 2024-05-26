/**@file   problem.c
 * @brief  This file contains routines specific for the problem and the functions loadInstance(), freeInstance, 
 * printInstance, and loadProblem must be implemented 
 *
 **/ 
#include<stdio.h>
#include<math.h>
#include "scip/scip.h"
#include "problem.h"
#include "probdata_mochila.h"

void freeInstance(instanceT* I)
{
   int i,j;
  
   if(I){
      for(i=0;i<I->n;i++){
         free(I->item[i].set);
      }
      free(I->item);
      for(j=0;j<I->nS;j++)
         free((I->S[j]).items);
      free(I->S);
      free(I);
      I = NULL;
   }
}
void createInstance(instanceT** I, int n, int nS, int C)
{
  int i;
  
  *I = (instanceT*) malloc(sizeof(instanceT));
  (*I)->item = (itemType*) malloc(sizeof(itemType)*n);
  (*I)->S = (forfeitType*) malloc(sizeof(forfeitType)*nS);
  (*I)->n = n;
  (*I)->nS = nS;
  (*I)->C = C;
}
void printInstance(instanceT* I)
{
  int i,j;
  printf("\nInstance with n=%d items, C=%d, nS=%d forfeits, and k=%d", I->n, I->C, I->nS, I->k);
  for(j=0;j<I->nS;j++){
     printf("\nForfeit %d h=%d d=%d = {", I->S[j].j, I->S[j].h, I->S[j].d);
     for(i=0;i<I->S[j].n;i++){
        printf("%d ", I->S[j].items[i]);
     }
     printf("\n");
  }
  printf("\nItems= \n");
  for(i=0;i<I->n;i++){
     printf("%d value=%d weight=%d sets={", I->item[i].label, I->item[i].value, I->item[i].weight);
     for(j=0;j<I->item[i].nsets;j++){
        printf("%d ", I->item[i].set[j]);
     }
     printf("}\n");
  }
}
int loadInstance(char* filename, instanceT** I)
{
  FILE* fin;
  int n, nS, i, j, C, ii;
  fin = fopen(filename, "r");
  if(!fin){
    printf("\nProblem to open file %s\n", filename);
    return 0;
  }

  fscanf(fin,"%d %d %d\n", &n, &nS, &C);
  createInstance(I, n, nS, C);
  for(i=0; i<n; i++){
     fscanf(fin, "%d\n", &((*I)->item[i].value));
     (*I)->item[i].label=i;
     (*I)->item[i].nsets=0;
  }
  for(i=0; i<n; i++){
     fscanf(fin, "%d\n", &((*I)->item[i].weight));
  }
  for(j=0; j<nS; j++){
     (*I)->S[j].j = j;
     fscanf(fin,"%d %d %d\n", &((*I)->S[j].h), &((*I)->S[j].d), &((*I)->S[j].n));
     (*I)->S[j].items = (int*)malloc(sizeof(int)*((*I)->S[j].n));
     for(i=0; i<(*I)->S[j].n; i++){
        fscanf(fin, "%d", &ii);
        ((*I)->S[j].items[i]) = ii;
        ((*I)->item[ii].nsets)++;
     }
  }
  for(i=0; i<n; i++){
     (*I)->item[i].set = (int*)malloc(sizeof(int)*((*I)->item[i].nsets));
     (*I)->item[i].nsets = 0;
  }
  for(j=0; j<nS; j++){
     for(i=0; i<(*I)->S[j].n; i++){
        ii = (*I)->S[j].items[i];
        (*I)->item[ii].set[((*I)->item[ii].nsets)++] = j;
     }
  }
  fscanf(fin, "\nk %d", &((*I)->k));
  fclose(fin);
  return 1;
}
// load instance problem into SCIP
int loadProblem(SCIP* scip, char* probname, instanceT* I, int relaxed, int* fixed, parametersT* param)
{
  SCIP_RETCODE ret_code;

  ret_code = SCIPprobdataCreate(scip, probname, I, relaxed, fixed, param);
  if(ret_code!=SCIP_OKAY)
    return 0;
  return 1;
}
