#ifndef __UTILS_H__
#define __UTILS_H__
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "problem.h"
#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_aleatoria.h"
#include "heur_rf.h"
#include "heur_lns.h"
SCIP_RETCODE SCIPtrySolMine(SCIP* scip, SCIP_SOL* sol, SCIP_Bool printreason, SCIP_Bool checkbounds, SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool *stored);
void removePath(char* fullfilename, char** filename);
void configOutputName(char* name, char* instance_filename, char* program, const parametersT* param);
SCIP_RETCODE printStatistic(SCIP* scip, double time, char* outputname);
void printSol(SCIP* scip, char* outputname);
SCIP_RETCODE configScip(SCIP** pscip, const parametersT* param);
//
/* sorteia um numero aleatorio entre [low,high] */
int RandomInteger(int low, int high);
/* put your local methods here, and declare them static */
SCIP_RETCODE SCIPtrySolMine(SCIP* scip, SCIP_SOL* sol, SCIP_Bool printreason, SCIP_Bool checkbounds, SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool *stored);
int comparador(const void *valor1, const void *valor2);
#endif
