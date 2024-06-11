/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_aleatoria.c
 * @brief  aleatoria primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_aleatoria.h"
#include "utils.h"
//#include "heur_problem.h"
/* configuracao da heuristica */
#define HEUR_NAME             "aleatoria"
#define HEUR_DESC             "aleatoria primal heuristic template"
#define HEUR_DISPCHAR         'a'
#define HEUR_PRIORITY         2 /**< heuristics of high priorities are called first */
#define HEUR_FREQ             1 /**< heuristic call frequency. 1 = in all levels of the B&B tree */
#define HEUR_FREQOFS          0 /**< starts of level 0 (root node) */
#define HEUR_MAXDEPTH         -1 /**< maximal level to be called. -1 = no limit */
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE /**< when the heuristic should be called? SCIP_HEURTIMING_DURINGLPLOOP or SCIP_HEURTIMING_AFTERNODE */
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#ifdef DEBUG
   #define PRINTF(...) printf(__VA_ARGS__)
#else
   #define PRINTF(...) 
#endif

/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
/*struct SCIP_HeurData
{
};
*/

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitAleatoria)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/**
 * @brief Core of the aleatoria heuristic: it builds one solution for the problem by aleatoria procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure where the solution wil be saved
 * @param heur pointer to the aleatoria heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int aleatoria(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur)
{
   int found, infeasible, nInSolution;
   unsigned int stored;
   int nvars;
   int *covered, n, custo, nCovered, *cand, nCands, selected, s, *forfeit, valor, violations;
   SCIP_VAR *var, **solution, **varlist;
   //  SCIP* scip_cp;
   SCIP_Real bestUb;
   SCIP_PROBDATA* probdata;
   int i, residual, j, peso, nS, ii, toBeViolated;
   instanceT* I;
   
   found = 0;
   infeasible = 0;
   
#ifdef DEBUG_ALEATORIA
   printf("\n============== New aleatoria heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recover the problem data*/
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   nvars = SCIPprobdataGetNVars(probdata);
   varlist = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
   nS = I->nS; // nS = total de forfeit sets
    
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*(n+nS));
   covered = (int*) calloc(n,sizeof(int)); // inicializa que nenhum item estah na solucao
   forfeit = (int*) calloc(nS,sizeof(int)); // inicializa que nenhum forfeit estah coberto na solucao
   cand = (int*) malloc(n*sizeof(int));
   nInSolution = 0;
   nCovered = 0;
   nCands = 0;
   custo = 0;
   residual = I->C;
   violations = 0;

   // first, select all variables already fixed in 1.0
   for(i=0;i<nvars;i++){
      var = varlist[i];
      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){ // var >= 1.0
        if (i<n){
          solution[nInSolution++]=var;
          // update residual capacity
          residual -= I->item[i].weight;
          covered[i]=1;
          // for each forfeit set that item i belongs to ...
          for(j=0;j<I->item[i].nsets;j++){
             ii = I->item[i].set[j];
             forfeit[ii]++; // update total of items selected from the forfeit set
          }
          // update solution value
          custo += I->item[i].value;
          infeasible = residual < 0?1:0;
#ifdef DEBUG_ALEATORIA
          printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, residual, infeasible);
#endif
        }
        else{
           // refers to variable v over forfeit set
           j = n - i;
           // nothing is necessary to do... 
        }
      }
      else{ // discard items fixed in 0.0
        if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
        }
        else{
          if (i < n){ // include item i in cand list
           cand[nCands++]=i;
          }
        }
      }
   }
   // update forfeits violated
   for(i=0;i<nS;i++){
      if(forfeit[i] > I->S[i].h){
         violations += forfeit[i] - I->S[i].h;
      }
   }
   // complete solution using items not fixed (not covered)
   for(i=0;i<n && !infeasible && nCands > 0 && residual>0;i++){
      s = RandomInteger (0, nCands-1);
      selected = cand[s]; // selected candidate
      cand[s] = cand[--nCands]; // remove selected candidate
      // only accept the item if not covered yet and not exceed the capacity
      if(!covered[selected] && I->item[selected].weight <= residual){
         // compute the real value
         toBeViolated = 0;
         valor = I->item[selected].value;
         for(j=0;j<I->item[selected].nsets;j++){
            ii = I->item[selected].set[j];
            // update the value if the item will exceed the maximum allowed for the set
            if(forfeit[ii] >= I->S[ii].h){
               valor -= I->S[ii].d;
               toBeViolated++;
            }
         }
         // if it worths
         if(valor>0 && toBeViolated+violations <= I->k){
            var = varlist[selected];
            // include selected var in the solution
            solution[nInSolution++]=var;
            // update residual capacity
            residual -= I->item[selected].weight;
            // update covered
            covered[selected] = 1;
            // update the solution value
            custo += valor;
            // update the total of elements in each set that are already in the solution
            for(j=0;j<I->item[selected].nsets;j++){
               ii = I->item[selected].set[j];
               forfeit[ii]++;
            }
            infeasible = residual<0?1:0;
            violations += toBeViolated;
#ifdef DEBUG_ALEATORIA
            printf("\n\nSelected var= %s. TotalItems=%d value item=%d toBeViolated=%d value = %d residual=%d violations=%d infeasible=%d\n", SCIPvarGetName(var), nInSolution, valor, toBeViolated, custo, residual, violations, infeasible);
#endif
         }
         else{
            // desconsidere o item
#ifdef DEBUG_ALEATORIA
            printf("\n\nNOT selected var= %s. TotalItems=%d value item=%d value=%d residual=%d violations=%d infeasible=%d\n", SCIPvarGetName(varlist[selected]), nInSolution, valor, custo, residual, violations, infeasible);
#endif
         }
      }
   }
   if(!infeasible){
      /* create SCIP solution structure sol */
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
      // save found solution in sol
      for(i=0;i<nInSolution;i++){
         var = solution[i];
         SCIP_CALL( SCIPsetSolVal(scip, *sol, var, 1.0) );
      }
      // update forfeit set variable
      for(j=0;j<nS;j++){
         valor = forfeit[j]>I->S[j].h?forfeit[j]-I->S[j].h:0;
         SCIP_CALL( SCIPsetSolVal(scip, *sol, varlist[I->n+j], (double) valor) );
      }
      bestUb = SCIPgetPrimalbound(scip);
#ifdef DEBUG_ALEATORIA
      printf("\nFound solution...\n");
      printf("\ninfeasible=%d value = %d > bestUb = %lf? %d\n\n", infeasible, custo, bestUb, custo > bestUb + EPSILON);
#endif
      if(!infeasible && custo > bestUb + EPSILON){
#ifdef DEBUG_ALEATORIA
         printf("\nBest solution found...\n");
         SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
         
         /* check if the solution is feasible */
         SCIP_CALL( SCIPtrySolMine(scip, *sol, TRUE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
#ifdef DEBUG_PRIMAL
            printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
            SCIPdebugMessage("found feasible aleatoria solution:\n");
            SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
            found = 1;
         }
         else{
            found = 0;
#ifdef DEBUG_ALEATORIA
            printf("\nCould not found\n. BestUb=%lf", bestUb);
#endif
         }
      }
   }
   free(cand);
   free(solution);
   free(forfeit);
   free(covered);
   return found;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAleatoria)
{  /*lint --e{715}*/
   SCIP_SOL*             sol;                /**< solution to round */
   int nlpcands;

   assert(result != NULL);
   //   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* continue only if the LP is finished */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* continue only if the LP value is less than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


   /* check if there exists integer variables with fractionary values in the LP */
   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL, NULL) );
   //Fractional implicit integer variables are stored at the positions *nlpcands to *nlpcands + *nfrac - 1
  
   /* stop if the LP solution is already integer   */
   if ( nlpcands == 0 )
     return SCIP_OKAY;

   /* solve aleatoria */
   if(aleatoria(scip, &sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
     *result = SCIP_DIDNOTFIND;
#ifdef DEBUG_PRIMAL
     printf("\nAleatoria could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the aleatoria_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurAleatoria(
   SCIP*                 scip,                /**< SCIP data structure */
   const parametersT* param
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;
   /* create aleatoria primal heuristic data */
   heurdata = NULL;
   heur = NULL;
   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param->heur_freq, param->heur_freqofs,
         param->heur_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyAleatoria, heurFreeAleatoria, heurInitAleatoria, heurExitAleatoria, heurInitsolAleatoria, heurExitsolAleatoria, heurExecAleatoria,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param->heur_freq, param->heur_freqofs,
         param->heur_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAleatoria, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAleatoria) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAleatoria) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAleatoria) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitAleatoria) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolAleatoria) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolAleatoria) );
#endif

   /* add aleatoria primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
