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

/**@file   heur_lns.c
 * @brief  lns primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "utils.h"
#include "heur_lns.h"

//#define DEBUG_LNS 1
/* configuracao da heuristica */
#define HEUR_NAME             "lns"
#define HEUR_DESC             "lns primal heuristic template"
#define HEUR_DISPCHAR         'l'
#define HEUR_PRIORITY         1 /**< heuristics of high priorities are called first */
#define HEUR_FREQ             1 /**< heuristic call frequency. 1 = in all levels of the B&B tree */
#define HEUR_FREQOFS          0 /**< starts of level 0 (root node) */
#define HEUR_MAXDEPTH         -1 /**< maximal level to be called. -1 = no limit */
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE /**< when the heuristic should be called? SCIP_HEURTIMING_DURINGLPLOOP or SCIP_HEURTIMING_AFTERNODE */
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

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
/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLns)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLns)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitLns)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLns)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolLns)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolLns)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/**
 * @brief Core of the lns heuristic: it builds one solution for the problem by lns procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure to be improved
 * @param heur pointer to the lns heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int lns(SCIP* scip, SCIP_SOL* initsol, SCIP_HEUR* heur)
{
   parametersT lnsparam;
   const parametersT* param;
   SCIP* subscip;
   int found, nInSolution;
   unsigned int stored;
   int custo, nFixed;
   SCIP_VAR *var, **vars, **vars2;
   SCIP_Real valor;
   SCIP_PROBDATA* probdata, *probdata2;
   SCIP_SOL *lnsSol, *sol;
   int i, ii, tobeViolated, violations;
   instanceT* I;
   double z, lnsZ, bestUb;
   int *fixed;
   itemType *cand;
   int nCands, capacRes, toRemove, perda, nRemoved;
#ifdef DEBUG_LNS
   int infeasible;
   unsigned int status;
#endif
   found = 0;
#ifdef DEBUG_LNS
   infeasible = 0;
   printf("\n============== New lns heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recover the problem data from original problem */
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);
   vars = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   param = SCIPprobdataGetParam(probdata);
   nFixed = 0;
   custo = 0;
   capacRes = I->C;
   // aloca candidatos
   cand = (itemType*)malloc(sizeof(itemType)*I->n);
   fixed = (int*)calloc(I->n, sizeof(int)); // fixed[i]=0, if item i is not fixed, fixed[i]=1 if item i is fixed in 1.0, fixed[i]=-1 if item i is fixed in 0.

   // first, select all variables already fixed in 1.0
   for(i=0;i<I->n;i++){
      var = vars[i];
      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){ // var >= 1.0
        // update residual capacity
        capacRes -= I->item[i].weight;
        fixed[i] = 1;
        nFixed++;
        custo += I->item[i].value;
#ifdef DEBUG_LNS
        infeasible = capacRes < 0?1:0;
        printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, capacRes, infeasible);
#endif
      }
      else{ // discard items fixed in 0.0
        if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
          fixed[i] = -1;
          nFixed++;
        }
      }
   }
   // Constroi lista de candidatos a remover da solucao atual e atualiza capacidade residual
   nCands = 0;
   for(i=0;i<I->n;i++){
     valor = SCIPgetSolVal(scip, initsol, vars[i]); // checa se o item esta na solucao atual
     if(valor > EPSILON && !fixed[i]){
       fixed[i]=1;
       cand[nCands].label = i;
       cand[nCands++].value = -(I->item[i].weight);
       capacRes -= I->item[i].weight;
    }
  }
  
  // Decide quem sairá da solução com base no peso de cada item
  qsort(cand, nCands, sizeof(itemType), comparador); // Ordena o vetor com base no peso dos itens selecionados
  toRemove= nCands*(param->lns_perc); // calcula total a ser destruida
#ifdef DEBUG
  printf("\nporcDestroy=%lf total a destruir=%d total de candidatos=%d\n", param->lns_perc, toRemove, nCands);
#endif
  perda = 0; 
  nRemoved = 0;
  for(i=0;i<toRemove;i++){
    ii = cand[i].label;
    capacRes += I->item[ii].weight;
    perda += I->item[ii].value;
#ifdef DEBUG
    printf("\nRemove %d (peso=%d valor=%d) da mochila (capac residual=%d)", I->item[ii].label, I->item[ii].weight, I->item[ii].value, capacRes);
#endif
    fixed[ii]=0;
    nRemoved++;
  }

#ifdef DEBUG
  printf("\nDestrui %.2lf%% (equivalente a um valor = %d)\n", (100.0*nRemoved/nCands), perda);
#endif

  // create scip and set scip configurations
  lnsparam.time_limit = param->lns_time;
  lnsparam.display_freq = -1;
  lnsparam.nodes_limit = -1;
  lnsparam.heur_rf = 0;
  lnsparam.heur_lns = 0;
  lnsparam.heur_aleatoria = 0;
  configScip(&subscip, &lnsparam);
  /* disable output to console */
  SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
  SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", 1) ); 
  
  // carga do lp
  // load problem into scip
  if(!loadProblem(subscip, "lns", I, 0, fixed, &lnsparam)){
     printf("\nProblem to load instance problem\n");
     return -1;
  }
  probdata2 = SCIPgetProbData(subscip);
  assert(probdata2 != NULL);
  // Recupera vars
  vars2 = SCIPprobdataGetVars(probdata2);
  // print problem
  SCIP_CALL( SCIPwriteOrigProblem(subscip, "lns.lp", "lp", FALSE) );
     
  // solve scip problem
  SCIP_CALL( SCIPsolve(subscip) );
#ifdef DEBUG_LNS
  SCIP_CALL( SCIPprintBestSol(subscip, NULL, FALSE) );
#endif
  
#ifdef DEBUG_LNS
  status = SCIPgetStatus(subscip);
  printf("\nstatus=%d", status);
#endif
     
  // Recupera solucao
  lnsSol = SCIPgetBestSol(subscip);
  lnsZ = SCIPgetPrimalbound(subscip); // o.f. for the solution found by LNS
  z = SCIPsolGetOrigObj(initsol); // objective function for the initial solution 
  if (lnsZ > z + EPSILON){
#ifdef DEBUG
     printf("\nSolucao do LNS:");
#endif
     /* create SCIP solution structure sol */
     SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
     nInSolution = 0;
     for (i = 0; i < I->n; i++)
     {
        valor = SCIPgetSolVal(subscip, lnsSol, vars2[i]);
        if(valor>EPSILON){
          // set found solution in sol (for original problem)
          SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], 1.0) );
#ifdef DEBUG
          printf("\nItem %d (peso=%d valor=%d)", I->item[i].label, I->item[i].weight, I->item[i].value);
#endif
          nInSolution++;
        }
     }
     // set o valor das variaveis de forfeit set
     for(i=0;i<I->nS;i++){
        valor = SCIPgetSolVal(subscip, lnsSol, vars2[I->n+i]); // recupera o valor da variavel yi
        if(valor>EPSILON){
           SCIP_CALL( SCIPsetSolVal(scip, sol, vars[I->n+i], valor) );
        }
     }
     // check if the solution found by LNS is better than the current bestsolution
     bestUb = SCIPgetPrimalbound(scip);
#ifdef DEBUG_LNS
     printf("\nFound solution...\n");
     //      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
     printf("\ninfeasible=%d value = %lf > bestUb = %lf? %d\n\n", infeasible, lnsZ, bestUb, lnsZ > bestUb + EPSILON);
#endif
     if(lnsZ > bestUb + EPSILON){
#ifdef DEBUG
      printf("\nBest solution found...\n");
      SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
      /* check if the solution is feasible */
      SCIP_CALL( SCIPtrySolMine(scip, sol, TRUE, TRUE, FALSE, TRUE, &stored) );
      if( stored )
      {
#ifdef DEBUG_PRIMAL
        printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
        SCIPdebugMessage("found feasible lns solution:\n");
        SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
      }
      else{
        found = 0;
#ifdef DEBUG_LNS
        printf("\nCould not found better solution\n. BestUb=%lf", bestUb);
#endif
      }
    }
  }
  // Destroi problema
  free(cand);
  free(fixed);
  // clear problem
  SCIP_CALL( SCIPfree(&subscip) );
  return found;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLns)
{  /*lint --e{715}*/
   SCIP_SOL*             sol;                /**< solution to improve */
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

   sol = SCIPgetBestSol(scip);
   if(sol==NULL)
      return SCIP_OKAY;
   /* solve lns */
   if(lns(scip, sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
     *result = SCIP_DIDNOTFIND;
#ifdef DEBUG_PRIMAL
     printf("\nLns could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific intelnsace methods
 */

/** creates the lns_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLns(
   SCIP*                 scip,                /**< SCIP data structure */
   const parametersT* param
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create lns primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyLns, heulnsreeLns, heurInitLns, heurExitLns, heurInitsolLns, heurExitsolLns, heurExecLns,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLns, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLns) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLns) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLns) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolLns) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolLns) );
#endif

   /* add lns primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
