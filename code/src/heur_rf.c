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

/**@file   heur_rf.c
 * @brief  rf primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "utils.h"
#include "heur_rf.h"

//#define DEBUG_RF 1
/* configuracao da heuristica */
#define HEUR_NAME             "rf"
#define HEUR_DESC             "rf primal heuristic template"
#define HEUR_DISPCHAR         'r'
#define HEUR_PRIORITY         2 /**< heuristics of high priorities are called first */
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
SCIP_DECL_HEURCOPY(heurCopyRf)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRf)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRf)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRf)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRf)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRf)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/**
 * @brief Core of the rf heuristic: it builds one solution for the problem by rf procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure where the solution wil be saved
 * @param heur pointer to the rf heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int rf(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur)
{
   parametersT rfparam;
   const parametersT* param;
   SCIP* subscip;
   int found, nInSolution;
   unsigned int stored, infeasible;
   int n, custo, nFixed;
   SCIP_VAR *var, **solution, **vars, **vars2;
   //  SCIP* scip_cp;
   SCIP_Real valor, bestUb;
   SCIP_PROBDATA* probdata, *probdata2;
   SCIP_SOL* bestSolution;
   int i;
   instanceT* I;
   double z;
   int k, *particao, K, parte, frac, *fixed, tam;
#ifdef DEBUG_RF
   char destaque='*', branco=' ';
   int status;
#endif
   itemType *cand;
   int nCand, capacRes;
#ifdef TESTE
   int ii;
#endif
   
   found = 0;
   infeasible = 0;
   z = 0;
#ifdef DEBUG_RF
   printf("\n============== New rf heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recover the problem data from original problem */
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);
   vars = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
   param = SCIPprobdataGetParam(probdata);
   
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*n);
   nInSolution = 0;
   nFixed = 0;
   custo = 0;
   capacRes = I->C;   
   // aloca candidatos
   cand = (itemType*)malloc(sizeof(itemType)*I->n);
   fixed = (int*)calloc(I->n, sizeof(int)); // fixed[i]=0, if item i is not fixed, fixed[i]=1 if item i is fixed in 1.0, fixed[i]=-1 if item i is fixed in 0.

   nCand = 0;
   // first, select all variables already fixed in 1.0
   for(i=0;i<I->n;i++){
      var = vars[i];
      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){ // var >= 1.0
        solution[nInSolution++]=var;        
        // update residual capacity
        capacRes -= I->item[i].weight;
        fixed[i] = 1;
        nFixed++;
        custo += I->item[i].value;
        infeasible = capacRes < 0?1:0;
#ifdef DEBUG_RF
        printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, capacRes, infeasible);
#endif
      }
      else{ // discard items fixed in 0.0
        if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
          fixed[i] = -1;
          nFixed++;
        }
        else{ // candidate to be selected
          cand[nCand++] = I->item[i];
        }
      }
   }

   //  qsort(cand, nCand, sizeof(itemType), comparador);
   
   // define tamanho de cada parte da particao
   tam = ceill(param->rf_perc * I->n);
   K = ceill(1.0/param->rf_perc); // total de partes
   // aloca vetor das particoes
   particao = (int*)calloc(I->n, sizeof(int));
   // particiona as variaveis
   for(parte=0;parte<K;parte++){
     for(i=0;i<tam && nCand>0;i++){
       k = RandomInteger(0,nCand-1); // sorteia um candidato
       particao[cand[k].label] = parte;
       // remove candidato
       cand[k] = cand[--nCand];
     }
   }
#ifdef DEBUG_RF2
   printf("\nTam de cada parte=%d. Total de partes=%d\n", tam, K);
   for(i=0;i<I->n;i++){
     printf("\ni=%d parte=%d (fixed=%d)", i, particao[i], fixed[i]);
   }
#endif
  // create scip and set scip configurations
  rfparam.time_limit = param->rf_time;
  rfparam.display_freq = -1;
  rfparam.nodes_limit = -1;
  rfparam.heur_rf = 0;
  rfparam.heur_lns = 0;
  rfparam.heur_aleatoria = 0;
  configScip(&subscip, &rfparam);
  /* disable output to console */
  SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
  frac = 1;
  for(parte=0;parte<K && frac && !infeasible && nFixed < I->n;parte++){ // itera para cada parte da particao, mas pode parar antes se a solucao ja eh inteira.
#ifdef DEBUG_RF
    printf("\n=====parte=%d K=%d\n", parte, K);
#endif
    // carga do lp
    // load problem into scip
    if(!loadProblem(subscip, "rf", I, 1, fixed, &rfparam)){ // relaxation
      printf("\nProblem to load instance problem\n");
      return 0;
    }
    /* recover the problem data from subproblem */
    probdata2=SCIPgetProbData(subscip);
    assert(probdata2 != NULL);    
    vars2 = SCIPprobdataGetVars(probdata2);

    for(i=0;i<I->n;i++){
      if(particao[i]==parte){ // torna as variaveis da parte atual como binarias
        //        printf("\nMuda para binario a var x[%d]", i);
        SCIPchgVarType(subscip, vars2[i], SCIP_VARTYPE_BINARY, &infeasible);
      }
    }
#ifdef DEBUG_RF
    // print problem
    printf("\n---LP gravado em relax_fix.lp");
    SCIP_CALL( SCIPwriteOrigProblem(subscip, "relaxFix.lp", "lp", FALSE) );
#endif
    // Executa Solver de PL
    SCIP_CALL( SCIPsolve(subscip) );
  
#ifdef DEBUG_RF
    status = SCIPgetStatus(subscip);
    printf("\nstatus=%d", status);
#endif
     
    // Recupera solucao
    bestSolution = SCIPgetBestSol(subscip);
    z = SCIPgetPrimalbound(subscip);//SCIPsolGetOrigObj(bestSolution);
    if(bestSolution==NULL){
       infeasible = 1;
       break; // para o RF se ficou inviavel
    }
#ifdef DEBUG_RF
    printf("\n*****\n\t z=%lf\n", z);
#endif
    frac = 0;
    for (i = 0; i < I->n; i++)
    {
       valor = SCIPgetSolVal(subscip, bestSolution, vars2[i]); // recupera o valor da variavel xi no subproblem
       if(particao[i]==parte){ // fixa variaveis da parte atual
#ifdef DEBUG_RF
         printf("x%2d = %.2lf (parte = %2d)%c%c\n", i+1, valor, particao[i], particao[i]<=parte?destaque:branco, particao[i]==parte?destaque:branco);
#endif
         if(valor > EPSILON){
           solution[nInSolution++]=vars[i]; // inclui a variavel do item no problema original
            fixed[i]=1; // fixa var xi em 1
            nFixed++;
            capacRes -= I->item[i].weight;
          }
          else{
            fixed[i]=-1; // fixa var xi em 0
            nFixed++;
          }
       }
       else if(valor>EPSILON){
#ifdef DEBUG_RF
          printf("x%2d = %.2lf (parte = %2d)%c%c\n", i+1, valor, particao[i], particao[i]<=parte?destaque:branco, particao[i]==parte?destaque:branco);
#endif
          if(valor < 1.0 - EPSILON){
             frac = 1;
          }
       }
    }
    //    getchar();
  } // next part

  if(!frac && !infeasible){
#ifdef DEBUG_RF
    printf("\nSol heur relax fix = %lf\n", z);
#endif
    if(parte<K){
      // falta copiar a solução das demais variaveis
       for (i = 0; i < I->n; i++)
       {
          valor = SCIPgetSolVal(subscip, bestSolution, vars2[i]); // recupera o valor da variavel xi
          if(valor>EPSILON && particao[i]>=parte){
            fixed[i] = 1;
            solution[nInSolution++]=vars[i];
#ifdef DEBUG_RF
             printf("x%2d = %.2lf (parte = %2d)%c%c\n", i+1, valor, particao[i], particao[i]<=parte?destaque:branco, particao[i]==parte?destaque:branco);
#endif
          }
       }
    }
  }  
  if(!infeasible){
    //    SCIP_CALL( SCIPprintSol(subscip, bestSolution, NULL, FALSE) );
    /* create SCIP solution structure sol */
    SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
    // save found solution in sol
    for(i=0;i<nInSolution;i++){
      var = solution[i];
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, 1.0) );
    }
    // set o valor das variaveis de forfeit set
    for(i=0;i<I->nS;i++){
      valor = SCIPgetSolVal(subscip, bestSolution, vars2[I->n+i]); // recupera o valor da variavel yi
      if(valor>EPSILON){
        SCIP_CALL( SCIPsetSolVal(scip, *sol, vars[I->n+i], valor) );
      }
    }
    bestUb = SCIPgetPrimalbound(scip);
#ifdef DEBUG_RF
    printf("\nFound solution...\n");
    //      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
    printf("\ninfeasible=%d value = %lf > bestUb = %lf? %d\n\n", infeasible, z, bestUb, z > bestUb + EPSILON);
#endif
    if(z > bestUb + EPSILON){
#ifdef DEBUG_RF
      printf("\nBest solution found...\n");
      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
      /* check if the solution is feasible */
      SCIP_CALL( SCIPtrySolMine(scip, *sol, TRUE, TRUE, FALSE, TRUE, &stored) );
      if( stored )
      {
#ifdef DEBUG_PRIMAL
        printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
        SCIPdebugMessage("found feasible rf solution:\n");
        SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
        //       *result = SCIP_FOUNDSOL;
      }
      else{
        found = 0;
#ifdef DEBUG_RF
        printf("\nCould not found\n. BestUb=%lf", bestUb);
#endif
      }
    }
  }
  //  getchar();
  // Destroi problema
  free(particao);
  free(cand);
  free(fixed);  
  // clear problem
  //  SCIP_CALL( SCIPfree(&subscip) );
  free(solution);
  return found;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRf)
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

   /* solve rf */
   if(rf(scip, &sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
     *result = SCIP_DIDNOTFIND;
#ifdef DEBUG_PRIMAL
     printf("\nRf could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the rf_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRf(
   SCIP*                 scip,                /**< SCIP data structure */
   const parametersT* param
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create rf primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyRf, heurFreeRf, heurInitRf, heurExitRf, heurInitsolRf, heurExitsolRf, heurExecRf,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRf, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRf) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRf) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRf) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRf) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRf) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRf) );
#endif

   /* add rf primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
