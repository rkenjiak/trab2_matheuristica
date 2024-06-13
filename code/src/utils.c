#include <time.h>
#include "utils.h"
//
SCIP_RETCODE SCIPtrySolMine(SCIP* scip, SCIP_SOL* sol, SCIP_Bool printreason, SCIP_Bool checkbounds, SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool *stored)
{
#if (defined SCIP_VERSION_MAJOR)  ||  (SCIP_VERSION_MAJOR==6)
   return SCIPtrySol(scip,sol, TRUE,printreason,checkbounds,checkintegrality,checklprows,stored);
#else
   return SCIPtrySol(scip,sol,printreason,checkbounds,checkintegrality,checklprows,stored);
#endif
}
//
SCIP_RETCODE printStatistic(SCIP* scip, double time, char* outputname)
{
  SCIP_Bool outputorigsol = TRUE;
  SCIP_SOL* bestSolution = NULL;
  SCIP_PROBDATA* probdata;
  char filename[SCIP_MAXSTRLEN];
  FILE* fout;
  SCIP_HEUR* heur_hdlr;
  instanceT* I;
  const parametersT* param;

  sprintf(filename, "%s.out", outputname);
  fout = fopen(filename,"w");
  if(!fout){
    printf("\nProblem to create file %s\n", filename);
    return 1;
  }
 
  /* I found the commands for those statistical information looking the scip source code at file scip.c (printPricerStatistics(), for instance)  */ 
  bestSolution = SCIPgetBestSol(scip); 
  if ( outputorigsol )
  {
    if ( bestSolution == NULL )
      printf("\nno solution available\n");
    else
    {
      SCIP_SOL* origsol;
      SCIP_CALL( SCIPcreateSolCopy(scip, &origsol, bestSolution) );
      SCIP_CALL( SCIPretransformSol(scip, origsol) );
      SCIP_CALL( SCIPprintSol(scip, origsol, NULL, FALSE) );
      SCIP_CALL( SCIPfreeSol(scip, &origsol) );
    }
  }
  else
  {
    SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
  }
  SCIPinfoMessage(scip, NULL, "\nStatistics\n");
  SCIPinfoMessage(scip, NULL, "==========\n\n");
  SCIP_CALL( SCIPprintStatistics(scip, NULL) );

  probdata = SCIPgetProbData(scip);
  assert(probdata != NULL);

  I = SCIPprobdataGetInstance(probdata);
  param = SCIPprobdataGetParam(probdata);
  
  if(fout!=NULL){
     fprintf(fout,"%s;%d;%d;%d;%lli;%lf;%lf;%lf;%lf;%lf;%lli;%d;%lf;%lf;%lli;%d;%d", SCIPgetProbName(scip),I->n,I->nS,I->C,SCIPgetNRootLPIterations(scip), time, SCIPgetDualbound(scip), SCIPgetPrimalbound(scip), SCIPgetGap(scip), SCIPgetDualboundRoot(scip), SCIPgetNTotalNodes(scip), SCIPgetNNodesLeft(scip), SCIPgetSolvingTime(scip),SCIPgetTotalTime(scip),SCIPgetMemUsed(scip),SCIPgetNLPCols(scip),SCIPgetStatus(scip));
    if(bestSolution!=NULL){
      fprintf(fout, ";bestsol in %lld;%lf;%d;%s", SCIPsolGetNodenum(bestSolution), SCIPsolGetTime(bestSolution), SCIPsolGetDepth(bestSolution), SCIPsolGetHeur(bestSolution) != NULL ? SCIPheurGetName(SCIPsolGetHeur(bestSolution)) : (SCIPsolGetRunnum(bestSolution) == 0 ? "initial" : "relaxation"));
    }
    else{
      fprintf(fout, ";bestsol in *;0;0;*");
    }
    if(param->heur_aleatoria){
       heur_hdlr = SCIPfindHeur(scip, "aleatoria");
       fprintf(fout, ";%lf;%lld;%lld;%lld;%s",SCIPheurGetTime(heur_hdlr),SCIPheurGetNCalls(heur_hdlr), SCIPheurGetNSolsFound(heur_hdlr), SCIPheurGetNBestSolsFound(heur_hdlr), SCIPheurGetName(heur_hdlr));
    }
    if(param->heur_rf){
       heur_hdlr = SCIPfindHeur(scip, "rf");
       fprintf(fout, ";%lf;%lld;%lld;%lld;%s",SCIPheurGetTime(heur_hdlr),SCIPheurGetNCalls(heur_hdlr), SCIPheurGetNSolsFound(heur_hdlr), SCIPheurGetNBestSolsFound(heur_hdlr), SCIPheurGetName(heur_hdlr));
    }
    if(param->heur_lns){
       heur_hdlr = SCIPfindHeur(scip, "lns");
       fprintf(fout, ";%lf;%lld;%lld;%lld;%s",SCIPheurGetTime(heur_hdlr),SCIPheurGetNCalls(heur_hdlr), SCIPheurGetNSolsFound(heur_hdlr), SCIPheurGetNBestSolsFound(heur_hdlr), SCIPheurGetName(heur_hdlr));
    }
    
    fprintf(fout, ";%s\n", param->parameter_stamp);
  }
  fclose(fout);
  return SCIP_OKAY;
}

/** 
 * creates a SCIP instance with default plugins, and set SCIP parameters 
 */
SCIP_RETCODE configScip(
   SCIP** pscip,
   const parametersT* param
   )
{
   SCIP* scip = NULL;
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) ); 
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   /* for column generation, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) ); 
   /* disable presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) ); // turn off
   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) ); // turn off
   /* disable heuristics */
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );  // turn off
   /* for column generation, usualy we prefer branching using pscost instead of relcost  */
   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority", 1000000) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", param->display_freq) );
   /* set time limit */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", param->time_limit) );
   // for only root, use 1
   SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", param->nodes_limit) );
   if(param->heur_aleatoria)
      SCIP_CALL( SCIPincludeHeurAleatoria(scip, param) );
   
   if(param->heur_rf)
      SCIP_CALL( SCIPincludeHeurRf(scip, param) );

   if(param->heur_lns)
      SCIP_CALL( SCIPincludeHeurLns(scip, param) );

   *pscip = scip;
   return SCIP_OKAY;
}
/**
 * set default+user parameters
 **/
int setParameters(int argc, char** argv, parametersT* param)
{
  typedef struct{
    const char* description;
    const char* param_name;
    void* param_var;
    enum {INT, DOUBLE, STRING} type;
    int ilb;
    int iub;
    double dlb;
    double dub;
    int idefault;
    double ddefault;    
  } settingsT;

  enum {time_limit,display_freq,nodes_limit,param_stamp, param_output_path, heur_freq, heur_maxdepth, heur_freqofs,heur_rf, rf_perc, rf_time, heur_lns, lns_time, lns_perc, heur_aleatoria, total_parameters};

  settingsT parameters[]={
            {"time limit", "--time", &(param->time_limit), INT, 0, 7200, 0,0,1800,0},
            {"display freq", "--display", &(param->display_freq), INT, -1, MAXINT, 0,0,50,0},
            {"nodes limit", "--nodes", &(param->nodes_limit), INT, -1, MAXINT, 0,0,-1,0},
            {"param stamp", "--param_stamp", &(param->parameter_stamp), STRING, 0,0,0,0,0,0},
            {"output path", "--output_path", &(param->output_path), STRING, 0,0,0,0,0,0},
            {"heur_freq", "--heur_freq", &(param->heur_freq), INT, 0,10,0,0,1,0},
            {"heur_maxdepth", "--heur_maxdepth", &(param->heur_maxdepth), INT, -1,10,0,0,-1,0},
            {"heur_freqofs", "--heur_freqofs", &(param->heur_freqofs), INT, 0,10,0,0,0,0},
            {"heur_rf", "--heur_rf", &(param->heur_rf), INT, 0,3,0,0,0,0},
            {"rf perc", "--rf_perc", &(param->rf_perc), DOUBLE, 0,0,0,1.0,0,0.3},
            {"rf_time", "--rf_time", &(param->rf_time), INT, 0, 3600, 0,0,30,0},
            {"heur lns", "--heur_lns", &(param->heur_lns), INT, 0,1,0,0,0,0},
            {"lns perc", "--lns_perc", &(param->lns_perc), DOUBLE, 0,0,0,1.0,0,0.3},
            {"lns_time", "--lns_time", &(param->lns_time), INT, 0, 3600, 0,0,30,0},
            {"heur aleatoria", "--heur_aleatoria", &(param->heur_aleatoria), INT, 0,1,0,0,0,0}
  };
  int i, j, ivalue, error;
  double dvalue;
  FILE *fin;
  char* current_path = ".";

  // total_parameters = sizeof(parameters)/sizeof(parameters[0]);
  
  
  if (param==NULL)
    return 0;
  
  // check arguments
  if(argc<2){
    printf("\nSintaxe: ./bin/program <instance-file> <parameters-setting>.\n\t or Use ./bin/program --options to show options to parameters settings.\nExample of usage:\n\t ./bin/program data/t50-10-1000-2-s-1.mochila\n\t ./bin/program data/t50-10-1000-2-s-1.mochila --param_stamp default.config\n\nIf no param_stamp is given by user, a new param stamp named dAAAAMMDDhHHMMSS will be created.\n\nIf the given param stamp is new (it does not exist in the current folder), it will be created to save all chosen parameters settings. Otherwise, if the param stamp file already exists, all saved parameters settings in the file must be the same as those parameters settings given in the command line.\n\nP.S.: To use a given param stamp file ""x.config"", please use the command xargs as follows:\n\n \t xargs -a x.config ./bin/program data/t50-10-1000-2-s-1.mochila\n");
    return 0;
  }
  else if(argc==2 && !strcmp(argv[1],"--options")){  // show options
     //    showOptions();
     printf("\noption                : default -   range        : description");
     for(i=0;i<total_parameters;i++){
        switch(parameters[i].type){
        case INT:
           printf("\n%-22s: %7d - [%3d,%8d] : %s", parameters[i].param_name, parameters[i].idefault, parameters[i].ilb, parameters[i].iub,parameters[i].description);
           break;
        case DOUBLE:
           printf("\n%-22s: %7.1lf - [%3.1lf,%8.1lf] : %s", parameters[i].param_name, parameters[i].ddefault, parameters[i].dlb, parameters[i].dub,parameters[i].description);
           break;
        case STRING:
           printf("\n%-22s:       * - [*,*]          : %s", parameters[i].param_name, parameters[i].description);
        }
     }
     printf("\n");
     return 0;
  }
  
  // check the existance of instance file
  fin = fopen(argv[1], "r");
  if(!fin){
      printf("\nInstance file not found: %s\n", argv[1]);
      return 0;
  }
  fclose(fin);  

  // set default parameters value
  for(i=0;i<total_parameters;i++){
    if(parameters[i].type==INT)
      *((int*)(parameters[i].param_var)) = parameters[i].idefault;
    else if (parameters[i].type==DOUBLE)
      *((double*)(parameters[i].param_var)) = parameters[i].ddefault;
    else
      *((char**) (parameters[i].param_var)) = NULL;
  }
  param->output_path = current_path;

  // set user parameters value
  error = 0;
  for(i=2;i<argc && !error;i+=2){
    for(j=0;j<total_parameters && strcmp(argv[i],parameters[j].param_name);j++)
      ;
    if(j>=total_parameters || i==argc-1){
      printf("\nParameter (%s) invalid or uncompleted.", argv[i]);
      error = 1;
    }
    else{
      switch(parameters[j].type){
      case INT:
        ivalue = atoi(argv[i+1]);
        if(ivalue < parameters[j].ilb || ivalue > parameters[j].iub){
          printf("\nParameter (%s) value (%d) out of range [%d,%d].", argv[i], ivalue, parameters[j].ilb, parameters[j].iub);
          error = 1;
          //          break;
        }
        else{
          *((int*)(parameters[j].param_var)) = ivalue;          
        }
        break;
      case DOUBLE:
        dvalue = atof(argv[i+1]);
        if(dvalue < parameters[j].dlb || dvalue > parameters[j].dub){
          printf("\nParameter (%s) value (%lf) out of range [%lf,%lf].", argv[i], dvalue, parameters[j].dlb, parameters[j].dub);
          error = 1;
          //          break;
        }
        else{
          *((double*)(parameters[j].param_var)) = dvalue;          
        }
        break;
      case STRING:
        *((char**) (parameters[j].param_var)) = argv[i+1];
      }
    }
  }

  // print parameters
  printf("\n\n----------------------------\nParameters settings");
  for(i=0;i<total_parameters;i++){
    printf("\nparameter %s (%s): default value=", parameters[i].param_name, parameters[i].description);
    switch(parameters[i].type){
    case INT:
      printf("%d - value = %d", parameters[i].idefault, *( (int*)parameters[i].param_var));
      break;
    case DOUBLE:
      printf("%lf - value = %lf", parameters[i].ddefault, *( (double*)parameters[i].param_var));
      break;
    case STRING:
      printf("(null) - value = %s", *( (char**)parameters[i].param_var));
      break;
    }
  }
  printf("\nerror = %d\n", error);
  if(!error){
    FILE* fout;
    char foutname[SCIP_MAXSTRLEN];

    if(param->parameter_stamp != NULL){
      // complete the fullname of the parameters stamp file
       sprintf(foutname, "%s/%s", param->output_path, param->parameter_stamp);
    }
    else{
      // define the parameters stamp file using stamp default = date-time
      struct tm * ct;
      const time_t t = time(NULL);
      param->parameter_stamp = (char*) malloc(sizeof(char)*100);
      ct = localtime(&t);
      snprintf(param->parameter_stamp, 100, "d%d%.2d%.2dh%.2d%.2d%.2d", ct->tm_year+1900, ct->tm_mon, ct->tm_mday, ct->tm_hour, ct->tm_min, ct->tm_sec);

      sprintf(foutname, "%s/%s", param->output_path, param->parameter_stamp);
    }
    // check if the stamp already exists
    fout = fopen(foutname, "r");
    // TODO: opendir() should be done first to avoid open a directory!
    if(!fout){
      // save parameters in the stamp file
      printf("\nwriting parameters in %s", foutname);
      fout = fopen(foutname, "w+");
      for(i=0;i<total_parameters;i++){
        fprintf(fout, "%s ", parameters[i].param_name);
        switch(parameters[i].type){
        case INT:
          fprintf(fout,"%d\n", *( (int*)parameters[i].param_var));
          break;
        case DOUBLE:
          fprintf(fout,"%lf\n", *( (double*)parameters[i].param_var));
          break;
        case STRING:
          fprintf(fout,"%s\n", *( (char**)parameters[i].param_var));
          break;
        }
      }
      fclose(fout);
    }
    else{
      // check if the stamp is valid
      char param_name[100];
      char svalue[100];
      while(!feof(fout)){
        fscanf(fout,"%s", param_name);
        for(j=0;j<total_parameters && strcmp(param_name,parameters[j].param_name);j++)
          ;
        if(j>=total_parameters){
          printf("\nParameter (%s) invalid or uncompleted.", param_name);
          error = 1;
          break;
        }
        else{
          switch(parameters[j].type){
          case INT:
            fscanf(fout, "%d\n", &ivalue);
            if(ivalue < parameters[j].ilb || ivalue > parameters[j].iub){
              printf("\nParameter (%s) value (%d) out of range [%d,%d].", param_name, ivalue, parameters[j].ilb, parameters[j].iub);
              error = 1;
            }
            else if(ivalue != *((int*)(parameters[j].param_var))){
              printf("\nParameter (%s) value (%d) differs to saved value = %d.", param_name, ivalue, *((int*)(parameters[j].param_var)));
              error = 1;
            }
            break;
          case DOUBLE:
            fscanf(fout, "%lf\n", &dvalue);
            if(dvalue < parameters[j].dlb || dvalue > parameters[j].dub){
              printf("\nParameter (%s) value (%lf) out of range [%lf,%lf].", param_name, dvalue, parameters[j].dlb, parameters[j].dub);
              error = 1;
            }
            else if(fabs(dvalue - *((double*)(parameters[j].param_var))) > EPSILON){
              printf("\nParameter (%s) value (%lf) differs to saved value = %lf.", param_name, dvalue, *((double*)(parameters[j].param_var)));
              error = 1;
            }
            break;
          case STRING:
            fscanf(fout, "%s\n", svalue);
            if(strcmp(svalue,*((char**)(parameters[j].param_var)))){
              printf("\nParameter (%s) value (%s) differs to saved value = %s.", param_name, svalue, *((char**)(parameters[j].param_var)));
              error = 1;
            }
            break;
          }
        } // each parameter
      } // while !feof
      fclose(fout);
    }
  }
  return !error;
}
// TODO: Get the best solution found and write the solution in a file. It depends on the problem!
void printSol(SCIP* scip, char* outputname)
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* bestSolution;
   SCIP_VAR** vars;
   SCIP_Real solval;
   FILE *file;
   int v;//, nvars;
   instanceT* I;
   char filename[SCIP_MAXSTRLEN];
   const parametersT* param;
   struct tm * ct;
   const time_t t = time(NULL);

   assert(scip != NULL);
   bestSolution = SCIPgetBestSol(scip);
   if( bestSolution == NULL )
     return;
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   I = SCIPprobdataGetInstance(probdata);
   //   nvars = SCIPprobdataGetNVars(probdata);
   vars = SCIPprobdataGetVars(probdata);
   param = SCIPprobdataGetParam(probdata);

   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.sol", outputname);
   file = fopen(filename, "w");
   if(!file)
     {
       printf("\nProblem to create solution file: %s", filename);
       return;
     } 
   fprintf(file, "\nValue: %lf\nItems: ", -SCIPsolGetOrigObj(bestSolution));

   for( v=0; v< I->n; v++ )
     {
       solval = SCIPgetSolVal(scip, bestSolution, vars[v]);
       if( solval > EPSILON )
	 {
	   fprintf(file, "%d ", I->item[v].label);
	 }
     }

   fprintf(file, "\n");
   //
   fprintf(file, "Parameters settings file=%s\n", param->parameter_stamp);
   fprintf(file, "Instance file=%s\n", SCIPgetProbName(scip));
   ct = localtime(&t);
   fprintf(file, "Date=%d-%.2d-%.2d\nTime=%.2d:%.2d:%.2d\n", ct->tm_year+1900, ct->tm_mon, ct->tm_mday, ct->tm_hour, ct->tm_min, ct->tm_sec);
   fclose(file);
}

void removePath(char* fullfilename, char** filename)
{
 // remove path, if there exists on fullfilename
 *filename = strrchr(fullfilename, '/');
 if(*filename==NULL){
   *filename = fullfilename;
 }
 else{
   (*filename)++; // discard /
 }
}
void configOutputName(char* name, char* instance_filename, char* program, const parametersT* param)
{
  char* program_filename, *filename;

  // remove path, if there exists on program name
  removePath(program, &program_filename);
  removePath(instance_filename, &filename);

 // append program name and parameter stamp
  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s/%s-%s-%s", param->output_path,filename, program_filename, param->parameter_stamp);
}
/* sorteia um numero aleatorio entre [low,high] */
int RandomInteger(int low, int high)
{
  int k;
  double d;
  //  srand(time(NULL));

  d = (double)rand() / ((double)RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}
// Função auxiliar de comparacao para o qsort
int comparador(const void *valor1, const void *valor2)
{
  int k=4; // trocar k para testar diferentes sorts value, weight
  switch (k)
  {
    case 0:
      if ((*(itemType *)valor1).value > (*(itemType *)valor2).value) return 1;
      else if ((*(itemType *)valor1).value == (*(itemType *)valor2).value) return 0;
      else return -1;
      break;
    case 1:
      if ((*(itemType *)valor1).value < (*(itemType *)valor2).value) return 1;
      else if ((*(itemType *)valor1).value == (*(itemType *)valor2).value) return 0;
      else return -1;
      break;
    case 2:
      if ((*(itemType *)valor1).weight > (*(itemType *)valor2).weight) return 1;
      else if ((*(itemType *)valor1).weight == (*(itemType *)valor2).weight) return 0;
      else return -1;
      break;
    case 3: // tira maior peso
      if ((*(itemType *)valor1).weight < (*(itemType *)valor2).weight) return 1;
      else if ((*(itemType *)valor1).weight == (*(itemType *)valor2).weight) return 0;
      else return -1;
      break; 
    case 4:
      if ((float)((*(itemType *)valor1).value)/((*(itemType *)valor1).weight) > ((float)(*(itemType *)valor2).value)/(*(itemType *)valor2).weight) return 1;
      else if ((float)((*(itemType *)valor1).value)/((*(itemType *)valor1).weight) == ((float)(*(itemType *)valor2).value)/(*(itemType *)valor2).weight) return 0;
      else return -1;
      break;
    case 5:
      if ((float)((*(itemType *)valor1).value)/((*(itemType *)valor1).weight) < ((float)(*(itemType *)valor2).value)/(*(itemType *)valor2).weight) return 1;
      else if ((float)((*(itemType *)valor1).value)/((*(itemType *)valor1).weight) == ((float)(*(itemType *)valor2).value)/(*(itemType *)valor2).weight) return 0;
      else return -1;
      break; 
  }  
}


