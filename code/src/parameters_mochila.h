#ifndef __SCIP_PARAMETERS_MOCHILA__
#define __SCIP_PARAMETERS_MOCHILA__

#define MAXINT 1000
typedef struct{
   // global settings
   int time_limit; /* limit of execution time (in sec). Default = 1800 (-1: unlimited) */
   int display_freq; /* frequency to display information about B&B enumeration. Default = 50 (-1: never) */
   int nodes_limit; /* limit of nodes to B&B procedure. Default = -1: unlimited (1: onlyrootnode) */

   // parameter stamp
   char* parameter_stamp;
   char* output_path;

   // primal heuristic
   int heur_freq;
   int heur_maxdepth;
   int heur_freqofs;
   int heur_rf;
   double rf_perc;
   int rf_time;
   int heur_lns;
   double lns_perc;
   int lns_time;  

   int heur_aleatoria;
} parametersT;

int setParameters(int argc, char** argv, parametersT* Param);

extern parametersT param;
extern const char* output_path;
#endif
