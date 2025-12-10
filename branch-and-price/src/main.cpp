/**@file   main.cpp
 * @brief  Main file of our algorithm for the k-vertex cut pricing problem
 * @author Fabio Ciccarelli
 *
 *  This file contains the \ref main() main function of the project. This includes all the default plugins of
 *  \SCIP and the ones which belong to the k-vertex cut project. 
 **/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "pricer_kvertexcut.h"
#include "reader_kvertexcut.h"
#include "probdata_kvertexcut.h"
<<<<<<< Updated upstream
#include "propagator.h"

#include "global_variables.h"

using namespace std;

// #define SERVER_CONDOR
=======
#include "plugins.h"

using namespace std;

>>>>>>> Stashed changes

string getStatusString(SCIP_STATUS status) {
    switch (status) {
        case SCIP_STATUS_OPTIMAL: return "Optimal";
        case SCIP_STATUS_INFEASIBLE: return "Infeasible";
        case SCIP_STATUS_UNBOUNDED: return "Unbounded";
        case SCIP_STATUS_TIMELIMIT: return "TimeLimit";
        case SCIP_STATUS_USERINTERRUPT: return "UserInterrupt";
        default: return "Unknown";
    }
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   
   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include k-vertex cut reader */
   SCIP_CALL( SCIPincludeReaderKvc(scip) );
   
   /* include k-vertex cut pricer */
   SCIP_CALL( SCIPincludePricerKvertexcut(scip) );
<<<<<<< Updated upstream
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
=======
   
   /* include SCIP plugins */
   SCIP_CALL( SCIPincludePlugins(scip) );
>>>>>>> Stashed changes

   SCIP_CALL(SCIPincludePropInterdiction(scip));
 
   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   //SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1LL) );


   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   
   /* disable presolving for now to see the original formulation */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );

<<<<<<< Updated upstream
=======
   /* disable output of original solution */
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/outputorigsol", FALSE) );

   
   /* PARAMETERS SPECIFIC FOR THE K-VERTEX CUT PROBLEM */

   SCIP_CALL( SCIPaddStringParam(scip,
      "output/results",                      
      "Path of the file for output results", 
      NULL,                        
      FALSE,                    
      ".",                                          
      NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "output/solution",                      
      "Path of the file for output solution", 
      NULL,                        
      FALSE,                    
      ".",                                          
      NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "output/plot",                      
      "Path of the file for output plot", 
      NULL,                        
      FALSE,                    
      ".",                                          
      NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/plot",                      
      "option to enable plotting of the k-vertex cut problem solution (1: enable, 0: disable)", 
      NULL,                        
      FALSE,                    
      0,                                                                   
      NULL, NULL) );
>>>>>>> Stashed changes

   SCIP_CALL( SCIPaddIntParam(scip,
      "params/k",                      
      "Value of k in the k-vertex cut problem", 
      NULL,                        
      FALSE,                    
      2,                         
      1,                        
      INT_MAX,                   
      NULL, NULL) );

<<<<<<< Updated upstream
   SCIP_CALL( SCIPsetBoolParam(scip, "misc/outputorigsol", FALSE) );
 
=======
   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/weighted",                      
      "option to solve the weighted version of the k-vertex cut problem (1: weighted, 0: unweighted)", 
      NULL,                        
      FALSE,                    
      0,                                                                   
      NULL, NULL) );


   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/connectivity/warmstart",                      
      "option to enable connectivity based warmstart (1: enable, 0: disable)", 
      NULL,                        
      FALSE,                    
      1,                                                                   
      NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/ILPwarmstart",                      
      "option to enable ILP based warmstart (1: enable, 0: disable)", 
      NULL,                        
      FALSE,                    
      0,                                                                   
      NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/connectivity/cut",                      
      "option to enable minimum connectivity cut (1: enable, 0: disable)", 
      NULL,                        
      FALSE,                    
      1,                                                                   
      NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "options/solveLP",                      
      "option to solve only the LP relaxation of the k-vertex cut problem (1: enable, 0: disable)", 
      NULL,                        
      FALSE,                    
      0,                                                                   
      NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, 
      "options/symmetrymethod",
      "symmetry handling method (0: none, 1: symresacks, 2: lex. red./orbital red.",
      NULL, 
      FALSE, 
      2, 
      0, 
      2, 
      NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, 
      "options/cliquegeneration",
      "clique generation method (0: single edges, 1: greedy maximal cliques",
      NULL, 
      FALSE, 
      2, 
      0, 
      2, 
      NULL, NULL) );

   SCIP_CALL(SCIPaddBoolParam(scip,
      "options/relaxconstraints",
      "whether to relax the equality constraints to >= (1: yes, 0: no)",
      NULL,
      FALSE,
      1,
      NULL, NULL));
   /**************************************************************************************************** */

   
   unsigned int LPonly;
   SCIP_CALL( SCIPgetBoolParam(scip, "options/solveLP", &LPonly) );

   if(LPonly) {
      SCIP_CALL( SCIPsetIntParam(scip, "limits/nodes", 1) );
   }

>>>>>>> Stashed changes
   /**********************************
    * Process command line arguments *
   **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );


    /**********************************
    * Write results in a file *
   **********************************/
<<<<<<< Updated upstream
   // check if a solution has been found
   // the results file has name SCIPgetProbName(scip)_results_k{k value}.txt

   char resultsFileName[100];
   int k_val;

   SCIPgetIntParam(scip, "k", &k_val);

#ifdef SERVER_CONDOR
   sprintf(resultsFileName, "/home/fciccarelli/BP_kvertexcut/branch-and-price/results/%s_k%d.res", SCIPgetProbName(scip), k_val);
#else
   sprintf(resultsFileName, "../results/%s_k%d.res", SCIPgetProbName(scip), k_val);
#endif
=======
   int k_val;
   int sym_method;
   int clique_gen_method;
   unsigned int plot_option;
   unsigned int option_conn_warmstart;
   unsigned int option_conn_cut;
   unsigned int weighted;
   char* resultsFileName;
   char* solutionsFileName;
   char* plotFileName;
   unsigned int relax_constraints;
   unsigned int option_ILP_warmstart;

   SCIPgetStringParam(scip, "output/results", &resultsFileName);
   SCIPgetStringParam(scip, "output/solution", &solutionsFileName);
   SCIPgetStringParam(scip, "output/plot", &plotFileName);
   SCIPgetIntParam(scip, "params/k", &k_val);
   SCIPgetIntParam(scip, "options/symmetrymethod", &sym_method);
   SCIPgetIntParam(scip, "options/cliquegeneration", &clique_gen_method);
   SCIPgetBoolParam(scip, "options/plot", &plot_option);
   SCIPgetBoolParam(scip, "options/connectivity/warmstart", &option_conn_warmstart);
   SCIPgetBoolParam(scip, "options/ILPwarmstart", &option_ILP_warmstart);
   SCIPgetBoolParam(scip, "options/connectivity/cut", &option_conn_cut);
   SCIPgetBoolParam(scip, "options/solveLP", &LPonly);
   SCIPgetBoolParam(scip, "options/weighted", &weighted);
   SCIPgetBoolParam(scip, "options/relaxconstraints", &relax_constraints);

   // add the problem name to the results file path
   if (strlen(resultsFileName) > 0) {
      string resultsFilePath = string(resultsFileName) + "/" + string(SCIPgetProbName(scip)) + "_k" + to_string(k_val) + ".res";
      resultsFileName = new char[resultsFilePath.length() + 1];
      strcpy(resultsFileName, resultsFilePath.c_str());
   } else {
      cout << "No results file path provided." << endl;
      exit(-1);
   }

   if (strlen(solutionsFileName) > 0) {
      string solutionsFilePath = string(solutionsFileName) + "/" + string(SCIPgetProbName(scip)) + "_k" + to_string(k_val) + ".sol";
      solutionsFileName = new char[solutionsFilePath.length() + 1];
      strcpy(solutionsFileName, solutionsFilePath.c_str());
   } else {
      cout << "No solutions file path provided." << endl;
      exit(-1);
   }

   if( strlen(plotFileName) > 0) {
      if(plot_option) {
         string plotFilePath = string(plotFileName) + "/" + string(SCIPgetProbName(scip)) + "_k" + to_string(k_val) + ".png"; 
         plotFileName = new char[plotFilePath.length() + 1];
         strcpy(plotFileName, plotFilePath.c_str());
      }
      else {
         string plotFilePath = string(plotFileName) + "/" + string(SCIPgetProbName(scip)) + "_k" + to_string(k_val); 
         plotFileName = new char[plotFilePath.length() + 1];
         strcpy(plotFileName, plotFilePath.c_str());
      }
   } else {
      cout << "No plot file path provided." << endl;
      exit(-1);
   }
   
   
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

   int tot_cut_cost = 0;
   Graph *orig_graph = SCIPprobdataGetOrigGraph(probdata);
   bool* cut = SCIPprobdataGetSolution(scip, probdata, tot_cut_cost);
   int preFixedCount = SCIPprobdataGetNFixed(probdata);
   bool* preFixed = SCIPprobdataGetPreFixed(probdata);

   int preFixedWeight = 0;
   for(int v = 0; v < orig_graph->nnodes; v++) {
      if(preFixed[v]) {
         preFixedWeight += orig_graph->vertex_weights[v];
      }
   }
>>>>>>> Stashed changes
   
   ofstream resultsFile(resultsFileName, ios::app);
   if (resultsFile.is_open())
   {

      const char* targetPricerName = "kvertexcut"; 
      SCIP_PRICER* targetPricer = SCIPfindPricer(scip, targetPricerName);

      /* Create residual graph for statistics about prefixing effect */
      vector<int> a;
      vector<int> b;
      cout << "\nCreating residual graph..." << endl;
      Graph* without_prefixed = buildResidualGraph(*orig_graph, preFixed, a, b);
      

      resultsFile 
                  << SCIPgetProbName(scip) << "\t"
                  << orig_graph->nnodes << "\t"
                  << orig_graph->nedges << "\t"
                  << without_prefixed->nnodes << "\t"
                  << without_prefixed->nedges << "\t"
                  << k_val << "\t"
<<<<<<< Updated upstream
=======
                  << sym_method << "\t"
                  << option_conn_warmstart << "\t"
                  << option_ILP_warmstart << "\t"
                  << option_conn_cut << "\t"
                  << clique_gen_method << "\t"
                  << relax_constraints << "\t"
                  << LPonly << "\t"
                  << weighted << "\t"
                  << preFixedCount << "\t"
                  << preFixedWeight << "\t"
>>>>>>> Stashed changes
                  << getStatusString(SCIPgetStatus(scip)) << "\t"
                  << SCIPgetPrimalbound(scip) << "\t"
                  << SCIPgetDualbound(scip) << "\t"
                  << tot_cut_cost << "\t"
                  << SCIPgetSolvingTime(scip) << "\t"
                  << SCIPpricerGetTime(targetPricer) << "\t"
                  << SCIPprobdataGetNAlphaVars(probdata) << "\t"
                  << SCIPgetNNodes(scip) << "\t"
                  << SCIPgetNOrigVars(scip) << "\t"
                  << SCIPgetNOrigConss(scip) << "\t"
                  << global_NGenColsRootNode << "\t"
                  << global_NGenFarkasCols << "\t"
                  << global_MaxDepth << "\t"
                  << global_RootNodeLowerbound << "\t"
                  << global_RootNodeUpperbound
                  << endl;
      resultsFile.close();

<<<<<<< Updated upstream
      cout << "\n\n" << SCIPprobdataGetNAlphaVars(probdata) << " variables added to the master problem" << endl;
      cout << SCIPpricerGetTime(targetPricer) << " seconds spent in pricing" << endl;
      cout << SCIPgetNNodes(scip) << " branch-and-bound nodes processed\n\n" << endl;

      // print the global variables
      cout << "Number of columns generated at the root node: " << global_NGenColsRootNode << endl;
      cout << "Number of Farkas columns generated: " << global_NGenFarkasCols << endl;
      cout << "Max depth: " << global_MaxDepth << endl;
      cout << "Root node lower bound: " << global_RootNodeLowerbound << endl;
      cout << "Root node upper bound: " << global_RootNodeUpperbound << endl;
      cout << "-----------------------------------------------------\n\n" << endl;
   }

   // write all the results in a single file named "all_results.res" in append mode
#ifdef SERVER_CONDOR
   ofstream allResultsFile("/home/fciccarelli/BP_kvertexcut/branch-and-price/bin/info_summary.txt", ios::app);
#else
   ofstream allResultsFile("../info_summary.txt", ios::app);
#endif
   if (allResultsFile.is_open())
   {
      const char* targetPricerName = "kvertexcut"; 
      SCIP_PRICER* targetPricer = SCIPfindPricer(scip, targetPricerName);

      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

      allResultsFile 
                  << SCIPgetProbName(scip) << "\t"
                  << SCIPprobdataGetNNodes(probdata) << "\t"
                  << SCIPprobdataGetNEdges(probdata) << "\t"
                  << k_val << "\t"
                  << getStatusString(SCIPgetStatus(scip)) << "\t"
                  << SCIPgetPrimalbound(scip) << "\t"
                  << SCIPgetDualbound(scip) << "\t"
                  << SCIPgetSolvingTime(scip) << "\t"
                  << SCIPpricerGetTime(targetPricer) << "\t"
                  << SCIPprobdataGetNAlphaVars(probdata) << "\t"
                  << SCIPgetNNodes(scip) << "\t"
                  << SCIPgetNOrigVars(scip) << "\t"
                  << SCIPgetNOrigConss(scip) << "\t"
                  << global_NGenColsRootNode << "\t"
                  << global_NGenFarkasCols << "\t"
                  << global_MaxDepth << "\t"
                  << global_RootNodeLowerbound << "\t"
                  << global_RootNodeUpperbound
                  << endl;
      allResultsFile.close();
   }

=======
      delete without_prefixed;
>>>>>>> Stashed changes

   }
   else
   {
      cout << "\nError: could not open results file " << resultsFileName << endl;
   }

   // Write the solution in the format tot_cut_cost \n vertices in cut separated by space
   ofstream solFile(solutionsFileName);
   if (solFile.is_open()) {
      solFile << tot_cut_cost << "\n";
      for (int v = 0; v < orig_graph->nnodes; ++v) {
         if (cut[v]) {
            solFile << v << " ";
         }
      }
      solFile << "\n";
      solFile.close();
      cout << "Solution written to " << solutionsFileName << endl;
   } else {
      cout << "\nError: could not open solution file " << solutionsFileName << endl;

      delete[] resultsFileName;
      delete[] solutionsFileName;
      delete[] plotFileName;
      delete[] cut;

      exit(-1);
   }

<<<<<<< Updated upstream
      int i,j;

      char *probname=new char[1000];
      char *probname2=new char[1000];
      sprintf(probname, "%s.gml",solutionsFileName);
      sprintf(probname2, "%s2.gml",solutionsFileName);

      cout << "Writing gml files " << probname << " and " << probname2 << endl;

      //cout << probname << endl;
      ofstream grafo(probname);
      ofstream grafo2(probname2);

      char *buf_=new char[1000];
      char *buf_COOR=new char[1000];


      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
      int nnodes = SCIPprobdataGetNNodes(probdata);
      int nedges = SCIPprobdataGetNEdges(probdata);

      int* tail = SCIPprobdataGetEdgeTails(probdata);
      int* head = SCIPprobdataGetEdgeHeads(probdata);

      SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);

      grafo  << "graph  [ hierarchic  1  directed  0 \n\n" << endl;
      grafo2 << "graph  [ hierarchic  1  directed  0 \n\n" << endl;

      for (i = 0; i < nnodes; ++i){
         //solo puntini
         sprintf(buf_COOR, " %d ",i);
         
         if(SCIPgetSolVal(scip, SCIPgetBestSol(scip), x_vars[i]) > 0.5)
            {
               sprintf(buf_, "node  [ id  %d label \"%d\" graphics  [ x %f   y %f  w 20 h 20 type \"roundrectangle\" fill	\"#da140dff\" ]  ]  ", i, i, 0.0, 0.0);
               
               grafo2 << buf_ << endl;
               grafo  << buf_ << endl;
               continue;
            }

         sprintf(buf_, "node  [ id  %d label \"%d\" graphics  [ x %f   y %f  w 20 h 20 type \"roundrectangle\" fill	\"#1b9c22ff\" ]  ]  ", i, i, 0.0, 0.0);

         grafo  << buf_ << endl;
         grafo2 << buf_ << endl;
      }

      
      
      for (j = 0; j < nedges; ++j){
         int u = head[j] - 1;
         int v = tail[j] - 1;

         if(SCIPgetSolVal(scip, SCIPgetBestSol(scip), x_vars[v]) > 0.5 ||
            SCIPgetSolVal(scip, SCIPgetBestSol(scip), x_vars[u]) > 0.5 )
            {
               sprintf(buf_, "edge [ source	%d target	%d graphics [ style \"dashed\" fill	\"#000000\" ] ]",
                        tail[j]-1, head[j]-1);
               grafo2 << buf_ << endl;
               continue;
            }
         //cout << vet_x[i] << " " << vet_y[i] << " "<< vet_x[j] << " " << vet_x[j] << endl; 
         sprintf(buf_, "edge [ source	%d target	%d graphics [ fill	\"#000000\" ] ]",
                        tail[j]-1, head[j]-1);

                        grafo  << buf_ << endl; 
                        grafo2 << buf_ << endl;
      }
      
      

      grafo  << "\n] \n\n" << endl;
      grafo2 << "\n] \n\n" << endl;

      grafo.close();
      grafo2.close();

      delete [] probname;
      delete [] probname2;
      delete [] buf_;
      delete [] buf_COOR;
=======
   // Plot the solution graph
   SCIP_CALL( SCIPprobdataPlotSolution(orig_graph, cut, plotFileName, false, plot_option) );
>>>>>>> Stashed changes


   delete[] resultsFileName;
   delete[] solutionsFileName;
   delete[] plotFileName;
   delete[] cut;


   


   /********************
    * Deinitialization *
   ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

