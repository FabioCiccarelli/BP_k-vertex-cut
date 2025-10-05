/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   our_code/src/cmain.c
 * @brief  Main file for k-vertex cut pricing example
 * @author Fabio Ciccarelli
 *
 *  This file contains the \ref main() main function of the project. This includes all the default plugins of
 *  \SCIP and the ones which belong to the k-vertex cut project. After that it starts the interactive shell of 
 *  \SCIP or processes the shell arguments if given.
 */
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
#include "propagator.h"

#include "global_variables.h"

using namespace std;

// #define SERVER_CONDOR

string getStatusString(SCIP_STATUS status) {
    switch (status) {
        case SCIP_STATUS_OPTIMAL: return "Optimal";
        case SCIP_STATUS_INFEASIBLE: return "Infeasible";
        case SCIP_STATUS_UNBOUNDED: return "Unbounded";
        case SCIP_STATUS_TIMELIMIT: return "TimeLimit";
        case SCIP_STATUS_USERINTERRUPT: return "UserInterrupt";
        // Aggiungi altri casi se necessario
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
   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL(SCIPincludePropInterdiction(scip));
 
   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   //SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1LL) );


   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   
   /* disable presolving for now to see the original formulation */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );


   SCIP_CALL( SCIPaddIntParam(scip,
      "k",                      
      "Value of k in the k-vertex cut problem", 
      NULL,                        
      FALSE,                    
      2,                         
      1,                        
      INT_MAX,                   
      NULL, NULL) );

   SCIP_CALL( SCIPsetBoolParam(scip, "misc/outputorigsol", FALSE) );
 
   /**********************************
    * Process command line arguments *
   **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );


    /**********************************
    * Write results in a file *
   **********************************/
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
   

   ofstream resultsFile(resultsFileName, ios::app);
   if (resultsFile.is_open())
   {

      const char* targetPricerName = "kvertexcut"; 
      SCIP_PRICER* targetPricer = SCIPfindPricer(scip, targetPricerName);

      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

      resultsFile 
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
      resultsFile.close();

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


   if(SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL)
   {
      char solutionsFileName[100];
   
   #ifdef SERVER_CONDOR
      sprintf(solutionsFileName, "/home/fciccarelli/BP_kvertexcut/branch-and-price/solutions/%s_k%d.sol", SCIPgetProbName(scip), k_val);
   #else
      sprintf(solutionsFileName, "../solutions/%s_k%d.sol", SCIPgetProbName(scip), k_val);
   #endif

      if (SCIPgetBestSol(scip) != NULL) {
         FILE* solFile = fopen(solutionsFileName, "w");
         if (solFile != NULL) {
            SCIP_CALL(SCIPprintBestSol(scip, solFile, FALSE));
            fclose(solFile);
            cout << "Solution written to " << solutionsFileName << endl;
         } else {
            cout << "Error: could not open solution file " << solutionsFileName << endl;
         }
      } 

      // Create a second solution file, with .xsol extension, containing only the indices of the x variables set to 1
      
   #ifdef SERVER_CONDOR
      sprintf(solutionsFileName, "/home/fciccarelli/BP_kvertexcut/branch-and-price/solutions/%s_k%d.xsol", SCIPgetProbName(scip), k_val);
   #else
      sprintf(solutionsFileName, "../solutions/%s_k%d.xsol", SCIPgetProbName(scip), k_val);
   #endif

      if (SCIPgetBestSol(scip) != NULL) {
         int dummy = 0;
         FILE* solFile = fopen(solutionsFileName, "w");
         if (solFile != NULL) {
            SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
            int nnodes = SCIPprobdataGetNNodes(probdata);
            SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);
            for (int v = 0; v < nnodes; ++v) {
               if (SCIPgetSolVal(scip, SCIPgetBestSol(scip), x_vars[v]) > 0.5) {
                  dummy++;
               }
            }
            fprintf(solFile, "%d\n", dummy);
            for (int v = 0; v < nnodes; ++v) {
               if (SCIPgetSolVal(scip, SCIPgetBestSol(scip), x_vars[v]) > 0.5) {
                  fprintf(solFile, "%d\t", v);
               }
            }
            fclose(solFile);
            cout << "Compact solution written to " << solutionsFileName << "\n\n";
         } else {
            cout << "Error: could not open compact solution file " << solutionsFileName << "\n\n";
         }
      }

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

         
   }



   


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

