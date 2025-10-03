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

/**@file   probdata_kvertexcut.cpp
 * @brief  Problem data for k-vertex cut problem
 * @author Fabio Ciccarelli
 *
 * This file handles the main problem data used in that project. For more details see \ref KVERTEXCUT_PROBLEMDATA page.
 *
 * @page KVERTEXCUT_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the k-vertex cut problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
*/

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "gurobi_c++.h"
#include "probdata_kvertexcut.h"
#include "vardata_kvertexcut.h"
#include "pricer_kvertexcut.h"

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"

using namespace std;

/** forward declarations */
static SCIP_RETCODE probdataFree(SCIP* scip, SCIP_PROBDATA** probdata);


/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the k-vertex-cut, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   int                   nnodes;             /**< number of nodes in the graph */
   int                   nedges;             /**< number of edges in the graph */
   int*                  tail;                /**< array of edge tail nodes */
   int*                  head;                /**< array of edge head nodes */
   int**                 adj_lists;          /**< adjacency lists for each node */
   int*                  adj_sizes;          /**< size of adjacency list for each node */
   SCIP_Bool**           adj_matrix;         /**< adjacency matrix */
   int*                  node_degree;       /**< degree of each node */
   int                   k;                  /**< parameter k for k-vertex cut */

   SCIP_VAR**            x_vars;             /**< array of node variables */
   SCIP_VAR**            z_v_vars;           /**< array of z_v variables */
   SCIP_VAR**            z_u_vars;           /**< array of z_u variables */
   SCIP_VAR**            beta_v_vars;        /**< array of beta_v variables */
   SCIP_VAR**            beta_u_vars;        /**< array of beta_u variables */

   SCIP_VAR**            alphavars;               /**< array of variables */
   int                   n_alphavars;             /**< number of variables */
   int                   alphavars_size;          /**< size of variables array */

   SCIP_CONS*            main_alpha_constr;          /**< main pricing constraint for α_S */
   SCIP_CONS**           alpha_constrs;              /**< array of pricing constraints for α_S */
   SCIP_CONS**           coverage_constrs;         /**< array of SCC coverage constraints */
   int                   nalphaconstrs;             /**< number of pricing constraints for α_S */


   SCIP_CONS**           other_contrs;        /**< array of other linear constraints*/
   int                   notherconstrs;        /**< number of other linear constraints */

   int                   root_generated_vars;              /**< number of variables generated at root */
   int                   farkas_generated_vars;             /**< number of generated variables with farkas pricing */
   double                root_lower_bound;                 /**< value of the LP relaxation */
   double                root_upper_bound;                  /**< best incumbent value at the root node */
};

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   SCIPdebugMsg(scip, "exec method of event handler for added variable to probdata\n");
   //cout << "exec method of event handler for added variable to probdata\n";

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event)) );

   return SCIP_OKAY;
}

// Il nostro gestore di eventi con logica personalizzata (versione corretta)
static
SCIP_DECL_EVENTEXEC(eventExecMyLogic)
{
    assert(eventhdlr != NULL);
    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), "DENTRO_NODO") == 0);
    assert(event != NULL);
    assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEFOCUSED);

    SCIP_NODE* node = SCIPgetCurrentNode(scip);
    int num_node = SCIPnodeGetNumber(node);

    // Update global_MaxDepth if needed
   int depth = SCIPnodeGetDepth(node);
   if( depth > global_MaxDepth )
   {
      global_MaxDepth = depth;
   }

   if((num_node == 2 || num_node == 3) && global_RootNodeLowerbound == -1){

      // Get the father node
      SCIP_NODE* parent_node = SCIPnodeGetParent(node);

      SCIP_Real node_lowerbound = SCIPnodeGetLowerbound(parent_node);
      SCIP_Real root_incumbent = SCIPgetPrimalbound(scip);

      global_RootNodeLowerbound = (double)node_lowerbound;
      global_RootNodeUpperbound = (double)root_incumbent;
   }

    // get local UB and LB for x variables
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

    // get alpha variables and their count
    SCIP_VAR** alpha_vars = SCIPprobdataGetAlphaVars(probdata);
    int n_alpha_vars = SCIPprobdataGetNAlphaVars(probdata);


    // print the number of the node we are at
    // cout << "FOCUSING ON NODE NUMBER: " << SCIPnodeGetNumber(node) << endl;

    // get domain changes for the node
    SCIP_DOMCHG* domchg = SCIPnodeGetDomchg(node);
    
    if( domchg != NULL )
    {
        
      SCIP_BOUNDCHG* boundchg = SCIPdomchgGetBoundchg(domchg, 0);
      SCIP_VAR* var = SCIPboundchgGetVar(boundchg);
      SCIP_Real newbound = SCIPboundchgGetNewbound(boundchg);

      const char* varname = SCIPvarGetName(var);
      
      if( strncmp(varname, "t_x_", 4) == 0 )
      {
         char* endptr;
         long vertex_index = strtol(varname + 4, &endptr, 10);

         // print which variable has been changed and its new bound
         // cout << "VARIABLE " << varname << " FIXED TO " << newbound << endl;
         
         // check if alpha variables contain the vertex whose variable has been fixed
         for( int i = 0; i < n_alpha_vars; ++i )
         {
               SCIP_VAR* alpha_var = alpha_vars[i];
               SCIP_VARDATA* vardata = SCIPvarGetData(alpha_var);
               
               if( vardata != NULL )
               {
                  // Get information about the subset
                  int* subset = SCIPvardataGetSubset(vardata);

                  // Check if vertex_index is present in the subset
                  if(newbound == 1.0 && subset[vertex_index] == 1 && SCIPisGT(scip, SCIPvarGetUbLocal(alpha_var), 0.0) )
                  {
                     SCIP_CALL( SCIPchgVarUb(scip, alpha_var, 0.0) );
                  }
                  
                  
               }
         }
        
      }
   }
     

   return SCIP_OKAY;
}



/** creates transformed problem data */
static
SCIP_RETCODE probdataTrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       targetdata,           /**< pointer to problem data */
   int                   nnodes,             /**< number of nodes in the graph */
   int                   nedges,             /**< number of edges in the graph */
   int*                  tail,                /**< array of edge tail nodes */
   int*                  head,                /**< array of edge head nodes */
   int**                 adj_lists,          /**< adjacency lists for each node */
   int*                  adj_sizes,          /**< size of adjacency list for each node */
   SCIP_Bool**           adj_matrix,         /**< adjacency matrix */
   int*                  node_degree,       /**< degree of each node */
   int                   k,                  /**< parameter k for k-vertex cut */

   SCIP_VAR**            x_vars,             /**< array of node variables */
   SCIP_VAR**            z_v_vars,           /**< array of z_v variables */
   SCIP_VAR**            z_u_vars,           /**< array of z_u variables */
   SCIP_VAR**            beta_v_vars,        /**< array of beta_v variables */
   SCIP_VAR**            beta_u_vars,        /**< array of beta_u variables */

   SCIP_VAR**            alphavars,               /**< array of variables */
   int                   n_alphavars,             /**< number of variables */
   int                   alphavars_size,          /**< size of variables array */

   SCIP_CONS*            main_alpha_constr,          /**< main pricing constraint for α_S */
   SCIP_CONS**           alpha_constrs,              /**< array of pricing constraints for α_S */
   SCIP_CONS**           coverage_constrs,         /**< array of SCC coverage constraints */
   int                   nalphaconstrs,             /**< number of pricing constraints for α_S */

   SCIP_CONS**           other_contrs,        /**< array of other linear constraints*/
   int                   notherconstrs        /**< number of other linear constraints */
   )
{
   int i;
   
   assert(scip != NULL);
   assert(targetdata != NULL);

   /* allocate memory for problem data */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );

   /* copy basic problem parameters */
   (*targetdata)->nnodes = nnodes;
   (*targetdata)->nedges = nedges;
   (*targetdata)->k = k;
   (*targetdata)->n_alphavars = n_alphavars;
   (*targetdata)->alphavars_size = alphavars_size;
   (*targetdata)->nalphaconstrs = nalphaconstrs;
   (*targetdata)->notherconstrs = notherconstrs;

   /* duplicate edge arrays */
   if( nedges > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->tail, tail, nedges) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->head, head, nedges) );
   }
   else
   {
      (*targetdata)->tail = NULL;
      (*targetdata)->head = NULL;
   }

   /* duplicate node arrays */
   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->node_degree, node_degree, nnodes) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->adj_sizes, adj_sizes, nnodes) );
   }
   else
   {
      (*targetdata)->node_degree = NULL;
      (*targetdata)->adj_sizes = NULL;
   }

   /* duplicate adjacency matrix */
   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*targetdata)->adj_matrix, nnodes) );
      for( i = 0; i < nnodes; ++i )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->adj_matrix[i], adj_matrix[i], nnodes) );
      }
   }
   else
   {
      (*targetdata)->adj_matrix = NULL;
   }

   /* duplicate adjacency lists */
   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*targetdata)->adj_lists, nnodes) );
      for( i = 0; i < nnodes; ++i )
      {
         if( adj_sizes[i] > 0 )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->adj_lists[i], adj_lists[i], adj_sizes[i]) );
         }
         else
         {
            (*targetdata)->adj_lists[i] = NULL;
         }
      }
   }
   else
   {
      (*targetdata)->adj_lists = NULL;
   }

   /* duplicate variable arrays */
   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->x_vars, x_vars, nnodes) );
   }
   else
   {
      (*targetdata)->x_vars = NULL;
   }

   if( nedges > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->beta_u_vars, beta_u_vars, nedges) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->beta_v_vars, beta_v_vars, nedges) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->z_u_vars, z_u_vars, nedges) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->z_v_vars, z_v_vars, nedges) );
   }
   else
   {
      (*targetdata)->beta_u_vars = NULL;
      (*targetdata)->beta_v_vars = NULL;
      (*targetdata)->z_u_vars = NULL;
      (*targetdata)->z_v_vars = NULL;
   }

   /* duplicate alpha variables array */
   if( n_alphavars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->alphavars, alphavars, alphavars_size) );
   }
   else
   {
      (*targetdata)->alphavars = NULL;
   }

   /* copy constraint pointers - COME NEL BPP */
   (*targetdata)->main_alpha_constr = main_alpha_constr;

   if( nalphaconstrs > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->alpha_constrs, alpha_constrs, nalphaconstrs) );
   }
   else
   {
      (*targetdata)->alpha_constrs = NULL;
   }

   if( notherconstrs > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->other_contrs, other_contrs, notherconstrs) );
   }
   else
   {
      (*targetdata)->other_contrs = NULL;
   }

   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->coverage_constrs, coverage_constrs, nnodes) );
   }
   else
   {
      (*targetdata)->coverage_constrs = NULL;
   }

   return SCIP_OKAY;
}

/**@} */


/** create initial columns */
static
SCIP_RETCODE createInitialColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   
   SCIP_CONS** coverage_constrs;

   int nnodes = probdata->nnodes;
   int nedges = probdata->nedges;

   coverage_constrs = probdata->coverage_constrs;
   
   SCIP_VAR* newvar;
   SCIP_VARDATA* vardata;
   
   int subsetsize;
   int* subset = new int[nnodes];
   
   char name[SCIP_MAXSTRLEN];
   int e, u, v;

   for(v = 0; v < nnodes; ++v)
   {
      subsetsize = 1;   

      for( int ind = 0; ind < nnodes; ++ind )
         subset[ind] = 0;

      subset[v] = 1;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "vertices_%d", v);

      SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, subsetsize) );

      /* create variable for a new column with objective function coefficient 0.0 */
      SCIP_CALL( SCIPcreateVarKvertexcut(scip, &newvar, name, 0.0, FALSE, FALSE, vardata));

      /* add the new variable to the pricer store */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, newvar) );

      SCIP_CALL( SCIPaddCoefLinear(scip, coverage_constrs[v], newvar, 1.0) );

      SCIPdebug(SCIPprintVar(scip, newvar, NULL) );
      SCIP_CALL(SCIPreleaseVar(scip, &newvar) );
   }

   for(e = 0; e < nedges; ++e)
   {
      u = probdata->tail[e] - 1;
      v = probdata->head[e] - 1;

      subsetsize = 2;   

      for( int ind = 0; ind < nnodes; ++ind )
         subset[ind] = 0;

      subset[u] = 1;
      subset[v] = 1;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "vertices_%d_%d", u, v);

      SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, subsetsize) );

      /* create variable for a new column with objective function coefficient 0.0 */
      SCIP_CALL( SCIPcreateVarKvertexcut(scip, &newvar, name, 0.0, FALSE, FALSE, vardata));

      /* add the new variable to the pricer store */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, newvar) );

      SCIP_CALL( SCIPaddCoefLinear(scip, coverage_constrs[u], newvar, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, coverage_constrs[v], newvar, 1.0) );
      
      SCIP_CALL( SCIPaddCoefLinear(scip, probdata->alpha_constrs[e], newvar, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, probdata->main_alpha_constr, newvar, 1.0) );

      SCIPdebug(SCIPprintVar(scip, newvar, NULL) );
      SCIP_CALL(SCIPreleaseVar(scip, &newvar) );
   }

   delete [] subset;

   return SCIP_OKAY;
}





/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigKvertexcut)
{
   SCIPdebugMsg(scip, "free original k-vertex cut problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data */
static
SCIP_DECL_PROBTRANS(probtransKvertexcut)
{

   SCIP_CALL( probdataTrans(scip, targetdata,
         sourcedata->nnodes, sourcedata->nedges, sourcedata->tail, sourcedata->head, sourcedata->adj_lists,
         sourcedata->adj_sizes, sourcedata->adj_matrix, sourcedata->node_degree, sourcedata->k,
         sourcedata->x_vars, sourcedata->z_v_vars, sourcedata->z_u_vars, sourcedata->beta_v_vars, sourcedata->beta_u_vars,
         sourcedata->alphavars, sourcedata->n_alphavars, sourcedata->alphavars_size,
         sourcedata->main_alpha_constr, sourcedata->alpha_constrs, sourcedata->coverage_constrs, sourcedata->nalphaconstrs,
         sourcedata->other_contrs, sourcedata->notherconstrs) );

   /* transform all variables separately */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nnodes, (*targetdata)->x_vars, (*targetdata)->x_vars) );

   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nedges, (*targetdata)->beta_u_vars, (*targetdata)->beta_u_vars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nedges, (*targetdata)->beta_v_vars, (*targetdata)->beta_v_vars) );
   
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nedges, (*targetdata)->z_u_vars, (*targetdata)->z_u_vars) );
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nedges, (*targetdata)->z_v_vars, (*targetdata)->z_v_vars) );

   if( (*targetdata)->n_alphavars > 0 )
   {
      SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->n_alphavars, (*targetdata)->alphavars, (*targetdata)->alphavars) );
   }

   /* transform all constraints separately */
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->main_alpha_constr, &(*targetdata)->main_alpha_constr) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nalphaconstrs, (*targetdata)->alpha_constrs, (*targetdata)->alpha_constrs) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nnodes, (*targetdata)->coverage_constrs, (*targetdata)->coverage_constrs) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->notherconstrs, (*targetdata)->other_contrs, (*targetdata)->other_contrs) );
   
   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransKvertexcut)
{
   SCIPdebugMsg(scip, "free transformed k-vertex cut problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );
   
   return SCIP_OKAY;
}

/** solving process initialization method of transformed data */
static
SCIP_DECL_PROBINITSOL(probinitsolKvertexcut)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   /* catch node focused event */
   eventhdlr = SCIPfindEventhdlr(scip, "DENTRO_NODO");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );


   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data */
static
SCIP_DECL_PROBEXITSOL(probexitsolKvertexcut)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

   /* drop node focused event */
   eventhdlr = SCIPfindEventhdlr(scip, "DENTRO_NODO");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   int                   nnodes,             /**< number of nodes in the graph */
   int                   nedges,             /**< number of edges in the graph */
   int*                  tail,              /**< array of edge tail nodes */
   int*                  head,              /**< array of edge head nodes */
   int                   k                   /**< parameter k for k-vertex cut */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** alpha_constrs;
   SCIP_CONS** coverage_constrs;

   char name[SCIP_MAXSTRLEN];
   int i, j;

   assert(scip != NULL);
   assert(nnodes > 0);
   assert(nedges >= 0);
   assert(tail != NULL || nedges == 0);
   assert(head != NULL || nedges == 0);

   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, "DENTRO_NODO", "Does stuff inside a node", eventExecMyLogic, NULL) );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigKvertexcut) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransKvertexcut) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransKvertexcut) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolKvertexcut) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolKvertexcut) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   /* allocate memory for problem data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &probdata) );


   /* First create variables for each node (x_v) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->x_vars, nnodes) );

   for( i = 0; i < nnodes; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      
      /* Create binary variable for each node (weight = 1 for unweighted case) */
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->x_vars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->x_vars[i]) );

      // SCIP_CALL(SCIPchgVarBranchPriority(scip, probdata->x_vars[i], 1));
   }

   /* Constraint (13): Main alpha constraint */
   /* expr: ∑_S (|S| - 1)α_S + ∑_V x_v + ∑_E (β^u_{uv} + β^v_{uv} - z^u_{uv} - z^v_{uv}) ≤ n - k */
   SCIP_CONS* main_alpha_constr;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &main_alpha_constr, "main_alpha_constraint", 0, NULL, NULL, -SCIPinfinity(scip), (SCIP_Real)(nnodes - k)) );

   /* Add +∑x_v terms to main constraint */
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, probdata->x_vars[i], 1.0) );
   }

   /* This constraint will be modifiable for adding α_S variables during pricing */
   SCIP_CALL( SCIPsetConsModifiable(scip, main_alpha_constr, TRUE) );
   SCIP_CALL( SCIPaddCons(scip, main_alpha_constr) );

   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->beta_u_vars, nedges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->beta_v_vars, nedges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->z_u_vars, nedges) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->z_v_vars, nedges) );

   /* Alloc memory for storing constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &alpha_constrs, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coverage_constrs, nnodes) );
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->other_contrs, nedges * 6) );
   probdata->notherconstrs = 0;

   for( i = 0; i < nedges; ++i )
   {
      int u = tail[i] - 1; /* convert to 0-indexed */
      int v = head[i] - 1;
      
      /* Create β^u_{uv} and β^v_{uv} variables */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "beta_u_%d_%d", u, v);
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->beta_u_vars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->beta_u_vars[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "beta_v_%d_%d", u, v);
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->beta_v_vars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->beta_v_vars[i]) );

      /* Create z^u_{uv} and z^v_{uv} variables for linearization */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "z_u_%d_%d", u, v);
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->z_u_vars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->z_u_vars[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "z_v_%d_%d", u, v);
      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->z_v_vars[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, probdata->z_v_vars[i]) );

      /* Constraint: ∑_S α_S + β^u_{uv} + β^v_{uv} ≥ 1  */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "second_alpha_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &alpha_constrs[i], name, 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
      
      SCIP_CALL( SCIPaddCoefLinear(scip, alpha_constrs[i], probdata->beta_u_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, alpha_constrs[i], probdata->beta_v_vars[i], 1.0) );
      
      /* This constraint will be modifiable for adding α_S variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, alpha_constrs[i], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, alpha_constrs[i]) );

      /* Add β variables to main constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, probdata->beta_u_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, probdata->beta_v_vars[i], 1.0) );

      /* Linearization constraints (15)-(18) for z^u_{uv} = β^u_{uv} * x_u */
      
      /* Constraint: z^u_{uv} ≤ x_u */
      SCIP_CONS* lin1_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin1_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin1_cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin1_cons, probdata->z_u_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin1_cons, probdata->x_vars[u], -1.0) ); 
      SCIP_CALL( SCIPaddCons(scip, lin1_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin1_cons;

      /* Constraint: z^u_{uv} ≤ β^u_{uv} */
      SCIP_CONS* lin2_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin2_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin2_cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin2_cons, probdata->z_u_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin2_cons, probdata->beta_u_vars[i], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, lin2_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin2_cons;

      /* Constraint: z^u_{uv} ≥ β^u_{uv} - 1 + x_u */
      SCIP_CONS* lin3_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin3_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin3_cons, name, 0, NULL, NULL, -1.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3_cons, probdata->z_u_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3_cons, probdata->beta_u_vars[i], -1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3_cons, probdata->x_vars[u], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, lin3_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin3_cons;

      /* Same linearization constraints for z^v_{uv} = β^v_{uv} * x_v */
      
      /* z^v_{uv} ≤ x_v */
      SCIP_CONS* lin1v_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin1v_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin1v_cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin1v_cons, probdata->z_v_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin1v_cons, probdata->x_vars[v], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, lin1v_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin1v_cons;

      /* z^v_{uv} ≤ β^v_{uv} */
      SCIP_CONS* lin2v_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin2v_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin2v_cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin2v_cons, probdata->z_v_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin2v_cons, probdata->beta_v_vars[i], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, lin2v_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin2v_cons;

      /* z^v_{uv} ≥ β^v_{uv} - 1 + x_v */
      SCIP_CONS* lin3v_cons;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin3v_%d_%d", u, v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lin3v_cons, name, 0, NULL, NULL, -1.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3v_cons, probdata->z_v_vars[i], 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3v_cons, probdata->beta_v_vars[i], -1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lin3v_cons, probdata->x_vars[v], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, lin3v_cons) );
      probdata->other_contrs[probdata->notherconstrs++] = lin3v_cons;

      /* Add z variables to main constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, probdata->z_u_vars[i], -1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, probdata->z_v_vars[i], -1.0) );
   }

   for(int v = 0; v < nnodes; ++v )
   {
      /* Constraint: ∑_S α_S + x_v = 0  */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "coverage_constr_%d", v);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &coverage_constrs[v], name, 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );

      SCIP_CALL( SCIPaddCoefLinear(scip, coverage_constrs[v], probdata->x_vars[v], 1.0) );

      /* This constraint will be modifiable for adding α_S variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, coverage_constrs[v], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, coverage_constrs[v]) );
   }

     
   /* store basic graph information */
   probdata->nnodes = nnodes;
   probdata->nedges = nedges;
   probdata->k = k;
   probdata->main_alpha_constr = main_alpha_constr;

   /* allocate and copy edge arrays */
   if( nedges > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->tail, nedges) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->head, nedges) );
      
      for( i = 0; i < nedges; ++i )
      {
         probdata->tail[i] = tail[i];
         probdata->head[i] = head[i];
      }
   }
   else
   {
      probdata->tail = NULL;
      probdata->head = NULL;
   }

   /* allocate memory for adjacency structures */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->adj_lists, nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->adj_sizes, nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->node_degree, nnodes) );

   /* allocate memory for adjacency matrix */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->adj_matrix, nnodes) );
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->adj_matrix[i], nnodes) );
   }

   /* initialize adjacency matrix to FALSE */
   for( i = 0; i < nnodes; ++i )
   {
      for( j = 0; j < nnodes; ++j )
      {
         probdata->adj_matrix[i][j] = FALSE;
      }
   }

   /* initialize node degrees to 0 */
   for( i = 0; i < nnodes; ++i )
   {
      probdata->node_degree[i] = 0;
   }

   /* count degrees and fill adjacency matrix */
   for( i = 0; i < nedges; ++i )
   {
      int u = tail[i];
      int v = head[i];
      
      /* assuming nodes are 1-indexed, convert to 0-indexed */
      if( u >= 1 && v >= 1 && u <= nnodes && v <= nnodes )
      {
         u--; /* convert to 0-indexed */
         v--; /* convert to 0-indexed */
         
         probdata->adj_matrix[u][v] = TRUE;
         probdata->adj_matrix[v][u] = TRUE; /* undirected graph */

         probdata->node_degree[u]++;
         probdata->node_degree[v]++;
      }
   }

   /* allocate and build adjacency lists */
   for( i = 0; i < nnodes; ++i )
   {
      probdata->adj_sizes[i] = probdata->node_degree[i];
      if( probdata->node_degree[i] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->adj_lists[i], probdata->node_degree[i]) );
      }
      else
      {
         probdata->adj_lists[i] = NULL;
      }
   }

   /* fill adjacency lists */
   int* counters;
   SCIP_CALL( SCIPallocBufferArray(scip, &counters, nnodes) );
   for( i = 0; i < nnodes; ++i )
      counters[i] = 0;

   for( i = 0; i < nedges; ++i )
   {
      int u = tail[i];
      int v = head[i];
      
      if( u >= 1 && v >= 1 && u <= nnodes && v <= nnodes )
      {
         u--; /* convert to 0-indexed */
         v--; /* convert to 0-indexed */
         
         probdata->adj_lists[u][counters[u]++] = v;
         probdata->adj_lists[v][counters[v]++] = u;
      }
   }
   
   SCIPfreeBufferArray(scip, &counters);

   /* store constraints */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &probdata->alpha_constrs, alpha_constrs, nedges) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &probdata->coverage_constrs, coverage_constrs, nnodes) );
   probdata->nalphaconstrs = nedges;

   /* initialize variables arrays */
   probdata->alphavars = NULL;
   probdata->n_alphavars = 0;
   probdata->alphavars_size = 0; 
   
   /* initialize singleton columns*/
   SCIP_CALL( createInitialColumns(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* activate pricer for k-vertex cut */
   SCIP_CALL( SCIPpricerKvertexcutActivate(scip, probdata->main_alpha_constr, probdata->alpha_constrs, probdata->coverage_constrs, nnodes, nedges) );

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &alpha_constrs);

   SCIPinfoMessage(scip, NULL, "Created k-vertex cut problem: %d nodes, %d edges, k=%d\n", nnodes, nedges, k);

   /* write problem formulation to .lp file for debugging */
   SCIP_CALL( SCIPwriteOrigProblem(scip, "new_formulation.lp", "lp", FALSE) );

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nedges; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->beta_u_vars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->beta_v_vars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->z_u_vars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->z_v_vars[i]) );
   }

   for( i = 0; i < (*probdata)->nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->x_vars[i]) );
   }
   
   for( i = 0; i < (*probdata)->n_alphavars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->alphavars[i]) );
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->nalphaconstrs; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->alpha_constrs[i]) );
   }

   for( i = 0; i < (*probdata)->nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->coverage_constrs[i]) );
   }

   for( i = 0; i < (*probdata)->notherconstrs; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->other_contrs[i]) );
   }

   /* release main constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->main_alpha_constr) );
   

   /* free memory of adjacency lists */
   for( i = 0; i < (*probdata)->nnodes; ++i )
   {
      if( (*probdata)->adj_lists[i] != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &(*probdata)->adj_lists[i], (*probdata)->adj_sizes[i]);
      }
   }

   /* free memory of adjacency matrix */
   for( i = 0; i < (*probdata)->nnodes; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->adj_matrix[i], (*probdata)->nnodes);
   }

   /* free memory of arrays */
   if( (*probdata)->n_alphavars > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->alphavars, (*probdata)->alphavars_size);
   }

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->x_vars, (*probdata)->nnodes);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->beta_u_vars, (*probdata)->nedges);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->beta_v_vars, (*probdata)->nedges);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->z_u_vars, (*probdata)->nedges);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->z_v_vars, (*probdata)->nedges);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->alpha_constrs, (*probdata)->nalphaconstrs);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->coverage_constrs, (*probdata)->nnodes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->other_contrs, (*probdata)->notherconstrs);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->adj_matrix, (*probdata)->nnodes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->adj_lists, (*probdata)->nnodes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->adj_sizes, (*probdata)->nnodes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->node_degree, (*probdata)->nnodes);
   
   if( (*probdata)->nedges > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->tail, (*probdata)->nedges);
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->head, (*probdata)->nedges);
   }

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/**@} */



/**@name Getter functions for k-vertex cut problem data
 *
 * @{
 */

  /** returns number of variables */
int SCIPprobdataGetNAlphaVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->n_alphavars;
}


/** returns array of all variables ordered in the way they got generated */
SCIP_VAR** SCIPprobdataGetAlphaVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->alphavars;
}

SCIP_VAR** SCIPprobdataGetXVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->x_vars;
}

SCIP_VAR** SCIPprobdataGetZUVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->z_u_vars;
}

SCIP_VAR** SCIPprobdataGetZVVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->z_v_vars;
}

SCIP_VAR** SCIPprobdataGetBetaUVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->beta_u_vars;
}


SCIP_VAR** SCIPprobdataGetBetaVVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->beta_v_vars;
}


/** returns number of nodes */
int SCIPprobdataGetNNodes(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->nnodes;
}

/** returns number of edges */
int SCIPprobdataGetNEdges(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->nedges;
}

/** returns array of edge tail nodes */
int* SCIPprobdataGetEdgeTails(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->tail;
}

/** returns array of edge head nodes */
int* SCIPprobdataGetEdgeHeads(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->head;
}

/** returns adjacency list for a given node */
int* SCIPprobdataGetAdjList(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(node >= 0 && node < probdata->nnodes);
   return probdata->adj_lists[node];
}

/** returns size of adjacency list for a given node */
int SCIPprobdataGetAdjSize(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(node >= 0 && node < probdata->nnodes);
   return probdata->adj_sizes[node];
}

/** returns degree of a given node */
int SCIPprobdataGetNodeDegree(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(node >= 0 && node < probdata->nnodes);
   return probdata->node_degree[node];
}

/** returns parameter k */
int SCIPprobdataGetK(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->k;
}

/** checks if two nodes are adjacent */
SCIP_Bool SCIPprobdataAreAdjacent(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node1,              /**< first node */
   int                   node2               /**< second node */
   )
{
   assert(probdata != NULL);
   assert(node1 >= 0 && node1 < probdata->nnodes);
   assert(node2 >= 0 && node2 < probdata->nnodes);
   return probdata->adj_matrix[node1][node2];
}

/** returns adjacency matrix */
SCIP_Bool** SCIPprobdataGetAdjMatrix(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->adj_matrix;
}

/** returns main alpha constraint */
SCIP_CONS* SCIPprobdataGetMainAlphaConstr(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->main_alpha_constr;
}

/** returns |E|-set of alpha constraints */
SCIP_CONS** SCIPprobdataGetAlphaConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->alpha_constrs;
}

/** returns |V|-set of coverage constraints */
SCIP_CONS** SCIPprobdataGetCoverageConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->coverage_constrs;
}


/**@} */
/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   /* check if enough memory is left */
   if( probdata->alphavars_size == probdata->n_alphavars )
   {
      int newsize;
      newsize = MAX(100, probdata->alphavars_size * 2);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->alphavars, probdata->alphavars_size, newsize) );
      probdata->alphavars_size = newsize;
   }

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   probdata->alphavars[probdata->n_alphavars] = var;
   probdata->n_alphavars++;

   return SCIP_OKAY;
}

/**@} */
