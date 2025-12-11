/**@file   probdata_kvertexcut.cpp
 * @brief  Problem data for k-vertex cut problem
 * @author Fabio Ciccarelli
 *
 * This file handles the main problem data used in that project. 
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
#include <stack>
#include <vector>
#include <algorithm>

#include <queue>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>


#include "gurobi_c++.h"
#include "probdata_kvertexcut.h"
#include "vardata_kvertexcut.h"
#include "pricer_kvertexcut.h"
#include "custom_compute_symmetry.h"
#include "graph.h"
#include "preprocessing.h"

#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip/symmetry_graph.h"
#include "scip/cons_symresack.h"
#include "scip/event_shadowtree.h"
#include "scip/symmetry_orbital.h"
#include "scip/symmetry_lexred.h"

using namespace std;

/** forward declaration */
static SCIP_RETCODE probdataFree(SCIP* scip, SCIP_PROBDATA** probdata);


/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the k-vertex-cut, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   Graph*                orig_graph;              /**< original graph structure */
   Graph*                res_graph;               /**< residual graph structure */
   int                   k;                      /**< parameter k for k-vertex cut */

   SCIP_VAR**            x_vars;             /**< array of node variables */

   SCIP_VAR**            alphavars;               /**< array of variables */
   int                   n_alphavars;             /**< number of variables */
   int                   alphavars_size;          /**< size of variables array */

   SCIP_ROW**           clique_cuts;            /**< additional clique cuts */
   int                   n_clique_cuts;          /**< number of additional clique cuts */
   int                   clique_cuts_size;      /**< size of clique cuts array */
   std::vector<int*>     clique_cuts_cliques;   /**< cliques associated with the clique cuts */

   SCIP_CONS**           symconss;           /**< symmetry handling constraints */
   int                   nsymconss;          /**< number of symmetry handling constraints */
   int**                 perms;              /**< permutations */
   int                   nperms;             /**< number of permutations */
   int                   lenperms;           /**< length of permutations */
   SCIP_ORBITALREDDATA*  orbitalreddata;     /**< container for the orbital reduction data */
   SCIP_LEXREDDATA*      lexreddata;         /**< container for the lexcographic reduction data */

   SCIP_CONS*            alpha_cardinality_constr;     /**< cardinality constraint for α variables */
   SCIP_CONS*            min_connectivity_constr;      /**< minimum connectivity constraint: sum(x) >= min_connectivity */
   SCIP_CONS**           vertex_cover_constrs;         /**< array of SCC coverage constraints */
   SCIP_CONS**           clique_constrs;               /**< array of clique constraints */
   int                   n_cliques;                    /**< number of clique constraints */
   int                   min_connectivity;             /**< minimum vertex connectivity of the residual graph */

   int                   root_generated_vars;              /**< number of variables generated at root */
   int                   farkas_generated_vars;            /**< number of generated variables with farkas pricing */
   double                root_lower_bound;                 /**< value of the LP relaxation */
   double                root_upper_bound;                 /**< best incumbent value at the root node */
   
   /* Preprocessing data */
   bool*                 preFixed;                /**< array indicating which nodes are fixed by preprocessing */
   int                   n_fixed;                 /**< number of fixed nodes */
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

/** creates transformed problem data */
static
SCIP_RETCODE probdataTrans(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       targetdata,         /**< pointer to target problem data */
   Graph*                orig_graph,         /**< original graph structure */
   Graph*                res_graph,          /**< residual graph structure */
   int                   k,                  /**< parameter k for k-vertex cut */

   SCIP_VAR**            x_vars,             /**< array of node variables */

   SCIP_VAR**            alphavars,               /**< array of variables */
   int                   n_alphavars,             /**< number of variables */
   int                   alphavars_size,          /**< size of variables array */

   SCIP_ROW**           clique_cuts,            /**< additional clique cuts */
   int                   n_clique_cuts,          /**< number of additional clique cuts */
   int                   clique_cuts_size,       /**< size of clique cuts array */
   std::vector<int*>     clique_cuts_cliques,   /**< cliques associated with the clique cuts */

   SCIP_CONS**           symconss,           /**< symmetry handling constraints */
   int                   nsymconss,          /**< number of symmetry handling constraints */
   SCIP_ORBITALREDDATA*  orbitalreddata,     /**< orbital reduction data */
   SCIP_LEXREDDATA*      lexreddata,         /**< lexicographic reduction data */

   SCIP_CONS*            alpha_cardinality_constr,        /**< cardinality constraint for α variables */
   SCIP_CONS*            min_connectivity_constr,         /**< minimum connectivity constraint */
   SCIP_CONS**           vertex_cover_constrs,            /**< array of SCC coverage constraints */
   SCIP_CONS**           clique_constrs,                  /**< array of clique constraints */
   int                   n_cliques,                       /**< number of clique constraints */
   int                   min_connectivity,                /**< minimum vertex connectivity */
   
   bool*                 preFixed,            /**< array of preprocessing fixed nodes */
   int                   n_fixed              /**< number of fixed nodes */
   )
{
   int nnodes = res_graph->nnodes;
   
   assert(scip != NULL);
   assert(targetdata != NULL);

   /* allocate memory for problem data */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );

   /* copy Graph structure (deep copy) */
   (*targetdata)->orig_graph = new Graph(*orig_graph);
   (*targetdata)->res_graph = new Graph(*res_graph);

   /* copy basic problem parameters */
   (*targetdata)->k = k;
   (*targetdata)->n_alphavars = n_alphavars;
   (*targetdata)->alphavars_size = alphavars_size;
   (*targetdata)->n_clique_cuts = n_clique_cuts;
   (*targetdata)->clique_cuts_size = clique_cuts_size;
   (*targetdata)->clique_cuts_cliques = clique_cuts_cliques;
   (*targetdata)->n_cliques = n_cliques;
   (*targetdata)->min_connectivity = min_connectivity;
   (*targetdata)->nsymconss = nsymconss;
   (*targetdata)->orbitalreddata = orbitalreddata;
   (*targetdata)->lexreddata = lexreddata;

   /* duplicate variable arrays */
   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->x_vars, x_vars, nnodes) );
   }
   else
   {
      (*targetdata)->x_vars = NULL;
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

   /* duplicate alpha variables array */
   if( n_clique_cuts > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->clique_cuts, clique_cuts, clique_cuts_size) );
   }
   else
   {
      (*targetdata)->clique_cuts = NULL;
   }

   /* copy constraint pointers */
   (*targetdata)->alpha_cardinality_constr = alpha_cardinality_constr;
   (*targetdata)->min_connectivity_constr = min_connectivity_constr;


   if( nnodes > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->vertex_cover_constrs, vertex_cover_constrs, nnodes) );
   }
   else
   {
      (*targetdata)->vertex_cover_constrs = NULL;
   }

   if( n_cliques > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->clique_constrs, clique_constrs, n_cliques) );
   }
   else
   {
      (*targetdata)->clique_constrs = NULL;
   }

   /* copy preprocessing data */
   (*targetdata)->n_fixed = n_fixed;
   if( orig_graph->nnodes > 0 && preFixed != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*targetdata)->preFixed, orig_graph->nnodes) );
      memcpy((*targetdata)->preFixed, preFixed, orig_graph->nnodes * sizeof(bool));
   }
   else
   {
      (*targetdata)->preFixed = NULL;
   }

   if( nsymconss > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->symconss, symconss, nsymconss) );
   }
   else
      (*targetdata)->symconss = NULL;

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
   
   SCIP_CONS* alpha_cardinality_constr;
   SCIP_CONS** vertex_cover_constrs;
   SCIP_CONS** clique_constrs;

   int nnodes = probdata->res_graph->nnodes;
   int n_cliques = probdata->n_cliques;
   const std::vector<std::vector<bool>>* cliques = probdata->res_graph->getCliques();

   alpha_cardinality_constr = probdata->alpha_cardinality_constr;
   vertex_cover_constrs = probdata->vertex_cover_constrs;
   clique_constrs = probdata->clique_constrs;
   
   SCIP_VAR* newvar;
   SCIP_VARDATA* vardata;
   
   int subsetsize;
   int* subset = new int[nnodes];
   
   char name[SCIP_MAXSTRLEN];

   for(int v = 0; v < nnodes; ++v)
   {
      subsetsize = 1;   

      for( int ind = 0; ind < nnodes; ++ind )
         subset[ind] = 0;

      subset[v] = 1;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "vertices_%d", v);

      SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, subsetsize, nnodes) );

      /* create variable for a new column with objective function coefficient 0.0 */
      SCIP_CALL( SCIPcreateVarKvertexcut(scip, &newvar, name, 0.0, FALSE, FALSE, vardata));

      /* add variable to the problem */
      SCIP_CALL( SCIPaddVar(scip, newvar) );
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, newvar) );

      SCIP_CALL( SCIPaddCoefLinear(scip, vertex_cover_constrs[v], newvar, 1.0) );

      for( int c = 0; c < n_cliques; ++c )
      {
         if( (*cliques)[c][v] == true )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, clique_constrs[c], newvar, 1.0) );
         }
      }

      SCIP_CALL( SCIPaddCoefLinear(scip, alpha_cardinality_constr, newvar, 1.0) );

      SCIPdebug(SCIPprintVar(scip, newvar, NULL) );

      SCIP_CALL(SCIPreleaseVar(scip, &newvar) );
   }

   delete [] subset;

   return SCIP_OKAY;
}



/** Warm start heuristic */
static
SCIP_RETCODE warmStartHeuristic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,            /**< problem data */
   bool                  ILP               /**< whether to use a ILP model */
   )
{  /*lint --e{715}*/
   SCIP_SOL* sol;
   SCIP_Bool success;
   
   assert(scip != NULL);
   assert(probdata != NULL);   
   
   /* get graph, k, and x variables using getter functions */
   Graph* res_graph = probdata->res_graph;
   int k = probdata->k;
   SCIP_VAR** xvars = probdata->x_vars;
   int n_nodes = res_graph->nnodes;
   
   /* allocate solution array */
   bool* heu_sol = new bool[n_nodes];
   for( int i = 0; i < n_nodes; i++ )
      heu_sol[i] = false;

   int heu_sol_weight = 0;

   if(!ILP)
   {
      /* compute iterative separator solution */
      heu_sol_weight = computeIterativeSeparatorSolution(*res_graph, k, heu_sol);
   }
   else{
      heu_sol_weight = ILPwarmStart(*res_graph, k, heu_sol);
   }
   
   
   // cout << "Separator size computed: " << heu_sol_weight << endl;
   
   if( heu_sol_weight == 0 )
   {
      cout << "Could not compute separator solution (size=0)" << endl;
      SCIPdebugMsg(scip, "Could not compute separator solution\n");
      delete[] heu_sol;
      return SCIP_OKAY;
   }
   
   /* compute connected components in the residual graph (after removing separator) */
   vector<bool> visited(n_nodes, false);
   vector<vector<int>> components;  // each component is a list of nodes
   
   for( int start = 0; start < n_nodes; start++ )
   {
      // Skip nodes in separator or already visited
      if( heu_sol[start] || visited[start] )
         continue;
      
      // Found a new component - do DFS
      vector<int> component;
      stack<int> st;
      st.push(start);
      visited[start] = true;
      
      while( !st.empty() )
      {
         int node = st.top();
         st.pop();
         component.push_back(node);
         
         // Visit all neighbors (that are not in separator)
         int num_neighbors;
         const int* neighbors = res_graph->getNeighbors(node, num_neighbors);
         for( int i = 0; i < num_neighbors; i++ )
         {
            int neighbor = neighbors[i];
            if( !heu_sol[neighbor] && !visited[neighbor] )
            {
               visited[neighbor] = true;
               st.push(neighbor);
            }
         }
      }
      
      components.push_back(component);
   }


   if( (int) components.size() < k)
   {
      // Graph is already connected
      SCIPdebugMsg(scip, "Separator solution does not disconnect the graph in enough components...\n");
      delete[] heu_sol;
      return SCIP_OKAY;
   }
   
   
   /* create solution structure */
   SCIP_CALL(SCIPcreateSol(scip, &sol, NULL));
   
   /* set values for x variables (node removal variables) */
   for( int i = 0; i < n_nodes; i++ )
   {
      SCIP_VAR* x_var = xvars[i];
      SCIP_Real val = heu_sol[i] ? 1.0 : 0.0;
      
      if( x_var != NULL )
      {
         SCIP_CALL(SCIPsetSolVal(scip, sol, x_var, val));
      }
   }
   
   /* Now create and set alpha variables for each connected component */
 
   for( size_t comp_idx = 0; comp_idx < components.size(); comp_idx++ )
   {
      vector<int>& component = components[comp_idx];
      int comp_size = component.size();
      
      // Sort the component for consistent naming
      sort(component.begin(), component.end());
      
      // Look for existing alpha variable for this subset
      SCIP_VAR** alphavars = SCIPprobdataGetAlphaVars(probdata);
      SCIP_VAR* alpha_var = NULL;
      bool found_existing = false;
      
      if(comp_size == 1)
      {
         // Singleton
         int node = component[0];
         alpha_var = alphavars[node];
         found_existing = true;
      }
      
      // If not found, create new alpha variable (ONLY for non-singletons!)
      if( !found_existing )
      {
      
         string name = "vertices";
      
         // Create vardata
         int* subset = new int[n_nodes];
         memset(subset, 0, n_nodes * sizeof(int));
         for( int j = 0; j < comp_size; j++ )
         {
            subset[component[j]] = 1;
            name += "_" + std::to_string(component[j]);
         }

         
         SCIP_VARDATA* vardata;
         SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, comp_size, n_nodes) );
         
         // Create variable (obj = 0, not initial, removable)
         SCIP_CALL( SCIPcreateVarKvertexcut(scip, &alpha_var, name.c_str(), 0.0, FALSE, TRUE, vardata) );
         
         // Add variable to problem
         SCIP_CALL( SCIPaddVar(scip, alpha_var) );

         int n_cliques = SCIPprobdataGetNCliques(probdata);
         const vector<vector<bool>>* cliques = res_graph->getCliques();

         SCIP_CONS** vertex_cover_constr = probdata->vertex_cover_constrs;
         SCIP_CONS** clique_constrs = probdata->clique_constrs;
         SCIP_CONS* alpha_cardinality_constr = probdata->alpha_cardinality_constr;


     
         for(int v = 0; v < n_nodes; ++v )
         {
            if( subset[v] == 1 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, vertex_cover_constr[v], alpha_var, 1.0) );
            }
         }

         for (int c = 0; c < n_cliques; c++)
         {
            for(int v = 0; v < n_nodes; v++)
            {
                  if((*cliques)[c][v] == true && subset[v] == 1)
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, clique_constrs[c], alpha_var, 1.0) );
                     break;
                  }
            }
         }

         SCIP_CALL( SCIPaddCoefLinear(scip, alpha_cardinality_constr, alpha_var, 1.0) );

         delete[] subset;
      }
      
      // Set alpha variable to 1 in the solution
      if( alpha_var != NULL )
      {
         if( found_existing )
         {
           SCIP_CALL( SCIPsetSolVal(scip, sol, alpha_var, 1.0) );
         }
         else
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, alpha_var, 1.0) );
            SCIP_CALL( SCIPreleaseVar(scip, &alpha_var) );
         }
      }
   }
   
   /* free the solution array */
   delete[] heu_sol;
      
   /* add solution to SCIP */
   SCIP_CALL( SCIPaddSol(scip, sol, &success) );

   // cout << "\nSCIPaddSol returned, success = " << success << endl;

   if( success )
   {
      // cout << "*** CONNECTIVITY HEURISTIC FOUND FEASIBLE SOLUTION! ***" << endl;
      SCIPdebugMsg(scip, "Connectivity heuristic found feasible solution\n");
   }
   else
   {
      // cout << "Connectivity heuristic solution was NOT accepted" << endl;
      SCIPdebugMsg(scip, "Connectivity heuristic solution was not accepted\n");
   }
   
   /* free solution structure */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
      
   return SCIP_OKAY;
}

/** print information about symmetry handling */
void printSymmetryInfo(
   int                   nperms,             /**< number of symmetry generators */
   int                   symmethod           /**< identifier of symmetry handling method */
   )
{
   printf("\n=== SYMMETRY ===\n");
   printf("number of generators:\t\t%d\n", nperms);
   if( nperms > 0 )
   {
      if( symmethod == 1 )
         printf("handle symmetries by symresacks\n");
      else if( symmethod == 2 )
         printf("handle symmetries by lexicographic reduction and orbital reductions\n");
   }
   printf("\n");
}

/** computes symmetries */
static
SCIP_RETCODE computeSymmetries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            xvars,              /**< x-variables */
   int***                perms,              /**< pointer to hold permutations */
   int*                  nperms,             /**< pointer to hold number of permutations */
   int*                  nmaxperms,          /**< pointer to hold length of perms */
   int                   nnodes,             /**< number of nodes in graph */
   int                   nedges,             /**< number of edges in graph */
   int*                  head,               /**< array storing the head of each egde */
   int*                  tail                /**< array storing the tail of each edge */
   )
{
   SYM_GRAPH* graph;
   SCIP_Real log10groupsize;
   SCIP_Real symcodetime;
   int idx;
   int i;

   assert(scip != NULL);
   assert(perms != NULL);
   assert(nperms != NULL);
   assert(nmaxperms != NULL);
   assert(nnodes > 0);
   assert(nedges > 0);
   assert(head != NULL);
   assert(tail != NULL);

   /* add some variables, we will ignore them after computing the symmetries */
   SCIP_CALL( SCIPcreateSymgraph(scip, SYM_SYMTYPE_PERM, &graph, xvars, 1, 0, nnodes, 0, nedges) );

   Graph* inst_graph = SCIPprobdataGetGraph(SCIPgetProbData(scip));

   /* create nodes of symmetry detection graph */
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPaddSymgraphValnode(scip, graph, inst_graph->vertex_weights[i], &idx) );
   }

   /* add edges */
   for( i = 0; i < nedges; ++i )
   {
      SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, head[i], tail[i], FALSE, 0.0) );
   }

   SCIP_CALL( SCIPcomputeSymgraphColors(scip, graph, 0) );

   SCIP_CALL( SYMcomputeSymmetryGeneratorsNode(scip, 0, graph, nperms, nmaxperms, perms,
         &log10groupsize, &symcodetime) );

   SCIP_CALL( SCIPfreeSymgraph(scip, &graph) );

   return SCIP_OKAY;
}

/** tries to add symmetry handling constraints */
static
SCIP_RETCODE tryAddSymmetryHandlingConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            xvars,              /**< x-variables */
   SCIP_CONS***          symconss,           /**< pointer to array for symmetry handling constraints */
   int*                  nsymconss,          /**< pointer to store number of symmetry handling constraints */
   int                   nnodes,             /**< number of nodes in graph */
   int                   nedges,             /**< number of edges in graph */
   int*                  head,               /**< array storing the head of each egde */
   int*                  tail                /**< array storing the tail of each edge */
   )
{
   char name[SCIP_MAXSTRLEN];
   int i;
   int nperms;
   int nmaxperms;
   int** perms;

   assert(scip != NULL);
   assert(symconss != NULL);

   assert(*symconss == NULL);
   assert(nsymconss != NULL);
   assert(*symconss == 0);
   assert(nnodes > 0);
   assert(nedges > 0);
   assert(head != NULL);
   assert(tail != NULL);

   SCIP_CALL( computeSymmetries(scip, xvars, &perms, &nperms, &nmaxperms, nnodes, nedges, head, tail) );

   printSymmetryInfo(nperms, 1);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, symconss, nperms) );
   for( i = 0; i < nperms; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symcons_%d", i);
      SCIP_CALL( SCIPcreateSymbreakCons(scip, &(*symconss)[i], name, perms[i], xvars, nnodes, FALSE,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, (*symconss)[i]) );
      // cout << "Added symmetry handling constraint " << name << endl;
   }
   *nsymconss = nperms;

   /* free symmetry information */
   for( i = 0; i < nperms; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &perms[i], 1 + nnodes);
   }
   SCIPfreeBlockMemoryArray(scip, &perms, nmaxperms);

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
         sourcedata->orig_graph, sourcedata->res_graph, sourcedata->k,
         sourcedata->x_vars, sourcedata->alphavars, sourcedata->n_alphavars, sourcedata->alphavars_size,
         sourcedata->clique_cuts, sourcedata->n_clique_cuts, sourcedata->clique_cuts_size, sourcedata->clique_cuts_cliques,
         sourcedata->symconss, sourcedata->nsymconss, sourcedata->orbitalreddata, sourcedata->lexreddata,
         sourcedata->alpha_cardinality_constr, sourcedata->min_connectivity_constr, sourcedata->vertex_cover_constrs, sourcedata->clique_constrs, sourcedata->n_cliques, sourcedata->min_connectivity,
         sourcedata->preFixed, sourcedata->n_fixed) );

   /* transform all variables separately */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->res_graph->nnodes, (*targetdata)->x_vars, (*targetdata)->x_vars) );

   if( (*targetdata)->n_alphavars > 0 )
   {
      SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->n_alphavars, (*targetdata)->alphavars, (*targetdata)->alphavars) );
   }

   /* transform all constraints separately */
   if( (*targetdata)->nsymconss > 0 )
   {
      SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nsymconss, (*targetdata)->symconss, (*targetdata)->symconss) );
   }

   /* transform all constraints separately */
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->alpha_cardinality_constr, &(*targetdata)->alpha_cardinality_constr) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->min_connectivity_constr, &(*targetdata)->min_connectivity_constr) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->res_graph->nnodes, (*targetdata)->vertex_cover_constrs, (*targetdata)->vertex_cover_constrs) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->n_cliques, (*targetdata)->clique_constrs, (*targetdata)->clique_constrs) );
   
   /* possibly add dynamic symmetry handling methods */
   if( sourcedata->perms != NULL )
   {
      SCIP_Bool success;
      int p;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*targetdata)->perms, sourcedata->nperms) );
      for( p = 0; p < sourcedata->nperms; ++p )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*targetdata)->perms[p], sourcedata->perms[p],
               sourcedata->res_graph->nnodes) );
      }
      (*targetdata)->nperms = sourcedata->nperms;
      (*targetdata)->lenperms = sourcedata->nperms;

      SCIP_CALL( SCIPorbitalReductionAddComponent(scip, (*targetdata)->orbitalreddata, (*targetdata)->x_vars,
            (*targetdata)->res_graph->nnodes, (*targetdata)->perms, (*targetdata)->nperms, &success) );
      assert(success);
      for( p = 0; p < (*targetdata)->nperms; ++p )
      {
         SCIP_CALL( SCIPlexicographicReductionAddPermutation(scip, (*targetdata)->lexreddata,
               (*targetdata)->x_vars, (*targetdata)->res_graph->nnodes, (*targetdata)->perms[p],
               SYM_SYMTYPE_PERM, NULL, TRUE, &success) );
         assert(success);
      }
   }

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
   int                   k,                   /**< parameter k for k-vertex cut */
   int*                  vertex_weights      /**< array of vertex weights */
   )
{
   SCIP_PROBDATA* probdata;
   probdata = NULL;

   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(nnodes > 0);
   assert(nedges >= 0);
   assert(tail != NULL || nedges == 0);
   assert(head != NULL || nedges == 0);
   assert(vertex_weights != NULL || nnodes == 0);

   /* allocate memory for problem data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &probdata) );

   /* store basic graph information */
   probdata->k = k;

   /* initialize symmetry informatio */
   probdata->symconss = NULL;
   probdata->nsymconss = 0;
   probdata->orbitalreddata = NULL;
   probdata->lexreddata = NULL;
   probdata->perms = NULL;
   probdata->nperms = -1;
   probdata->lenperms = -1;

   /* Create graph - IMPORTANT: pass true to make Graph copy the arrays! */
   cout << "Creating Graph structure with " << nnodes << " nodes and " << nedges << " edges..." << endl;
   probdata->orig_graph = new Graph(nnodes, nedges, tail, head, vertex_weights, true);
   
   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   // SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, "DENTRO_NODO", "Does stuff inside a node", eventExecMyLogic, NULL) );

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
   unsigned int LPonly;
   SCIPgetBoolParam(scip, "options/solveLP", &LPonly);
   
   if( !LPonly )
   {
      SCIP_CALL( SCIPsetObjIntegral(scip) );
   }
   

   unsigned int relax_constraints;
   SCIPgetBoolParam(scip, "options/relaxconstraints", &relax_constraints);

   /* Allocate memory for preprocessing data */
   probdata->preFixed = new bool[nnodes];
   
   /* Perform preprocessing */
   cout << "\n=== Starting Preprocessing ===" << endl;
   probdata->n_fixed = performPreprocessing(*probdata->orig_graph, k, probdata->preFixed);
   
   /* Count connected components in the residual graph after preprocessing */
   cout << "\nCounting connected components after preprocessing..." << endl;
   int num_SCC = countConnectedComponents(*probdata->orig_graph, probdata->preFixed);
   cout << "Found " << num_SCC << " connected components\n";

   cout << "\n\nNUMBER OF VERTICES FIXED IN PREPROCESSING: " << probdata->n_fixed << "\n\n";

   if( k <= num_SCC)
   {
      cout<<"-----------------------------------------------------------------------------------------\n";
      cout << "The instance is trivial since k <= number of connected components after pre-fixing (" << num_SCC << ")." << endl;
      cout << "Exiting..." << endl;
      cout << "-----------------------------------------------------------------------------------------\n";

      /* Create residual graph */
      vector<int> old_to_new;
      vector<int> new_to_old;
      cout << "\nCreating residual graph..." << endl;
      probdata->res_graph = buildResidualGraph(*probdata->orig_graph, probdata->preFixed, old_to_new, new_to_old);
      
      /**********************************
       * Write results in a file *
      **********************************/
      int k_val;
      int sym_method;
      unsigned int plot_option;
      char* resultsFileName;
      char* solutionsFileName;
      char* plotFileName;
      unsigned int option_conn_warmstart;
      unsigned int option_conn_cut;
      unsigned int weighted;
      int clique_gen_method;
      unsigned int option_ILP_warmstart;

      SCIPgetStringParam(scip, "output/results", &resultsFileName);
      SCIPgetStringParam(scip, "output/solution", &solutionsFileName);
      SCIPgetStringParam(scip, "output/plot", &plotFileName);
      SCIPgetIntParam(scip, "params/k", &k_val);
      SCIPgetIntParam(scip, "options/symmetrymethod", &sym_method);
      SCIPgetBoolParam(scip, "options/plot", &plot_option);
      SCIPgetBoolParam(scip, "options/connectivity/warmstart", &option_conn_warmstart);
      SCIPgetBoolParam(scip, "options/connectivity/cut", &option_conn_cut);
      SCIPgetBoolParam(scip, "options/weighted", &weighted);
      SCIPgetIntParam(scip, "options/cliquegeneration", &clique_gen_method);
      SCIPgetBoolParam(scip, "options/ILPwarmstart", &option_ILP_warmstart);

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
      
      int tot_cut_cost = 0;
      bool* cut = new bool[probdata->orig_graph->nnodes];

      // Initialize cut to 0
      for (int i = 0; i < probdata->orig_graph->nnodes; ++i) {
         cut[i] = false;
         if(probdata->preFixed[i]) {
            tot_cut_cost += probdata->orig_graph->vertex_weights[i];
            cut[i] = true;
         }
      }
      
      ofstream resultsFile(resultsFileName, ios::app);
      if (resultsFile.is_open())
      {

         resultsFile 
                  << SCIPgetProbName(scip) << "\t"
                  << probdata->orig_graph->nnodes << "\t"
                  << probdata->orig_graph->nedges << "\t"
                  << probdata->res_graph->nnodes << "\t"
                  << probdata->res_graph->nedges << "\t"
                  << k_val << "\t"
                  << sym_method << "\t"
                  << option_conn_warmstart << "\t"
                  << option_ILP_warmstart << "\t"
                  << option_conn_cut << "\t"
                  << clique_gen_method << "\t"
                  << relax_constraints << "\t"
                  << LPonly << "\t"
                  << weighted << "\t"
                  << probdata->n_fixed << "\t"
                  << tot_cut_cost << "\t"
                  << "Trivial" << "\t"
                  << tot_cut_cost << "\t"
                  << tot_cut_cost << "\t"
                  << tot_cut_cost << "\t"
                  << 0.01 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << -1 << "\t"
                  << endl;
      resultsFile.close();
      }
      else
      {
         cout << "\nError: could not open results file " << resultsFileName << endl;
      }

      // Write the solution in the format tot_cut_cost \n vertices in cut separated by space
      ofstream solFile(solutionsFileName);
      if (solFile.is_open()) {
         solFile << tot_cut_cost << "\n";
         for (int v = 0; v < probdata->orig_graph->nnodes; ++v) {
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

      // Plot the solution graph
      SCIP_CALL( SCIPprobdataPlotSolution(probdata->orig_graph, cut, plotFileName, false, plot_option) );
   

      delete[] resultsFileName;
      delete[] solutionsFileName;
      delete[] plotFileName;
      delete[] cut;
         
      exit(0);
   }

   /* Identify isolated cliques (includes singletons) */
   cout << "Identifying isolated cliques..." << endl;
   bool* is_in_isolated_clique = new bool[nnodes];
   int num_isolated_cliques = identifyIsolatedCliques(*probdata->orig_graph, probdata->preFixed, is_in_isolated_clique);
   cout << "Number of isolated cliques found: " << num_isolated_cliques << endl;
   cout << "=========================================================================\n";
      
   // Update k by subtracting the number of nodes in isolated cliques
   probdata->k -= num_isolated_cliques;
   cout << "Updated k = " << probdata->k << "\n" << endl;

   // Build to_remove array combining preFixed and isolated cliques
   bool* to_remove = new bool[nnodes];
   for(int v = 0; v < nnodes; ++v)
      to_remove[v] = probdata->preFixed[v] || is_in_isolated_clique[v];
   
   // Clean up
   delete[] is_in_isolated_clique;

   /* Create residual graph */
   vector<int> old_to_new;
   vector<int> new_to_old;
   cout << "\nCreating residual graph..." << endl;
   probdata->res_graph = buildResidualGraph(*probdata->orig_graph, to_remove, old_to_new, new_to_old);
   delete [] to_remove;


   unsigned int option_conn_cut;
   SCIPgetBoolParam(scip, "options/connectivity/cut", &option_conn_cut);

   int min_connectivity = 1;

   /* Compute vertex connectivity of the residual graph (LB on the optimal value) */
   if(option_conn_cut) {
      bool* separator = new bool[probdata->res_graph->nnodes];
      min_connectivity = computeVertexConnectivity(*probdata->res_graph, true, separator);
      delete [] separator;
   }

   int clique_gen_method;
   SCIPgetIntParam(scip, "options/cliquegeneration", &clique_gen_method);

   /* Build edge covering cliques */
   cout << "\nBuilding edge covering cliques..." << endl;
   int n_cliques = probdata->res_graph->computeEdgeCoveringCliques(clique_gen_method);
   
   cout<<"\n-------------------------------------------------------------------------\n";
   cout << "NUMBER OF CLIQUES ADDED FOR EDGE COVERING: " << n_cliques << endl;
   cout << "MIN_CONNECTIVITY of residual graph: " << min_connectivity << endl;
   cout << "-------------------------------------------------------------------------\n\n";
   
   /* Store min_connectivity in probdata */
   probdata->min_connectivity = min_connectivity;

   /* FORMULATION */
   SCIP_CONS** vertex_cover_constrs;
   SCIP_CONS** clique_constrs;
   
   /* First create variables for each node (x_v) */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->x_vars, probdata->res_graph->nnodes) );

   for(int v = 0; v < probdata->res_graph->nnodes; ++v )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", v);

      if( !LPonly )
      {
         /* Create binary variable x_v */
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->x_vars[v], name, 0.0, 1.0, (SCIP_Real)probdata->res_graph->vertex_weights[v], SCIP_VARTYPE_BINARY) );
      }
      else
      {
         /* Create continuous variable x_v in [0,1] */
         SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->x_vars[v], name, 0.0, 1.0, (SCIP_Real)probdata->res_graph->vertex_weights[v], SCIP_VARTYPE_CONTINUOUS) );
      }

      SCIP_CALL( SCIPaddVar(scip, probdata->x_vars[v]) );

      // Uncomment for giving branching priority based on node degree (very ineffective)
      // SCIP_CALL(SCIPchgVarBranchPriority(scip, probdata->x_vars[i], probdata->res_graph->node_degree[i]));
   }

   /* Constraint (14): Cardinality constraint for α variables */
   /* expr: ∑_S α_S >= k */
   SCIP_CONS* alpha_cardinality_constr;
   if(relax_constraints)
   {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &alpha_cardinality_constr, "alpha_cardinality_constraint", 0, NULL, NULL, (SCIP_Real)probdata->k, SCIPinfinity(scip)) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &alpha_cardinality_constr, "alpha_cardinality_constraint", 0, NULL, NULL, (SCIP_Real)probdata->k, (SCIP_Real)probdata->k) );
   }

   SCIP_CALL(SCIPsetConsModifiable(scip, alpha_cardinality_constr, TRUE));
   SCIP_CALL(SCIPaddCons(scip, alpha_cardinality_constr));
   
   /* Constraint: Minimum connectivity lower bound */
   /* expr: ∑_v x_v >= min_connectivity */
   /* Notice that this constraint will be useful when connectivity enhancements are enabled, otherwise it will be trivially satisfied by any feasible solution */
   SCIP_CONS* min_connectivity_constr;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &min_connectivity_constr, "min_connectivity_constraint", 0, NULL, NULL, (SCIP_Real)min_connectivity, SCIPinfinity(scip)) );
   
   /* Add all x variables to the constraint */
   for(int v = 0; v < probdata->res_graph->nnodes; ++v )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, min_connectivity_constr, probdata->x_vars[v], probdata->res_graph->vertex_weights[v]) );
   }
   
   SCIP_CALL(SCIPaddCons(scip, min_connectivity_constr));

   /* Alloc memory for storing constraints */
   probdata->n_cliques = n_cliques;
   SCIP_CALL( SCIPallocBufferArray(scip, &vertex_cover_constrs, probdata->res_graph->nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &clique_constrs, n_cliques) );


   for(int v = 0; v < probdata->res_graph->nnodes; ++v )
   {
      /* Constraint: ∑_S(v) α_S + x_v >= 1  */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "coverage_constr_%d", v);

      if(relax_constraints)
      {
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &vertex_cover_constrs[v], name, 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
      }
      else
      {
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &vertex_cover_constrs[v], name, 0, NULL, NULL, 1.0, 1.0) );
      }

      SCIP_CALL( SCIPaddCoefLinear(scip, vertex_cover_constrs[v], probdata->x_vars[v], 1.0) );

      /* This constraint will be modifiable for adding α_S variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, vertex_cover_constrs[v], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, vertex_cover_constrs[v]) );
   }


   for(int c = 0; c < n_cliques; ++c )
   {
      /* Constraint: ∑_S(C) α_S <= 1  */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "clique_constr_%d", c);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &clique_constrs[c], name, 0, NULL, NULL, -SCIPinfinity(scip), 1.0) );

      /* This constraint will be modifiable for adding α_S variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, clique_constrs[c], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, clique_constrs[c]) );
   }


   probdata->alpha_cardinality_constr = alpha_cardinality_constr;
   probdata->min_connectivity_constr = min_connectivity_constr;

   /* store constraints */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &probdata->vertex_cover_constrs, vertex_cover_constrs, probdata->res_graph->nnodes) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &probdata->clique_constrs, clique_constrs, n_cliques) );

   /* initialize variables arrays */
   probdata->alphavars = NULL;
   probdata->n_alphavars = 0;
   probdata->alphavars_size = 0; 

   probdata->clique_cuts = NULL;
   probdata->n_clique_cuts = 0;
   probdata->clique_cuts_size = 0;

   
   /* initialize singleton columns*/
   SCIP_CALL(createInitialColumns(scip, probdata));

   unsigned int option_conn_warmstart;
   SCIPgetBoolParam(scip, "options/connectivity/warmstart", &option_conn_warmstart);

   unsigned int option_ILP_warmstart;
   SCIPgetBoolParam(scip, "options/ILPwarmstart", &option_ILP_warmstart);

   /* If connectivity warm start is enabled, run the connectivity heuristic */
   if(option_conn_warmstart && !option_ILP_warmstart)
   {
      SCIP_CALL(warmStartHeuristic(scip, probdata, false));
   }
   /* If the ILP warm start is enabled, run the ILP heuristic */
   else if(option_ILP_warmstart)
   {  
      SCIP_CALL(warmStartHeuristic(scip, probdata, true));
   }

   /* Notice that the ILP heuristic "dominates" the connectivity heuristic. If both are enabled, only the ILP heuristic is run. */

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* activate pricer for k-vertex cut */
   SCIP_CALL( SCIPpricerKvertexcutActivate(scip, 
                                           probdata->alpha_cardinality_constr, probdata->vertex_cover_constrs, probdata->clique_constrs, 
                                           probdata->res_graph->nnodes, probdata->res_graph->nedges, probdata->n_cliques) );


   /* add symmetry handling constraints */
   int symmetrymethod;
   SCIP_CALL( SCIPgetIntParam(scip, "options/symmetrymethod", &symmetrymethod) );
   if( symmetrymethod == 1 )
   {
      SCIP_CALL( tryAddSymmetryHandlingConss(scip, probdata->x_vars, &probdata->symconss, &probdata->nsymconss,
            probdata->res_graph->nnodes, probdata->res_graph->nedges,
            probdata->res_graph->head, probdata->res_graph->tail) );
   }
   else if( symmetrymethod == 2 )
   {
      SCIP_EVENTHDLR* shadowtree = NULL;

      SCIP_CALL( SCIPincludeEventHdlrShadowTree(scip, &shadowtree) );
      assert(shadowtree != NULL );

      SCIP_CALL( SCIPincludeOrbitalReduction(scip, &probdata->orbitalreddata, shadowtree) );
      SCIP_CALL( SCIPincludeLexicographicReduction(scip, &probdata->lexreddata, shadowtree) );

      SCIP_CALL( computeSymmetries(scip, probdata->x_vars, &probdata->perms, &probdata->nperms,
            &probdata->lenperms, probdata->res_graph->nnodes, probdata->res_graph->nedges,
            probdata->res_graph->head, probdata->res_graph->tail) );
      printSymmetryInfo(probdata->nperms, 2);
   }

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &vertex_cover_constrs);
   SCIPfreeBufferArray(scip, &clique_constrs);

   SCIPinfoMessage(scip, NULL, "Created k-vertex cut problem: %d nodes, %d edges, k=%d\n", 
                                probdata->res_graph->nnodes, probdata->res_graph->nedges, probdata->k);

   /* write problem formulation to .lp file for debugging */
   // SCIP_CALL( SCIPwriteOrigProblem(scip, "formulation.lp", "lp", FALSE) );

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

   /* Get nnodes from graph before deleting it */
   int nnodes = (*probdata)->res_graph->nnodes;

   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->x_vars[i]) );
   }
   
   for( i = 0; i < (*probdata)->n_alphavars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->alphavars[i]) );
   }

   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->vertex_cover_constrs[i]) );
   }

   for (i = 0; i < (*probdata)->n_cliques; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->clique_constrs[i]) );
   }

   /* release symmetry constraints */
   for( i = 0; i < (*probdata)->nsymconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->symconss[i]) );
   }

   /* release cardinality constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->alpha_cardinality_constr) );
   
   /* release minimum connectivity constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->min_connectivity_constr) );
   
   /* free memory of arrays */
   if( (*probdata)->n_alphavars > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->alphavars, (*probdata)->alphavars_size);
   }

   /* possibly free symmetry data */
   if( (*probdata)->perms != NULL )
   {
      for( i = 0; i < (*probdata)->nperms; ++i )
      {
         SCIPfreeBlockMemoryArray(scip, &(*probdata)->perms[i], (*probdata)->res_graph->nnodes);
      }
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->perms, (*probdata)->lenperms);

      SCIP_CALL( SCIPlexicographicReductionFree(scip, &(*probdata)->lexreddata) );
      SCIP_CALL( SCIPorbitalReductionFree(scip, &(*probdata)->orbitalreddata) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->x_vars, nnodes);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->vertex_cover_constrs, nnodes);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->clique_constrs, (*probdata)->n_cliques);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->clique_cuts, (*probdata)->clique_cuts_size);
   
   /* free preprocessing data */
   if( (*probdata)->preFixed != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->preFixed, (*probdata)->orig_graph->nnodes);
   }
   
   /* delete graph structure (calls destructor) */
   delete (*probdata)->orig_graph;
   delete (*probdata)->res_graph;

   if( (*probdata)->nsymconss > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->symconss, (*probdata)->nsymconss);
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


SCIP_VAR** SCIPprobdataGetAlphaVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->alphavars;
}

/** returns array of clique cuts */
SCIP_ROW** SCIPprobdataGetCliqueCuts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->clique_cuts;
}

/**returns number of clique cuts */
int SCIPprobdataGetNCliqueCuts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->n_clique_cuts;
}


SCIP_VAR** SCIPprobdataGetXVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->x_vars;
}


/** returns number of nodes */
int SCIPprobdataGetNNodes(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->nnodes;
}

/** returns number of edges */
int SCIPprobdataGetNEdges(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->nedges;
}

/** returns array of edge tail nodes */
int* SCIPprobdataGetEdgeTails(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->tail;
}

/** returns array of edge head nodes */
int* SCIPprobdataGetEdgeHeads(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->head;
}

/** returns adjacency list for a given node */
int* SCIPprobdataGetAdjList(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   assert(node >= 0 && node < probdata->res_graph->nnodes);
   return probdata->res_graph->adj_lists[node];
}


/** returns all adjacency lists */
int** SCIPprobdataGetAdjLists(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->adj_lists;
}

/** returns size of adjacency list for a given node */
int SCIPprobdataGetAdjSize(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   assert(node >= 0 && node < probdata->res_graph->nnodes);
   return probdata->res_graph->node_degree[node];
}

/** returns size of adjacency list for all nodes */
int* SCIPprobdataGetAdjSizes(
   SCIP_PROBDATA*        probdata           /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph->node_degree;
}

/** returns degree of a given node */
int SCIPprobdataGetNodeDegree(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   assert(node >= 0 && node < probdata->res_graph->nnodes);
   return probdata->res_graph->node_degree[node];
}

/** returns parameter k */
int SCIPprobdataGetk(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->k;
}

/** returns minimum connectivity value */
int SCIPprobdataGetMinConnectivity(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->min_connectivity;
}

/** returns graph structure */
Graph* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   return probdata->res_graph;
}


/** returns graph structure */
Graph* SCIPprobdataGetOrigGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->orig_graph != NULL);
   return probdata->orig_graph;
}

/** checks if two nodes are adjacent */
SCIP_Bool SCIPprobdataAreAdjacent(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node1,              /**< first node */
   int                   node2               /**< second node */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);
   assert(node1 >= 0 && node1 < probdata->res_graph->nnodes);
   assert(node2 >= 0 && node2 < probdata->res_graph->nnodes);
   return (probdata->res_graph->adj_matrix[node1][node2] ? TRUE : FALSE);
}

/** returns adjacency matrix */
bool** SCIPprobdataGetAdjMatrix(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   assert(probdata->res_graph != NULL);

   return probdata->res_graph->adj_matrix;
}


/** returns alpha cardinality constraint */
SCIP_CONS* SCIPprobdataGetAlphaCardinalityConstr(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->alpha_cardinality_constr;
}


/** returns |V|-set of coverage constraints */
SCIP_CONS** SCIPprobdataGetVertexCoverConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->vertex_cover_constrs;
}

/** returns |C|-set of clique constraints */
SCIP_CONS** SCIPprobdataGetCliqueConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->clique_constrs; 
}

/** returns preprocessing fixed nodes array */
bool* SCIPprobdataGetPreFixed(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->preFixed;
}

/** returns number of fixed nodes in preprocessing */
int SCIPprobdataGetNFixed(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->n_fixed;
}

/** returns number of cliques */
int SCIPprobdataGetNCliques(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->n_cliques;
}


/** returns lexicographic reduction data */
SCIP_LEXREDDATA* SCIPgetProbdataLexreddata(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->lexreddata;
}

/** returns orbital reduction data */
SCIP_ORBITALREDDATA* SCIPgetProbdataOrbitalreddata(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);
   return probdata->orbitalreddata;
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


/** adds given cut to the problem data */
SCIP_RETCODE SCIPprobdataAddCliqueCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_ROW*             cut,              /**< cut to add */
   int*                  clique             /**< clique associated with the cut */
   )
{

   /* check if enough memory is left */
   if( probdata->clique_cuts_size == probdata->n_clique_cuts )
   {
      int newsize;
      newsize = MAX(100, probdata->clique_cuts_size * 2);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->clique_cuts, probdata->clique_cuts_size, newsize) );
      probdata->clique_cuts_size = newsize;
   }

   probdata->clique_cuts[probdata->n_clique_cuts] = cut;
   probdata->clique_cuts_cliques.push_back(clique);

   probdata->n_clique_cuts++;

   
   return SCIP_OKAY;
}


/* retrieves solution from SCIP and maps it back to original graph */
bool* SCIPprobdataGetSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int&                  tot_cut_cost                 /**< reference to store total cut cost */
   )
{

   bool* is_in_isolated_clique = new bool[probdata->orig_graph->nnodes];
   identifyIsolatedCliques(*probdata->orig_graph, probdata->preFixed, is_in_isolated_clique);
      
   // Build to_remove array combining preFixed and isolated cliques
   bool* to_remove = new bool[probdata->orig_graph->nnodes];
   for(int v = 0; v < probdata->orig_graph->nnodes; ++v)
      to_remove[v] = probdata->preFixed[v] || is_in_isolated_clique[v];
   
   // Clean up
   delete[] is_in_isolated_clique;

   bool* cut = new bool[probdata->orig_graph->nnodes];
   for(int v = 0; v < probdata->orig_graph->nnodes; ++v)
      cut[v] = false;


   /* Create residual graph */
   vector<int> old_to_new;
   vector<int> new_to_old;

   Graph *res_graph = buildResidualGraph(*probdata->orig_graph, to_remove, old_to_new, new_to_old);
   
   tot_cut_cost = 0;

   delete [] to_remove;
   
   if (SCIPgetBestSol(scip) != NULL)
   {
      for (int v = 0; v < res_graph->nnodes; ++v)
      {
         cut[new_to_old[v]] = (SCIPgetSolVal(scip, SCIPgetBestSol(scip), probdata->x_vars[v]) > 0.5);
         if (cut[new_to_old[v]])
            tot_cut_cost += res_graph->vertex_weights[v];
      }

      for(int v = 0; v < probdata->orig_graph->nnodes; ++v)
      {
         if (probdata->preFixed[v])
         {
            cut[v] = true;
            tot_cut_cost += probdata->orig_graph->vertex_weights[v];
         }
      }
      
      return cut;
   }
   else
   {
      SCIPerrorMessage("No best solution found.\n");
      return cut;
   }
}


/* Given a cut, plots the residual graph and saves the png in the specified location */
SCIP_RETCODE SCIPprobdataPlotSolution(
   Graph*                graph,           /**< instance graph */
   bool*                 cut,                /**< cut array */
   const char*           filename,            /**< output filename */
   bool                  draw_attacked_edges,       /**< whether to highlight cut nodes */
   bool                  plot_graph            /**< whether to plot the graph */
   )
{
   int nnodes = graph->nnodes;
   int nedges = graph->nedges;

   std::vector<std::vector<int>> adj(nnodes);
   for (int e = 0; e < nedges; ++e) {
      adj[graph->tail[e]].push_back(graph->head[e]);
      adj[graph->head[e]].push_back(graph->tail[e]);
   }

   // --- 1: Identify components (Green vs Red) ---
   
   std::vector<bool> visited(nnodes, false);
   std::vector<std::vector<int>> green_components;
   std::vector<int> red_nodes;

   // Isolate red nodes
   for(int i = 0; i < nnodes; ++i) {
      if(cut[i]) {
         visited[i] = true; 
         red_nodes.push_back(i);
      }
   }

   // Find green connected components
   for (int i = 0; i < nnodes; ++i)
   {
      if (!visited[i])
      {
         std::vector<int> component;
         std::queue<int> q;
         visited[i] = true;
         q.push(i);
         component.push_back(i);

         while(!q.empty())
         {
               int u = q.front(); q.pop();
               for(int v : adj[u]) {
                  if(!visited[v]) {
                     visited[v] = true;
                     q.push(v);
                     component.push_back(v);
                  }
               }
         }
         green_components.push_back(component);
      }
   }

   // Sort green components from largest to smallest
   std::sort(green_components.begin(), green_components.end(), 
            [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });

   // --- 2: Calculate Coordinates ---
   
   struct Point { double x, y; };
   std::vector<Point> node_pos(nnodes);
   
   // Graphic constants
   const double NODE_SPACING = 70.0; // Space between nodes in the same component
   const double COMP_PADDING = 100.0; // Extra space around each component
   const double ORBIT_PADDING = 150.0; // Space between orbits

   double max_y_reached = 0.0; // To know where to place red nodes at the end

   // 2.1 Placement of Green Components
   if (!green_components.empty())
   {
      // --- THE SUN (Largest component) ---
      // Place it at (0,0)
      double sun_radius = (green_components[0].size() * NODE_SPACING) / (2 * M_PI);
      if (green_components[0].size() == 1) sun_radius = 0.0;
      else if (sun_radius < 50.0) sun_radius = 50.0;

      for (size_t k = 0; k < green_components[0].size(); ++k) {
         int node_idx = green_components[0][k];
         if (green_components[0].size() == 1) {
               node_pos[node_idx] = {0.0, 0.0};
         } else {
               double angle = 2.0 * M_PI * k / green_components[0].size();
               node_pos[node_idx].x = sun_radius * cos(angle);
               node_pos[node_idx].y = sun_radius * sin(angle);
         }
      }
      max_y_reached = sun_radius;

      // --- THE PLANETS (Subsequent components) ---
      // Orbital algorithm
      double current_orbit_radius = sun_radius + ORBIT_PADDING;
      double current_angle = 0.0;
      double max_radius_in_this_orbit = 0.0;

      for (size_t i = 1; i < green_components.size(); ++i)
      {
         int n_comp = green_components[i].size();
         
         // Local radius of the component
         double local_radius = (n_comp > 1) ? (n_comp * NODE_SPACING) / (2 * M_PI) : 0.0;
         if (local_radius < 30.0 && n_comp > 1) local_radius = 30.0;
         
         // Update the maximum radius found in this orbit (to calculate the next one)
         if (local_radius > max_radius_in_this_orbit) max_radius_in_this_orbit = local_radius;

         // Calculate how much angular space this component occupies on the current orbit
         // Angle = arc length / orbit radius
         // Add a safety margin
         double diam_with_padding = (local_radius * 2) + COMP_PADDING;
         double angular_width = diam_with_padding / current_orbit_radius;

         // If the component doesn't fit in the remaining angle (exceeds 2PI), create a new orbit
         if (current_angle + angular_width > 2 * M_PI)
         {
               // Move to the next orbit
               current_orbit_radius += (max_radius_in_this_orbit * 2) + ORBIT_PADDING;
               current_angle = 0.0;
               max_radius_in_this_orbit = local_radius; // Reset for the new orbit
               // Recalculate angular width on the new radius (which is larger, so angle is smaller)
               angular_width = diam_with_padding / current_orbit_radius;
         }

         // Center of the component (Planet)
         double planet_center_x = current_orbit_radius * cos(current_angle);
         double planet_center_y = current_orbit_radius * sin(current_angle);

         // Place the nodes of the planetary component
         for (int k = 0; k < n_comp; ++k)
         {
               int node_idx = green_components[i][k];
               if (n_comp == 1) {
                  node_pos[node_idx] = {planet_center_x, planet_center_y};
               } else {
                  double local_angle = 2.0 * M_PI * k / n_comp;
                  node_pos[node_idx].x = planet_center_x + local_radius * cos(local_angle);
                  node_pos[node_idx].y = planet_center_y + local_radius * sin(local_angle);
               }
         }

         // Update the angle for the next component
         current_angle += angular_width;
         
         if (current_orbit_radius + local_radius > max_y_reached) 
               max_y_reached = current_orbit_radius + local_radius;
      }
   }

   // 2.2 Placement of Red Nodes (Asteroid Belt)
   // We place them at the top, separated by a wide margin
   if (!red_nodes.empty())
   {
      double red_start_y = max_y_reached + 100.0; // Lots of vertical space
      double red_grid_spacing = 60.0;
      int n_red = red_nodes.size();
      
      // Calculate grid width to center it
      int cols = (int)ceil(sqrt(n_red)) + 5; // A bit wide rectangular
      double total_width = cols * red_grid_spacing;
      double start_x = -(total_width / 2.0); // Centered with respect to x=0

      for (int k = 0; k < n_red; ++k)
      {
         int r = k / cols;
         int c = k % cols;
         
         int node_idx = red_nodes[k];
         node_pos[node_idx].x = start_x + c * red_grid_spacing;
         node_pos[node_idx].y = red_start_y + r * red_grid_spacing;
      }
   }

   // --- 3: Writing DOT File ---

   string dot_filename = string(filename) + ".dot";
   ofstream dotfile(dot_filename);
   
   if (!dotfile.is_open()) return SCIP_ERROR;
   
   dotfile << "graph G {" << endl;
   dotfile << "  splines=false;" << endl; 
   dotfile << "  pad=0.5;" << endl; 
   dotfile << "  outputorder=edgesfirst;" << endl;
   dotfile << "  dpi=300;" << endl;
   
   // Standard nodes and edges style
   dotfile << "  node [shape=circle, style=filled, fontsize=10, fixedsize=true, width=0.4];" << endl;
   dotfile << "  edge [penwidth=1.0];" << endl;
   dotfile << endl;
   dotfile << std::fixed << std::setprecision(2);

   for (int v = 0; v < nnodes; ++v)
   {
      string color = cut[v] ? "lightcoral" : "palegreen";
      
      dotfile << "  " << v << " [label=\"" << (v + 1) << "\", "
               << "fillcolor=" << color 
               << ", pos=\"" << node_pos[v].x << "," << node_pos[v].y << "!\""
               << "];" << endl;
   }
   
   dotfile << endl;
   
   // Scrittura Archi
   for (int e = 0; e < nedges; ++e)
   {
      int u = graph->tail[e];
      int v = graph->head[e];
      
      bool is_attacked = cut[u] || cut[v];
      
      if (is_attacked) {
            if (draw_attacked_edges) {
               // Very light gray to not disturb the view
               dotfile << "  " << u << " -- " << v << " [style=dashed, color=\"#BBBBBB\"];" << endl;
            }
      } else {
            dotfile << "  " << u << " -- " << v << " [color=black];" << endl;
      }
   }
   
   dotfile << "}" << endl;
   dotfile.close();
   
   // --- 4: Rendering ---
   if(!plot_graph) return SCIP_OKAY;

   string command = "neato -n2 -Tpng " + dot_filename + " -o " + string(filename);
   int result = system(command.c_str());
   
   if (result != 0) {
      remove(dot_filename.c_str());
      return SCIP_ERROR;
   }
   remove(dot_filename.c_str());
   return SCIP_OKAY;
}

/**@} */
