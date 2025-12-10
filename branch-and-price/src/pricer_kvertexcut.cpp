/**@file   pricer_kvertexcut.cpp
 * @brief  k-vertex cut variable pricer
 * @author Fabio Ciccarelli
 *
 * This file implements the variable pricer which check if variables exist with negative reduced cost. 
 * If such variables exist, they are added to the LP.
 *
**/

#include <assert.h>
#include <string.h>

#include "scip/cons_varbound.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"

#include "pricer_kvertexcut.h"
#include "probdata_kvertexcut.h"
#include "vardata_kvertexcut.h"
#include "min_cut_pricer.h"

#include <iostream>
#include <algorithm>
using namespace std;

/** Pricer properties
 *
 * 
 **/
#define PRICER_NAME            "kvertexcut"
#define PRICER_DESC            "pricer for k-vertex cut tours"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           FALSE     /* only call pricer if all problem variables have non-negative reduced costs */

/**@} */


/*
 * Data structures
 */

/** @brief Variable pricer data used in the \ref pricer_kvertexcut.cpp "pricer" */
struct SCIP_PricerData
{
   SCIP_CONS*            alpha_cardinality_constr;   /**< cardinality constraint for α variables */
   SCIP_CONS**           vertex_cover_constrs;         /**< array of vertex cover constraints */
   SCIP_CONS**           clique_constrs;           /**< array of clique constraints */
   int                   nnodes;                     /**< number of nodes in the graph */
   int                   nedges;                     /**< number of edges in the graph */
   int                   n_cliques;                  /**< number of cliques */
   MinCutPricer*         min_cut_pricer;            /**< min cut based pricer */
   SCIP_NODE*            lastnode;                   /**< last processed node */
   int                   nodecnt;                    /**< node counter */
   int                   n_gen_cols_root_node;       /**< number of generated columns at root node */
   int                   n_gen_farkas_cols;          /**< number of generated Farkas columns */
};



/*
 * Data structures
 *//** perform pricing using min-cut algorithm */
static
SCIP_RETCODE doPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< variable pricer structure */
   SCIP_Bool             isfarkas,           /**< whether to use Farkas pricing */
   SCIP_Real* lowerbound, 
   SCIP_Bool* stopearly,
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS*            alpha_cardinality_constr;   /**< cardinality constraint for α variables */
   SCIP_CONS**           vertex_cover_constrs;         /**< array of vertex cover constraints for SCCs */
   SCIP_CONS**           clique_constrs;           /**< array of clique constraints */
   SCIP_PROBDATA* probdata;

  
   int nnodes, nedges;

   assert(scip != NULL);
   assert(pricer != NULL);
   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* get the pricer data */
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   if( SCIPgetCurrentNode(scip) != pricerdata->lastnode )
   {
      pricerdata->lastnode = SCIPgetCurrentNode(scip);
      pricerdata->nodecnt = 0;
   }
   // else if( pricerdata->nodecnt > 100 )
   //    return SCIP_OKAY;

   /* get the problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);

   int n_cliques = SCIPprobdataGetNCliques(probdata);
   const std::vector<std::vector<bool>>* cliques = SCIPprobdataGetGraph(probdata)->getCliques();

   nnodes = pricerdata->nnodes;
   nedges = pricerdata->nedges;
   alpha_cardinality_constr = pricerdata->alpha_cardinality_constr;
   vertex_cover_constrs = pricerdata->vertex_cover_constrs;
   clique_constrs = pricerdata->clique_constrs;

   SCIP_NODE* node = SCIPgetCurrentNode(scip);
   int num_node = SCIPnodeGetNumber(node);

   // SCIP_NODE* parentnode;

   // if(num_node > 1)
   //    parentnode = SCIPnodeGetParent(pricerdata->lastnode);
   
   
   SCIP_Real mu_val = isfarkas ? SCIPgetDualfarkasLinear(scip, alpha_cardinality_constr) : SCIPgetDualsolLinear(scip, alpha_cardinality_constr);
   
   // cout << "Pricing at node " << num_node << " with π = " << pi_val << endl;
   // cin.get();


   double* node_dual_coeffs = new double[nnodes];
   for( int i = 0; i < nnodes; ++i )
   {
      node_dual_coeffs[i] = isfarkas ? SCIPgetDualfarkasLinear(scip, vertex_cover_constrs[i]) : SCIPgetDualsolLinear(scip, vertex_cover_constrs[i]);
      // cout << "Node " << i << " dual coeff = " << node_dual_coeffs[i] << endl;
   }

   double* cl_dual_coeffs = new double[n_cliques];
   for( int c = 0; c < n_cliques; ++c )
   {
      cl_dual_coeffs[c] = isfarkas ? SCIPgetDualfarkasLinear(scip, clique_constrs[c]) : SCIPgetDualsolLinear(scip, clique_constrs[c]);
      // cout << "dual coeff = " << cl_dual_coeffs[c] << endl;
   }

   vector<int> down_branched_vertices;
   vector<int> up_branched_vertices;

   /* Update coefficients of x variables (node variables) with π value */
   for( int v = 0; v < nnodes; ++v )
   {  
      SCIP_Real current_lb = SCIPvarGetLbLocal(x_vars[v]);
      SCIP_Real current_ub = SCIPvarGetUbLocal(x_vars[v]);
                 
      
      if( current_lb > 0.5)
      {
         up_branched_vertices.push_back(v);
      }
      else if( current_ub < 0.5 )
      {
         down_branched_vertices.push_back(v);
      }
   }

   int **adj_lists = SCIPprobdataGetAdjLists(probdata);
   int *adj_sizes = SCIPprobdataGetAdjSizes(probdata);


   /* --- EARLY BRANCHING: UNCOMMENT FOR FUTURE TESTING --- */

   // if(!isfarkas)
   // {
   //    if( pricerdata->nodecnt >= 1000 )
   //    {
   //       // cout << "EARLY BRANCHING at node " << num_node << " after " << pricerdata->nodecnt << " generated columns" << endl;

   //       PoolResult pricing_result = pricerdata->min_cut_pricer->solve(nnodes, nedges, pi_val, mu_val,
   //                                                               edge_dual_coeffs, node_dual_coeffs, cl_dual_coeffs,
   //                                                               adj_lists, adj_sizes,
   //                                                               down_branched_vertices, up_branched_vertices, true);

   //       double lagrangian_bound = 0.0;
         
   //       if( pricing_result.num_solutions > 0 )
   //       {
   //          lagrangian_bound = SCIPgetLPObjval(scip) - pricing_result.obj_vals[0]*SCIPprobdataGetk(probdata);
   //       }

   //       // cout << "Lagrangian bound = " << lagrangian_bound << endl;

   //       (*result) = SCIP_SUCCESS;
   //       if(num_node > 1)
   //          (*lowerbound) = (SCIP_Real) max(max((double)SCIPnodeGetLowerbound(parentnode), 1.0), lagrangian_bound);
   //       else
   //          (*lowerbound) = (SCIP_Real) max(1.0, lagrangian_bound);

   //       (*stopearly) = TRUE;
   //       return SCIP_OKAY;
   //    }
   // }


   /* Solve the pricing problem */
   PoolResult pricing_result = pricerdata->min_cut_pricer->solve(nnodes, nedges, n_cliques, mu_val,
                                                                 node_dual_coeffs, cl_dual_coeffs,
                                                                 adj_lists, adj_sizes,
                                                                 down_branched_vertices, up_branched_vertices, false);


   /* Check if we found solutions with positive reduced cost */
   if( pricing_result.num_solutions > 0 )
   {    
      for( int s = 0; s < pricing_result.num_solutions; ++s )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* vardata;
         int* subset;
         int subsetsize;
         int i;

         // Use pre-built subset from min_cut_pricer
         subset = pricing_result.subsets[s];
         subsetsize = pricing_result.subset_sizes[s];
         
         string name = "vertices";

         // Build variable name from subset
         for( i = 0; i < nnodes; ++i )
         {
            if( subset[i] == 1 )
            {
               name += "_" + std::to_string(i);
            }
         }

         // cout << name << " variable being generated!" << endl;
         // cin.get();

         if( subsetsize < 2 )  /* skip invalid subsets (singletons are generated initially) */
         {
            continue;
         }

         SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, subsetsize, nnodes) );

         /* create variable for a new column with objective function coefficient 0.0 */
         SCIP_CALL( SCIPcreateVarKvertexcut(scip, &newvar, name.c_str(), 0.0, FALSE, FALSE, vardata));

         /* add the new variable to the pricer store */
         SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
         ++pricerdata->nodecnt;


         for( i = 0; i < nnodes; ++i )
         {
            if( subset[i] == 1 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, vertex_cover_constrs[i], newvar, 1.0) );
            }
         }


         for (int c = 0; c < n_cliques; c++)
         {
            for(i = 0; i < nnodes; i++)
            {
                  if((*cliques)[c][i] == true && subset[i] == 1)
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, clique_constrs[c], newvar, 1.0) );
                     break;
                  }
            }
         }

         SCIP_CALL( SCIPaddCoefLinear(scip, alpha_cardinality_constr, newvar, 1.0) );

         SCIPdebug(SCIPprintVar(scip, newvar, NULL) );

         SCIP_CALL(SCIPreleaseVar(scip, &newvar) );

         // cout << "Added variable " << name << endl << " with reduced cost " << overall_val - pi_val << endl;
         // cin.get();

         if(isfarkas){
            pricerdata->n_gen_farkas_cols++;
         }

         if(num_node == 1){
            pricerdata->n_gen_cols_root_node++;
         }      
      }
   }

   /* Clean up subset arrays allocated by min_cut_pricer */
   for (size_t i = 0; i < pricing_result.subsets.size(); ++i) {
      delete[] pricing_result.subsets[i];
   }

   /* Clean up */
   delete[] node_dual_coeffs;

   (*result) = SCIP_SUCCESS;
   
   // write LP for debugging
   // SCIP_CALL( SCIPwriteTransProblem(scip, "afterpricing.lp", "lp", FALSE) );
   // cin.get();
   
   
   return SCIP_OKAY;

}

/**@} */

/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeKvertexcut)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   if( pricerdata != NULL)
   {
      if( pricerdata->min_cut_pricer != NULL )
      {
         delete pricerdata->min_cut_pricer;
      }

      /* free memory */
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->vertex_cover_constrs, pricerdata->nnodes);
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->clique_constrs, pricerdata->n_cliques);

      SCIPfreeBlockMemory(scip, &pricerdata);
   }

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitKvertexcut)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS* cons;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   for( i = 0; i < pricerdata->nnodes; ++i )
   {
      cons = pricerdata->vertex_cover_constrs[i];

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->vertex_cover_constrs[i]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->vertex_cover_constrs[i]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->vertex_cover_constrs[i]) );
   }

   for(int c = 0; c < pricerdata->n_cliques; ++c )
   {
      cons = pricerdata->clique_constrs[c];

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->clique_constrs[c]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->clique_constrs[c]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->clique_constrs[c]) );
   }


   /* get transformed alpha cardinality constraint */
   cons = pricerdata->alpha_cardinality_constr;
   SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->alpha_cardinality_constr) );
   SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->alpha_cardinality_constr) );
   SCIP_CALL( SCIPcaptureCons(scip, pricerdata->alpha_cardinality_constr) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolKvertexcut)
{
   SCIP_PRICERDATA* pricerdata;
   int i;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);


   for( i = 0; i < pricerdata->nnodes; ++i )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->vertex_cover_constrs[i])) );
   }

   for( i = 0; i < pricerdata->n_cliques; ++i )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->clique_constrs[i])) );
   }


   /* release alpha cardinality constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->alpha_cardinality_constr)) );

   return SCIP_OKAY;
}

/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostKvertexcut)
{  /*lint --e{715}*/


   SCIP_CALL( doPricing(scip, pricer, FALSE, lowerbound, stopearly, result) );
   
   return SCIP_OKAY;
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasKvertexcut)
{  /*lint --e{715}*/

   SCIP_CALL( doPricing(scip, pricer, TRUE, NULL, NULL, result) );
   
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the k-vertex-cut variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerKvertexcut(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata = NULL;
   SCIP_PRICER* pricer = NULL;

   /* create k-vertex-cut variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );

   pricerdata->alpha_cardinality_constr = NULL;
   pricerdata->vertex_cover_constrs = NULL;
   pricerdata->clique_constrs = NULL;
   pricerdata->nnodes = 0;
   pricerdata->nedges = 0;
   pricerdata->n_cliques = 0;
   pricerdata->min_cut_pricer = NULL;
   pricerdata->lastnode = NULL;
   pricerdata->nodecnt = 0;
   pricerdata->n_gen_cols_root_node = 0;
   pricerdata->n_gen_farkas_cols = 0;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostKvertexcut, pricerFarkasKvertexcut, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeKvertexcut) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitKvertexcut) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolKvertexcut) );

   return SCIP_OKAY;
}


/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerKvertexcutActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            alpha_cardinality_constr,  /**< alpha cardinality constraint */
   SCIP_CONS**           vertex_cover_constrs,   /**< vertex cover constraints */
   SCIP_CONS**           clique_constrs,     /**< clique constraints */
   int                   nnodes,             /**< number of nodes in the graph */
   int                   nedges,             /**< number of edges in the graph */
   int                   n_cliques           /**< number of cliques */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_PROBDATA* probdata;
   int c;

   assert(scip != NULL);
   assert(alpha_cardinality_constr != NULL);
   assert(vertex_cover_constrs != NULL);
   assert(clique_constrs != NULL);
   assert(nnodes > 0);
   assert(nedges > 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);


   /* copy arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->vertex_cover_constrs, vertex_cover_constrs, nnodes) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->clique_constrs, clique_constrs, n_cliques) );

   pricerdata->nnodes = nnodes;
   pricerdata->nedges = nedges;
   pricerdata->n_cliques = n_cliques;
   pricerdata->alpha_cardinality_constr = alpha_cardinality_constr;

   bool** adj_matr = SCIPprobdataGetAdjMatrix(probdata);
   const std::vector<std::vector<bool>>* cliques = SCIPprobdataGetGraph(probdata)->getCliques();
   pricerdata->min_cut_pricer = new MinCutPricer(nnodes, adj_matr, n_cliques, cliques);

   if( pricerdata->min_cut_pricer == NULL )
   {
      SCIPerrorMessage("Failed to initialize MinCutPricer\n");
      return SCIP_ERROR;
   }
      
   SCIPdebugMsg(scip, "   nnodes: %d nedges: %d  \n", nnodes, nedges);


   for( c = 0; c < nnodes; ++c )
   {
      SCIP_CALL( SCIPcaptureCons(scip, vertex_cover_constrs[c]) );
   }

   for( c = 0; c < n_cliques; ++c )
   {
      SCIP_CALL( SCIPcaptureCons(scip, clique_constrs[c]) );
   }

   SCIP_CALL( SCIPcaptureCons(scip, alpha_cardinality_constr) );
   
   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/** returns number of columns generated at root node */
int SCIPpricerKvertexcutGetNGenColsRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->n_gen_cols_root_node;
}

/** returns number of Farkas columns generated */
int SCIPpricerKvertexcutGetNGenFarkasCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->n_gen_farkas_cols;
}

/**@} */
