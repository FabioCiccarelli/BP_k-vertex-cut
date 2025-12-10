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
#include "pricing_algo.h"

#include "global_variables.h"

#include <iostream>
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
<<<<<<< Updated upstream
   SCIP_CONS*            main_alpha_constr;          /**< main pricing constraint for α_S */
   SCIP_CONS**           alpha_constrs;              /**< array of pricing constraints for α_S */
   SCIP_CONS**           coverage_constrs;         /**< array of coverage constraints for SCCs */
=======
   SCIP_CONS*            alpha_cardinality_constr;   /**< cardinality constraint for α variables */
   SCIP_CONS**           vertex_cover_constrs;         /**< array of vertex cover constraints */
   SCIP_CONS**           clique_constrs;           /**< array of clique constraints */
>>>>>>> Stashed changes
   int                   nnodes;                     /**< number of nodes in the graph */
   int                   nedges;                     /**< number of edges in the graph */
   GRBModel*             gurobi_model;               /**< Gurobi model for pricing problem */
   SCIP_NODE* lastnode;
   int nodecnt;
};



/*
 * Data structures
 *//** perform pricing using Gurobi model */
static
SCIP_RETCODE doPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< variable pricer structure */
   SCIP_Bool             isfarkas,           /**< whether to use Farkas pricing */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_PRICERDATA* pricerdata;
<<<<<<< Updated upstream
   SCIP_CONS*            main_alpha_constr;          /**< main pricing constraint for α_S */
   SCIP_CONS**           alpha_constrs;              /**< array of pricing constraints for α_S */
   SCIP_CONS**           coverage_constrs;         /**< array of coverage constraints for SCCs */
=======
   SCIP_CONS*            alpha_cardinality_constr;   /**< cardinality constraint for α variables */
   SCIP_CONS**           vertex_cover_constrs;         /**< array of vertex cover constraints for SCCs */
   SCIP_CONS**           clique_constrs;           /**< array of clique constraints */
>>>>>>> Stashed changes
   SCIP_PROBDATA* probdata;

  
   int nnodes, nedges;
   SCIP_Bool addvar;

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

<<<<<<< Updated upstream
   tail = SCIPprobdataGetEdgeTails(probdata);
   head = SCIPprobdataGetEdgeHeads(probdata);

   nnodes = pricerdata->nnodes;
   nedges = pricerdata->nedges;
   main_alpha_constr = pricerdata->main_alpha_constr;
   alpha_constrs = pricerdata->alpha_constrs;
   coverage_constrs = pricerdata->coverage_constrs;
=======
   int n_cliques = SCIPprobdataGetNCliques(probdata);
   const std::vector<std::vector<bool>>* cliques = SCIPprobdataGetGraph(probdata)->getCliques();

   nnodes = pricerdata->nnodes;
   nedges = pricerdata->nedges;
   alpha_cardinality_constr = pricerdata->alpha_cardinality_constr;
   vertex_cover_constrs = pricerdata->vertex_cover_constrs;
   clique_constrs = pricerdata->clique_constrs;
>>>>>>> Stashed changes

   SCIP_NODE* node = SCIPgetCurrentNode(scip);
   int num_node = SCIPnodeGetNumber(node);


   SCIP_Real pi_val = isfarkas ? SCIPgetDualfarkasLinear(scip, main_alpha_constr) : SCIPgetDualsolLinear(scip, main_alpha_constr);



   // Check if Gurobi model exists
   if( pricerdata->gurobi_model == NULL )
   {
      SCIPwarningMessage(scip, "Gurobi model not initialized in pricer\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

<<<<<<< Updated upstream
   /* Prepare dual values array for edge variables (ψ values) */
   double* edge_dual_coeffs = new double[nedges];
   for( int i = 0; i < nedges; ++i )
   {
      edge_dual_coeffs[i] = isfarkas ? SCIPgetDualfarkasLinear(scip, alpha_constrs[i]) : SCIPgetDualsolLinear(scip, alpha_constrs[i]);
      // cout << "Edge " << i << " dual coeff = " << edge_dual_coeffs[i] << endl;
   }
=======
>>>>>>> Stashed changes

   double* node_dual_coeffs = new double[nnodes];
   for( int i = 0; i < nnodes; ++i )
   {
      node_dual_coeffs[i] = isfarkas ? SCIPgetDualfarkasLinear(scip, vertex_cover_constrs[i]) : SCIPgetDualsolLinear(scip, vertex_cover_constrs[i]);
      // cout << "Node " << i << " dual coeff = " << node_dual_coeffs[i] << endl;
   }


   /* Update objective coefficients in Gurobi model */
   set_coeffs(pricerdata->gurobi_model, pi_val, edge_dual_coeffs, node_dual_coeffs, nedges, nnodes);

   /* Update coefficients of x variables (node variables) with π value */
   for( int v = 0; v < nnodes; ++v )
   {  
      SCIP_Real current_lb = SCIPvarGetLbLocal(x_vars[v]);
      SCIP_Real current_ub = SCIPvarGetUbLocal(x_vars[v]);

      try {
         GRBVar var = pricerdata->gurobi_model->getVarByName("x_" + std::to_string(v));
         var.set(GRB_DoubleAttr_Obj, pi_val);
         
         var.set(GRB_DoubleAttr_LB, 0.0);
         var.set(GRB_DoubleAttr_UB, 1.0);

         /* if variable is fixed to 1 in the current node, set its upper bound to 0 in the pricing problem */
         if( current_lb > 0.5 || preFixed[v])
         {
            var.set(GRB_DoubleAttr_UB, 0.0);
         }
      
      } catch (GRBException &e) {
         SCIPwarningMessage(scip, "Error setting x variable coefficient: %s\n", e.getMessage().c_str());
      }
   }


   /* Solve the pricing problem with Gurobi */
   PoolResult pricing_result = solve_pricing(pricerdata->gurobi_model, nnodes, nedges, pi_val);

<<<<<<< Updated upstream
   addvar = FALSE;
=======
   /* Solve the pricing problem */
   PoolResult pricing_result = pricerdata->min_cut_pricer->solve(nnodes, nedges, n_cliques, mu_val,
                                                                 node_dual_coeffs, cl_dual_coeffs,
                                                                 adj_lists, adj_sizes,
                                                                 down_branched_vertices, up_branched_vertices, false);
>>>>>>> Stashed changes

   // cout << "Pricing result: obj val = " << pricing_result.obj_val << endl;

   /* Check if we found solutions with positive reduced cost */
   if( pricing_result.num_solutions > 0 && pricing_result.obj_val > 1e-9 )
   {      
      for( int s = 0; s < pricing_result.num_solutions; ++s )
      {
         SCIP_VAR* newvar;
         SCIP_VARDATA* vardata;
         int* subset;
         int subsetsize;
         char name[SCIP_MAXSTRLEN];
         char strtmp[SCIP_MAXSTRLEN];
         int i;

         std::vector<double>& sol = pricing_result.solutions[s];

         subsetsize = 0;   

         double overall_val = 0.0;
         
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "vertices");

         SCIP_CALL( SCIPallocBufferArray(scip, &subset, nnodes) );

         for( i = 0; i < nnodes; ++i )
         {
            if( sol[i] > 0.5 )
            {
               subset[i] = 1;
               ++subsetsize;
               (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", i);
               strcat(name, strtmp);
               overall_val += pi_val;

            }
            else
            {
               subset[i] = 0;
            }
         }

<<<<<<< Updated upstream
         if( subsetsize < 2 )  /* skip invalid subsets */
=======
         // cout << name << " variable being generated!" << endl;
         // cin.get();

         if( subsetsize < 2 )  /* skip invalid subsets (singletons are generated initially) */
>>>>>>> Stashed changes
         {
            SCIPfreeBufferArray(scip, &subset);
            continue;
         }

         SCIP_CALL( SCIPvardataCreateKvertexcut(scip, &vardata, subset, subsetsize) );

         /* create variable for a new column with objective function coefficient 0.0 */
         SCIP_CALL( SCIPcreateVarKvertexcut(scip, &newvar, name, 0.0, FALSE, FALSE, vardata));

         /* add the new variable to the pricer store */
         SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
         addvar = TRUE;
         ++pricerdata->nodecnt;

<<<<<<< Updated upstream
         // add the new variable to all the constraints where it appears
         for( i = 0; i < nedges; ++i )
         {
            if( subset[tail[i]-1] == 1 && subset[head[i]-1] == 1 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, alpha_constrs[i], newvar, 1.0) );
               overall_val += edge_dual_coeffs[i];
            }
         }
=======
>>>>>>> Stashed changes

         for( i = 0; i < nnodes; ++i )
         {
            if( subset[i] == 1 )
            {
<<<<<<< Updated upstream
               SCIP_CALL( SCIPaddCoefLinear(scip, coverage_constrs[i], newvar, 1.0) );
               overall_val += node_dual_coeffs[i];
=======
               SCIP_CALL( SCIPaddCoefLinear(scip, vertex_cover_constrs[i], newvar, 1.0) );
>>>>>>> Stashed changes
            }
         }

         SCIP_CALL( SCIPaddCoefLinear(scip, main_alpha_constr, newvar, subsetsize - 1) );

         SCIPdebug(SCIPprintVar(scip, newvar, NULL) );
         SCIP_CALL(SCIPreleaseVar(scip, &newvar) );

         // cout << "Added variable " << name << endl << " with reduced cost " << overall_val - pi_val << endl;
         // cin.get();

         SCIPfreeBufferArray(scip, &subset);

         if(isfarkas){global_NGenFarkasCols++;}

   

         if(num_node == 1){
            global_NGenColsRootNode++;
         }      
      }
   }

   /* Clean up */
   delete[] node_dual_coeffs;

   (*result) = SCIP_SUCCESS;
   
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
      /* free Gurobi model */
      if( pricerdata->gurobi_model != NULL )
      {
         free_LP_pricer(pricerdata->gurobi_model);
         pricerdata->gurobi_model = NULL;
      }

      /* free memory */
<<<<<<< Updated upstream
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->alpha_constrs, pricerdata->nedges);
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->coverage_constrs, pricerdata->nnodes);
=======
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->vertex_cover_constrs, pricerdata->nnodes);
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->clique_constrs, pricerdata->n_cliques);
>>>>>>> Stashed changes

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

   /* get transformed main constraint */
   cons = pricerdata->main_alpha_constr;
   SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->main_alpha_constr) );
   SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->main_alpha_constr) );
   SCIP_CALL( SCIPcaptureCons(scip, pricerdata->main_alpha_constr) );

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

   /* release main constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->main_alpha_constr)) );

   return SCIP_OKAY;
}

/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostKvertexcut)
{  /*lint --e{715}*/
   SCIP_CALL( doPricing(scip, pricer, FALSE, result) );

<<<<<<< Updated upstream
   if( *result == SCIP_DIDNOTRUN )
   {
      *stopearly = TRUE;
   }
  *lowerbound = 1.0;
=======

   SCIP_CALL( doPricing(scip, pricer, FALSE, lowerbound, stopearly, result) );
>>>>>>> Stashed changes
   
   return SCIP_OKAY;
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasKvertexcut)
{  /*lint --e{715}*/

<<<<<<< Updated upstream
   SCIP_CALL( doPricing(scip, pricer, TRUE, result) );
   // (*result) = SCIP_SUCCESS;

=======
   SCIP_CALL( doPricing(scip, pricer, TRUE, NULL, NULL, result) );
   
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
   pricerdata->main_alpha_constr = NULL;
   pricerdata->alpha_constrs = NULL;
   pricerdata->coverage_constrs = NULL;
=======
   pricerdata->alpha_cardinality_constr = NULL;
   pricerdata->vertex_cover_constrs = NULL;
   pricerdata->clique_constrs = NULL;
>>>>>>> Stashed changes
   pricerdata->nnodes = 0;
   pricerdata->nedges = 0;
   pricerdata->gurobi_model = NULL;
   pricerdata->lastnode = NULL;
   pricerdata->nodecnt = 0;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostKvertexcut, pricerFarkasKvertexcut, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeKvertexcut) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitKvertexcut) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolKvertexcut) );

   /* add k-vertex-cut variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}


/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerKvertexcutActivate(
   SCIP*                 scip,               /**< SCIP data structure */
<<<<<<< Updated upstream
   SCIP_CONS*            main_alpha_constr,  /**< main alpha constraint */
   SCIP_CONS**           alpha_constrs,      /**< other alpha constraints */
   SCIP_CONS**           coverage_constrs,   /**< coverage constraints */ 
=======
   SCIP_CONS*            alpha_cardinality_constr,  /**< alpha cardinality constraint */
   SCIP_CONS**           vertex_cover_constrs,   /**< vertex cover constraints */
   SCIP_CONS**           clique_constrs,     /**< clique constraints */
>>>>>>> Stashed changes
   int                   nnodes,             /**< number of nodes in the graph */
   int                   nedges              /**< number of edges in the graph */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   SCIP_PROBDATA* probdata;
   int c;

   assert(scip != NULL);
<<<<<<< Updated upstream
   assert(main_alpha_constr != NULL);
   assert(alpha_constrs != NULL);
   assert(coverage_constrs != NULL);
=======
   assert(alpha_cardinality_constr != NULL);
   assert(vertex_cover_constrs != NULL);
   assert(clique_constrs != NULL);
>>>>>>> Stashed changes
   assert(nnodes > 0);
   assert(nedges > 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);


   /* copy arrays */
<<<<<<< Updated upstream
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->alpha_constrs, alpha_constrs, nedges) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->coverage_constrs, coverage_constrs, nnodes) );

   pricerdata->nnodes = nnodes;
   pricerdata->nedges = nedges;
   pricerdata->main_alpha_constr = main_alpha_constr;
=======
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->vertex_cover_constrs, vertex_cover_constrs, nnodes) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->clique_constrs, clique_constrs, n_cliques) );

   pricerdata->nnodes = nnodes;
   pricerdata->nedges = nedges;
   pricerdata->n_cliques = n_cliques;
   pricerdata->alpha_cardinality_constr = alpha_cardinality_constr;

   bool** adj_matr = SCIPprobdataGetAdjMatrix(probdata);
   const std::vector<std::vector<bool>>* cliques = SCIPprobdataGetGraph(probdata)->getCliques();
   pricerdata->min_cut_pricer = new MinCutPricer(nnodes, adj_matr, n_cliques, cliques);
>>>>>>> Stashed changes

   /* Initialize Gurobi model for pricing problem */
   pricerdata->gurobi_model = init_LP_pricer(nedges, nnodes, tail, head);
   
   if( pricerdata->gurobi_model == NULL )
   {
      SCIPerrorMessage("Failed to initialize Gurobi model for pricing problem\n");
      return SCIP_ERROR;
   }

   SCIPdebugMsg(scip, "   nnodes: %d nedges: %d  \n", nnodes, nedges);
   SCIPdebugMsg(scip, "   Gurobi model initialized successfully\n");


   for( c = 0; c < nnodes; ++c )
   {
      SCIP_CALL( SCIPcaptureCons(scip, vertex_cover_constrs[c]) );
   }

   SCIP_CALL( SCIPcaptureCons(scip, main_alpha_constr) );
   
   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/**@} */
