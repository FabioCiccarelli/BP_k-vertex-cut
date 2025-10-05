/**@file   prop_cardinality.c
 * @brief  cardinality propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "propagator.h"
#include "probdata_kvertexcut.h"

/* fundamental propagator properties */
#define PROP_NAME              "interdiction"
#define PROP_DESC              "updates bound on beta and zed variables"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/* optional propagator properties */
#define PROP_PRESOL_PRIORITY         -1 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
   int*                  tail;               /**< array of edge tail nodes */
   int*                  head;               /**< array of edge head nodes */
   
   SCIP_VAR**            x_vars;             /**< array of node variables */
   SCIP_VAR**            z_v_vars;           /**< array of z_v variables */
   SCIP_VAR**            z_u_vars;           /**< array of z_u variables */
   SCIP_VAR**            beta_v_vars;        /**< array of beta_v variables */
   SCIP_VAR**            beta_u_vars;        /**< array of beta_u variables */
};


/*
 * Local methods
 */

/** frees data of propagator */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to propdata */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   // SCIPfreeBlockMemoryArray(scip, &propdata->tail, propdata->nedges);
   // SCIPfreeBlockMemoryArray(scip, &propdata->head, propdata->nedges);

   // SCIPfreeBlockMemoryArray(scip, &propdata->x_vars, propdata->nnodes);
   // SCIPfreeBlockMemoryArray(scip, &propdata->z_v_vars, propdata->nedges);
   // SCIPfreeBlockMemoryArray(scip, &propdata->z_u_vars, propdata->nedges);
   // SCIPfreeBlockMemoryArray(scip, &propdata->beta_v_vars, propdata->nedges);
   // SCIPfreeBlockMemoryArray(scip, &propdata->beta_u_vars, propdata->nedges);

   SCIPfreeBlockMemory(scip, propdata);

   return SCIP_OKAY;
}


/** creates data structure of propagator */
static
SCIP_RETCODE propdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to store propagator data */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, propdata) );

   return SCIP_OKAY;
}


/** set data structure of propagator */
static
SCIP_RETCODE propdataSet(
   SCIP_PROPDATA*        propdata,
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   int*                  tail,               /**< array of edge tail nodes */
   int*                  head,               /**< array of edge head nodes */
   SCIP_VAR**            x_vars,             /**< array of node variables */
   SCIP_VAR**            z_v_vars,           /**< array of z_v variables */
   SCIP_VAR**            z_u_vars,           /**< array of z_u variables */
   SCIP_VAR**            beta_v_vars,        /**< array of beta_v variables */
   SCIP_VAR**            beta_u_vars         /**< array of beta_u variables */
   )
{
   assert( propdata != NULL );
   assert( nnodes > 0 );
   assert( nedges > 0 );
   assert( tail != NULL );
   assert( head != NULL );
   assert( x_vars != NULL );
   assert( beta_u_vars != NULL );
   assert( beta_v_vars != NULL );
   assert( z_u_vars != NULL );
   assert( z_v_vars != NULL );

   propdata->nnodes = nnodes;
   propdata->nedges = nedges;

   propdata->tail = tail;
   propdata->head = head;

   propdata->x_vars = x_vars;
   
   propdata->beta_u_vars = beta_u_vars;
   propdata->beta_v_vars = beta_v_vars;
   propdata->z_u_vars = z_u_vars;
   propdata->z_v_vars = z_v_vars;

   return SCIP_OKAY;
}


/** performs the convex hull and special cone fixings */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nreductions,        /**< pointer to store number of rerductions found by propagator */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   SCIP_VAR** x_vars;
   SCIP_VAR** z_v_vars;
   SCIP_VAR** z_u_vars;
   SCIP_VAR** beta_v_vars;
   SCIP_VAR** beta_u_vars;
   SCIP_Bool isfixed;

   int nnodes;
   int nedges;
   int *tail;
   int *head;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( propdata->nnodes > 0 );
   assert( propdata->nedges > 0 );
   assert( propdata->tail != NULL );
   assert( propdata->head != NULL );
   assert( propdata->x_vars != NULL );
   assert( propdata->beta_u_vars != NULL );
   assert( propdata->beta_v_vars != NULL );
   assert( propdata->z_u_vars != NULL );
   assert( propdata->z_v_vars != NULL );

   *nreductions = 0;
   *infeasible = FALSE;

   x_vars = propdata->x_vars;
   z_v_vars = propdata->z_v_vars;
   z_u_vars = propdata->z_u_vars;
   beta_v_vars = propdata->beta_v_vars;
   beta_u_vars = propdata->beta_u_vars;
   nnodes = propdata->nnodes;
   nedges = propdata->nedges;
   tail = propdata->tail;
   head = propdata->head;


   for(int v = 0; v < nnodes; ++v)
   {
      if(SCIPvarGetLbLocal(x_vars[v]) > 0.5) // x_v is fixed to 1
      {
         for(int e = 0; e < nedges; ++e)
         {
            if(tail[e] - 1 == v) // edge e is outgoing from v
            {
               if( SCIPisLE(scip, SCIPvarGetLbLocal(beta_v_vars[e]), 0.0) && SCIPisLE(scip, SCIPvarGetLbLocal(beta_u_vars[e]), 0.0)
                  && SCIPisGE(scip, SCIPvarGetUbLocal(beta_v_vars[e]), 1.0) && SCIPisGE(scip, SCIPvarGetUbLocal(beta_u_vars[e]), 1.0)
                  && SCIPisLE(scip, SCIPvarGetLbLocal(z_v_vars[e]), 0.0) && SCIPisLE(scip, SCIPvarGetLbLocal(z_u_vars[e]), 0.0)
                  && SCIPisGE(scip, SCIPvarGetUbLocal(z_v_vars[e]), 1.0) && SCIPisGE(scip, SCIPvarGetUbLocal(z_u_vars[e]), 1.0))
               {
                  SCIP_CALL( SCIPfixVar(scip, beta_u_vars[e], 1.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;
                  
                  SCIP_CALL( SCIPfixVar(scip, beta_v_vars[e], 0.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;

                  SCIP_CALL( SCIPfixVar(scip, z_u_vars[e], 1.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;
                  
                  SCIP_CALL( SCIPfixVar(scip, z_v_vars[e], 0.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;

                  *nreductions += 4;
               }
            }
            else if(head[e] - 1 == v) // edge e is outgoing from v
            {
               if( SCIPisLE(scip, SCIPvarGetLbLocal(beta_v_vars[e]), 0.0) && SCIPisLE(scip, SCIPvarGetLbLocal(beta_u_vars[e]), 0.0)
                  && SCIPisGE(scip, SCIPvarGetUbLocal(beta_v_vars[e]), 1.0) && SCIPisGE(scip, SCIPvarGetUbLocal(beta_u_vars[e]), 1.0)
                  && SCIPisLE(scip, SCIPvarGetLbLocal(z_v_vars[e]), 0.0) && SCIPisLE(scip, SCIPvarGetLbLocal(z_u_vars[e]), 0.0)
                  && SCIPisGE(scip, SCIPvarGetUbLocal(z_v_vars[e]), 1.0) && SCIPisGE(scip, SCIPvarGetUbLocal(z_u_vars[e]), 1.0))
               {
                  SCIP_CALL( SCIPfixVar(scip, beta_v_vars[e], 1.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;
                  
                  SCIP_CALL( SCIPfixVar(scip, beta_u_vars[e], 0.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;

                  SCIP_CALL( SCIPfixVar(scip, z_v_vars[e], 1.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;
                  
                  SCIP_CALL( SCIPfixVar(scip, z_u_vars[e], 0.0, infeasible, &isfixed) );
                  if( *infeasible )
                     return SCIP_OKAY;

                  *nreductions += 4;
               }
            }
         }
      }
      
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitInterdiction)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
   int*                  tail;               /**< array of edge tail nodes */
   int*                  head;               /**< array of edge head nodes */
   
   SCIP_VAR**            x_vars;             /**< array of node variables */
   SCIP_VAR**            z_v_vars;           /**< array of z_v variables */
   SCIP_VAR**            z_u_vars;           /**< array of z_u variables */
   SCIP_VAR**            beta_v_vars;        /**< array of beta_v variables */
   SCIP_VAR**            beta_u_vars;        /**< array of beta_u variables */

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   nnodes = SCIPprobdataGetNNodes(probdata);
   nedges = SCIPprobdataGetNEdges(probdata);
   tail = SCIPprobdataGetEdgeTails(probdata);
   head = SCIPprobdataGetEdgeHeads(probdata);
   x_vars = SCIPprobdataGetXVars(probdata);
   z_v_vars = SCIPprobdataGetZVVars(probdata);
   z_u_vars = SCIPprobdataGetZUVars(probdata);
   beta_v_vars = SCIPprobdataGetBetaVVars(probdata);
   beta_u_vars = SCIPprobdataGetBetaUVars(probdata);

   /* create convexity propagator data */
   SCIP_CALL( propdataSet(propdata, nnodes, nedges, tail, head, x_vars, z_v_vars, z_u_vars, beta_v_vars, beta_u_vars) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeInterdiction)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( propdataFree(scip, &propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolInterdiction)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nreductions = 0;
   SCIP_Bool infeasible = FALSE;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nreductions, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nreductions > 0 )
   {
      *result = SCIP_SUCCESS;
      *nchgbds += nreductions;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecInterdiction)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   int nfixings = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;


   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nfixings > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropInterdiction)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the cardinality propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropInterdiction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( propdataCreate(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecInterdiction, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitInterdiction) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeInterdiction) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolInterdiction, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropInterdiction) );

   return SCIP_OKAY;
}
