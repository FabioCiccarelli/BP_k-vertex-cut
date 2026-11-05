/**@file   cons_cutcheck.cpp
 * @brief  constraint handler checking that integral solutions are k-vertex cuts
 * @author Fabio Ciccarelli
 *
 * With relaxed coverage constraints (>=), an integral x can be supported by columns containing cut vertices,
 * so it is not necessarily a k-vertex cut. Integral solutions are therefore checked by counting the connected
 * components of the residual graph without the cut vertices. If an LP solution fails the check, the coverage
 * constraints of the cut vertices that are covered by some column are switched to equalities and the LP is resolved.
 **/

#include <assert.h>

#include "scip/cons_linear.h"

#include "cons_cutcheck.h"
#include "probdata_kvertexcut.h"
#include "graph.h"

#define CONSHDLR_NAME          "cutcheck"
#define CONSHDLR_DESC          "checks that integral solutions are k-vertex cuts"
#define CONSHDLR_ENFOPRIORITY  -2000000
#define CONSHDLR_CHECKPRIORITY -2000000
#define CONSHDLR_EAGERFREQ     -1
#define CONSHDLR_NEEDSCONS     FALSE

struct SCIP_ConshdlrData
{
   int                   nrejected;          /**< number of solutions rejected because they are not k-vertex cuts */
   int                   nfallbacks;         /**< number of times coverage constraints have been switched to equalities */
   int                   nswitchedrows;      /**< number of coverage constraints switched to equalities */
};


/** checks whether the x-part of the solution is integral and, if so, whether it is a k-vertex cut of the residual graph */
static
void checkCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to check, or NULL for the current solution */
   SCIP_Bool*            integral,           /**< pointer to store whether x is integral */
   SCIP_Bool*            iscut               /**< pointer to store whether x is a k-vertex cut */
   )
{
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

   *integral = FALSE;
   *iscut = TRUE;

   if( probdata == NULL )
      return;

   SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);
   Graph* graph = SCIPprobdataGetGraph(probdata);

   /* LP-only mode */
   if( SCIPvarGetType(x_vars[0]) == SCIP_VARTYPE_CONTINUOUS )
      return;

   bool* removed = new bool[graph->nnodes];

   for( int v = 0; v < graph->nnodes; ++v )
   {
      SCIP_Real xval = SCIPgetSolVal(scip, sol, x_vars[v]);

      if( !SCIPisFeasIntegral(scip, xval) )
      {
         delete[] removed;
         return;
      }
      removed[v] = (xval > 0.5);
   }

   *integral = TRUE;
   *iscut = (countConnectedComponents(*graph, removed) >= SCIPprobdataGetk(probdata));

   delete[] removed;
}

/** returns whether all x variables are fixed in the current node */
static
SCIP_Bool allXFixed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);
   int nnodes = SCIPprobdataGetGraph(probdata)->nnodes;

   for( int v = 0; v < nnodes; ++v )
   {
      if( SCIPvarGetLbLocal(x_vars[v]) < 0.5 && SCIPvarGetUbLocal(x_vars[v]) > 0.5 )
         return FALSE;
   }

   return TRUE;
}

/** enforces the cut check for a solution without changing the LP */
static
SCIP_RETCODE enforceSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol,                /**< solution to enforce, or NULL for the current solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_Bool integral;
   SCIP_Bool iscut;

   *result = SCIP_FEASIBLE;

   checkCut(scip, sol, &integral, &iscut);

   if( integral && !iscut )
   {
      SCIPconshdlrGetData(conshdlr)->nrejected++;
      *result = allXFixed(scip) ? SCIP_CUTOFF : SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data */
static
SCIP_DECL_CONSFREE(consFreeCutcheck)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCutcheck)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIP_Bool integral;
   SCIP_Bool iscut;

   *result = SCIP_FEASIBLE;

   checkCut(scip, NULL, &integral, &iscut);

   if( !integral || iscut )
      return SCIP_OKAY;

   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   SCIP_VAR** x_vars = SCIPprobdataGetXVars(probdata);
   SCIP_CONS** coverage_constrs = SCIPprobdataGetVertexCoverConstrs(probdata);
   int nnodes = SCIPprobdataGetGraph(probdata)->nnodes;
   int nswitched = 0;

   for( int v = 0; v < nnodes; ++v )
   {
      SCIP_Real xval = SCIPgetSolVal(scip, NULL, x_vars[v]);

      assert(SCIPconsIsTransformed(coverage_constrs[v]));

      if( xval < 0.5 || !SCIPisInfinity(scip, SCIPgetRhsLinear(scip, coverage_constrs[v])) )
         continue;

      if( SCIPisPositive(scip, SCIPgetActivityLinear(scip, coverage_constrs[v], NULL) - xval) )
      {
         SCIP_CALL( SCIPchgRhsLinear(scip, coverage_constrs[v], 1.0) );
         ++nswitched;
      }
   }

   if( nswitched > 0 )
   {
      conshdlrdata->nfallbacks++;
      conshdlrdata->nswitchedrows += nswitched;

      /* the row changes only mark the LP as unflushed; flushing marks it as unsolved, as required by SCIP_SOLVELP */
      SCIP_CALL( SCIPflushLP(scip) );
      *result = SCIP_SOLVELP;
   }
   else
   {
      conshdlrdata->nrejected++;
      *result = allXFixed(scip) ? SCIP_CUTOFF : SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxCutcheck)
{  /*lint --e{715}*/
   SCIP_CALL( enforceSol(scip, conshdlr, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCutcheck)
{  /*lint --e{715}*/
   SCIP_CALL( enforceSol(scip, conshdlr, NULL, result) );

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
static
SCIP_DECL_CONSCHECK(consCheckCutcheck)
{  /*lint --e{715}*/
   SCIP_Bool integral;
   SCIP_Bool iscut;

   *result = SCIP_FEASIBLE;

   checkCut(scip, sol, &integral, &iscut);

   if( integral && !iscut )
   {
      SCIPconshdlrGetData(conshdlr)->nrejected++;
      *result = SCIP_INFEASIBLE;

      if( printreason )
         SCIPinfoMessage(scip, NULL, "removing the cut vertices leaves fewer than k connected components\n");
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCutcheck)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** creates the cut check constraint handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCutcheck(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->nrejected = 0;
   conshdlrdata->nfallbacks = 0;
   conshdlrdata->nswitchedrows = 0;

   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpCutcheck, consEnfopsCutcheck, consCheckCutcheck, consLockCutcheck, conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeCutcheck) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxCutcheck) );

   return SCIP_OKAY;
}

/** returns the constraint handler data */
static
SCIP_CONSHDLRDATA* getConshdlrData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr);
}

/** returns number of solutions rejected because they are not k-vertex cuts */
int SCIPconshdlrCutcheckGetNRejected(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return getConshdlrData(scip)->nrejected;
}

/** returns number of times the coverage constraints have been switched to equalities */
int SCIPconshdlrCutcheckGetNFallbacks(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return getConshdlrData(scip)->nfallbacks;
}

/** returns number of coverage constraints switched to equalities */
int SCIPconshdlrCutcheckGetNSwitchedRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return getConshdlrData(scip)->nswitchedrows;
}
