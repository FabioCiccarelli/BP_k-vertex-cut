/**@file   prop_problemsymmetry.c
 * @brief  problem symmetry propagator
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_problemsymmetry.h"
#include "probdata_kvertexcut.h"

#include "scip/symmetry_orbital.h"
#include "scip/symmetry_lexred.h"

/* fundamental propagator properties */
#define PROP_NAME              "probsym"
#define PROP_DESC              "propagator for orbital and lexicographic reduction"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_ORBITALREDDATA*  orbitalreddata;     /**< data needed for orbital reduction */
   SCIP_LEXREDDATA*      lexreddata;         /**< data needed for lexicographic reduction */
   SCIP_Bool             active;             /**< whether the propagator is active */
};

/*
 * Callback methods of propagator
 */

/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitProblemsymmetry)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->orbitalreddata = SCIPgetProbdataOrbitalreddata(probdata);
   propdata->lexreddata = SCIPgetProbdataLexreddata(probdata);
   if( propdata->orbitalreddata == NULL )
   {
      assert(propdata->lexreddata == NULL);
      propdata->active = FALSE;
   }
   else
      propdata->active = TRUE;

   return SCIP_OKAY;
}

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitProblemsymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecProblemsymmetry)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool didrun1 = FALSE;
   SCIP_Bool didrun2 = FALSE;
   SCIP_Bool infeasible1 = FALSE;
   SCIP_Bool infeasible2 = FALSE;
   int nred1 = 0;
   int nred2 = 0;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !propdata->active || SCIPgetNNodes(scip) <= 1 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPorbitalReductionPropagate(scip, propdata->orbitalreddata, &infeasible1, &nred1, &didrun1) );

   if( infeasible1 )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPlexicographicReductionPropagate(scip, propdata->lexreddata, &infeasible2, &nred2, &didrun2) );

   if( !didrun1 && !didrun2 )
      return SCIP_OKAY;

   if( infeasible2 )
      *result = SCIP_CUTOFF;
   else if( nred1 + nred2 > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the problem symmetry propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropProblemsymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create propagator data */
   propdata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   prop = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecProblemsymmetry, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitProblemsymmetry) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitProblemsymmetry) );

   return SCIP_OKAY;
}
