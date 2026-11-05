/**@file   cons_cutcheck.h
 * @brief  constraint handler checking that integral solutions are k-vertex cuts
 * @author Fabio Ciccarelli
 **/

#ifndef __CONS_CUTCHECK__
#define __CONS_CUTCHECK__

#include "scip/scip.h"

/** creates the cut check constraint handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCutcheck(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns number of solutions rejected because they are not k-vertex cuts */
int SCIPconshdlrCutcheckGetNRejected(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns number of times the coverage constraints have been switched to equalities */
int SCIPconshdlrCutcheckGetNFallbacks(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns number of coverage constraints switched to equalities */
int SCIPconshdlrCutcheckGetNSwitchedRows(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
