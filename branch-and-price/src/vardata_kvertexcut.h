/**@file   vardata_kvertexcut.h
 * @brief  Variable data for k-vertex cut problem containing subset information
 * @author Fabio Ciccarelli
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_KVERTEXCUT__
#define __SCIP_VARDATA_KVERTEXCUT__

#include "scip/scip.h"

/** create variable data for k-vertex cut problem */
SCIP_RETCODE SCIPvardataCreateKvertexcut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
<<<<<<< Updated upstream
   int*                  subset,             /**< array of nodes in the subset S */
   int                   subsetsize         /**< size of the subset S */
=======
   int*                  subset,             /**< array of vertices in the subset S */
   int                   subsetsize,         /**< size of the subset S */
   int                   nnodes              /**< number of vertices in the graph */
>>>>>>> Stashed changes
   );

/** get size of the subset S */
int SCIPvardataGetSubsetSize(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns array of nodes in the subset S */
int* SCIPvardataGetSubset(
   SCIP_VARDATA*         vardata             /**< variable data */
   );


/** check if edge {u, v} is internal to subset S */
SCIP_Bool SCIPvardataIsEdgeInternal(
   SCIP_VARDATA*         vardata,
   int                   u,
   int                   v
   );

/** creates α_S variable for k-vertex cut problem */
SCIP_RETCODE SCIPcreateVarKvertexcut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** prints vardata to file stream */
void SCIPvardataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data */
   FILE*                 file                /**< the text file to store the information into */
   );

#endif