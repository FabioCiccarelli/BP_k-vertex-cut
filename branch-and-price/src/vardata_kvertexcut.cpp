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

/**@file   vardata_kvertexcut.cpp
 * @brief  Variable data for k-vertex cut problem containing subset information
 * @author Fabio Ciccarelli
 *
 * This file implements the handling of the variable data which is attached to each α_S variable
 * in the k-vertex cut problem. Each variable represents a subset S and stores information
 * about which nodes are in the subset and which constraints it appears in.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "probdata_kvertexcut.h"
#include "vardata_kvertexcut.h"

/** Variable data which is attached to each α_S variable.
 *
 *  This variable data stores information about the subset S that the variable represents,
 *  including which nodes are in the subset and in which constraints this variable appears.
 */
struct SCIP_VarData
{
   int*                  subset;             /**< array of nodes in the subset S */
   int                   subsetsize;         /**< size of the subset S */
};

/**@name Local methods
 *
 * @{
 */

/** create a vardata */
static
SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  subset,             /**< array of nodes in the subset S */
   int                   subsetsize         /**< size of the subset S */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, vardata) );

   /* copy subset array */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*vardata)->subset, subset, subsetsize) );
   //SCIPsortInt((*vardata)->subset, subsetsize);
   (*vardata)->subsetsize = subsetsize;

   return SCIP_OKAY;
}

/** frees user data of variable */
static
SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   SCIPfreeBlockMemoryArray(scip, &(*vardata)->subset, (*vardata)->subsetsize);
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}


/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees user data of transformed variable (called when the transformed variable is freed) */
static
SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
}/*lint !e715*/

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** create variable data for k-vertex cut problem */
SCIP_RETCODE SCIPvardataCreateKvertexcut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  subset,             /**< array of nodes in the subset S */
   int                   subsetsize        /**< size of the subset S */
   )
{
   SCIP_CALL( vardataCreate(scip, vardata, subset, subsetsize) );

   return SCIP_OKAY;
}

/** get size of the subset S */
int SCIPvardataGetSubsetSize(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   return vardata->subsetsize;
}

/** returns array of nodes in the subset S */
int* SCIPvardataGetSubset(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   return vardata->subset;
}


/** check if edge {u, v} is internal to subset S */
SCIP_Bool SCIPvardataIsEdgeInternal(
   SCIP_VARDATA*         vardata,
   int                   u,
   int                   v
   )
{
   assert(vardata != NULL);
   assert(vardata->subsetsize > 0);

   for(int i = 0; i < vardata->subsetsize; i++) {
       if(vardata->subset[i] == u) {
           for(int j = i+1; j < vardata->subsetsize; j++) {
               if(vardata->subset[j] == v) {
                   return TRUE; // both u and v are in the subset
               }
           }
       }
   }

   return FALSE;
}

/** creates α_S variable for k-vertex cut problem */
SCIP_RETCODE SCIPcreateVarKvertexcut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   /* create a basic variable object - α_S variables are continuous, non-negative */
   SCIP_CALL( SCIPcreateVarBasic(scip, var, name, 0.0, SCIPinfinity(scip), obj, SCIP_VARTYPE_CONTINUOUS) );
   assert(*var != NULL);

   /* set callback functions */
   SCIPvarSetData(*var, vardata);
   SCIPvarSetDeltransData(*var, vardataDelTrans);

   /* set initial and removable flag */
   SCIP_CALL( SCIPvarSetInitial(*var, initial) );
   SCIP_CALL( SCIPvarSetRemovable(*var, removable) );

   SCIPvarMarkDeletable(*var);

   SCIPdebug( SCIPprintVar(scip, *var, NULL) );

   return SCIP_OKAY;
}

/** prints vardata to file stream */
void SCIPvardataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA*         vardata,            /**< variable data */
   FILE*                 file                /**< the text file to store the information into */
   )
{
   int i;

   assert(vardata != NULL);

   SCIPinfoMessage(scip, file, "subset S = {");

   for( i = 0; i < vardata->subsetsize; ++i )
   {
      SCIPinfoMessage(scip, file, "%d", vardata->subset[i]);

      if( i < vardata->subsetsize - 1 )
         SCIPinfoMessage(scip, file, ",");
   }

   SCIPinfoMessage(scip, file, "}, |S| = %d\n", vardata->subsetsize);

}

/**@} */
