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

/**@file   pricer_kvertexcut.h
 * @brief  k-vertex cut variable pricer
 * @author Fabio Ciccarelli
 **/

 
#ifndef __PRICER_KVERTEXCUT__
#define __PRICER_KVERTEXCUT__

#include "scip/scip.h"
#include "global_variables.h"


/** creates the k-vertex cut variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerKvertexcut(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** activates k-vertex cut pricer */
SCIP_RETCODE SCIPpricerKvertexcutActivate(  
   SCIP*                 scip,                       /**< SCIP data structure */
<<<<<<< Updated upstream
   SCIP_CONS*            main_alpha_constr,          /**< main constraint */
   SCIP_CONS**           alpha_constrs,              /**< edge constraints array */
=======
   SCIP_CONS*            alpha_cardinality_constr,   /**< alpha cardinality constraint */
>>>>>>> Stashed changes
   SCIP_CONS**           coverage_constrs,         /**< coverage constraints array */
   int                   nnodes,                     /**< number of nodes */
   int                   nedges                     /**< number of edges */
   );

#endif
