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

/**@file   probdata_kvertexcut.h
 * @brief  Problem data for k-vertex cut problem
 * @author Fabio Ciccarelli
 *
 * This file handles the main problem data used in that project. For more details see \ref KVERTEXCUT_PROBLEMDATA page.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_KVERTEXCUT__
#define __SCIP_PROBDATA_KVERTEXCUT__

#include "scip/scip.h"
#include "global_variables.h"

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   int                   nnodes,             /**< number of nodes in the graph */
   int                   nedges,             /**< number of edges in the graph */
   int*                  tail,              /**< array of edge tail nodes */
   int*                  head,              /**< array of edge head nodes */
   int                   k,                   /**< parameter k for k-vertex cut */
   int*                  vertex_weights     /**< array of vertex weights */
   );

/** returns number of nodes */
int SCIPprobdataGetNNodes(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of edges */
int SCIPprobdataGetNEdges(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of edge tail nodes */
int* SCIPprobdataGetEdgeTails(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of edge head nodes */
int* SCIPprobdataGetEdgeHeads(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns adjacency list for a given node */
int* SCIPprobdataGetAdjList(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   );

/** returns size of adjacency list for a given node */
int SCIPprobdataGetAdjSize(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   );

/** returns degree of a given node */
int SCIPprobdataGetNodeDegree(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node                /**< node index */
   );

/** returns parameter k */
int SCIPprobdataGetk(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

<<<<<<< Updated upstream
=======
/** returns graph structure */
Graph* SCIPprobdataGetGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns graph structure */
Graph* SCIPprobdataGetOrigGraph(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

>>>>>>> Stashed changes
/** returns adjacency matrix */
SCIP_Bool** SCIPprobdataGetAdjMatrix(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** checks if two nodes are adjacent */
SCIP_Bool SCIPprobdataAreAdjacent(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   node1,              /**< first node */
   int                   node2               /**< second node */
   );

/** returns array of all variables ordered in the way they got generated */
SCIP_VAR** SCIPprobdataGetAlphaVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of x variables (node variables) */
SCIP_VAR** SCIPprobdataGetXVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

SCIP_VAR** SCIPprobdataGetZUVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

<<<<<<< Updated upstream
SCIP_VAR** SCIPprobdataGetZVVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

SCIP_VAR** SCIPprobdataGetBetaUVars(
=======
/** returns array of clique cuts */
SCIP_ROW** SCIPprobdataGetCliqueCuts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/**returns number of clique cuts */
int SCIPprobdataGetNCliqueCuts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


/** returns array of coverage constraints */
SCIP_CONS** SCIPprobdataGetVertexCoverConstrs(
>>>>>>> Stashed changes
   SCIP_PROBDATA*        probdata            /**< problem data */
   );


SCIP_VAR** SCIPprobdataGetBetaVVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns number of variables */
int SCIPprobdataGetNAlphaVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of set pricing constraints */
SCIP_CONS** SCIPprobdataGetAlphaConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns array of coverage constraints */
SCIP_CONS** SCIPprobdataGetCoverageConstrs(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** returns main pricing constraint */
SCIP_CONS* SCIPprobdataGetMainAlphaConstr(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   );

SCIP_RETCODE SCIPprobdataAddCliqueCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_ROW*             cut,              /**< cut to add */
   int*                  clique             /**< clique associated with the cut */
   );

/** returns solution value as cut array and updates total cut cost */
bool* SCIPprobdataGetSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int&                  tot_cut_cost        /**< reference to store total cut cost */
   );

/** plots the solution to a PNG file */
SCIP_RETCODE SCIPprobdataPlotSolution(
   Graph*                graph,           /**< instance graph */
   bool*                 cut,                /**< cut array */
   const char*           filename,            /**< output filename */
   bool                  draw_attacked_edges,       /**< whether to highlight cut nodes */
   bool                  plot_graph            /**< whether to plot the graph */
   );

#endif
