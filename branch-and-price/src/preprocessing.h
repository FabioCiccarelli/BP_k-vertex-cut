/**@file   preprocessing.h
 * @brief  Preprocessing routines for k-vertex cut problem
 * @author Fabio Ciccarelli
 */

#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "graph.h"

/**
 * @brief Performs preprocessing on the graph to fix variables
 * 
 * Uses maximum stable set computation to determine which nodes can be fixed.
 * A node v can be fixed if the maximum stable set in its anti-neighborhood is less than k.
 * 
 * @param graph Graph structure
 * @param k Parameter k for k-vertex cut
 * @param preFixed Output array (size nnodes) indicating which nodes are fixed
 * @return Number of fixed nodes
 */
int performPreprocessing(const Graph& graph, int k, bool* preFixed);

/**
 * @brief Computes a feasible k-vertex cut solution using iterative separator extraction
 * 
 * This function iteratively:
 * 1. Computes the minimum weight vertex separator of the current graph
 * 2. Removes the separator vertices
 * 3. Checks if the resulting graph has at least k connected components
 * 4. If not, repeats on the residual graph
 * 
 * The process continues until the graph has at least k connected components.
 * The union of all computed separators forms a feasible k-vertex cut solution.
 * 
 * @param graph Original graph structure
 * @param k Target number of connected components
 * @param solution Output array (size nnodes) indicating which nodes are in the solution
 * @return Value of the computed solution (weight of nodes to remove)
 */
int computeIterativeSeparatorSolution(const Graph& graph, int k, bool* solution);


/**
 * @brief Computes a feasible k-vertex cut solution using an ILP formulation
 * 
 * This function formulates and solves an ILP model to find a k-vertex cut.
 * The ILP uses binary variables to assign vertices to components while ensuring
 * that the removal of selected vertices results in at least k pairwise disjoint
 * connected components.
 * 
 * @param graph Graph structure
 * @param k Target number of connected components
 * @param solution Output array (size nnodes) indicating which nodes are in the solution
 * @return Value of the computed solution (weight of nodes to remove)
 */
int ILPwarmStart(const Graph& graph, int k, bool* solution);

#endif // PREPROCESSING_H
