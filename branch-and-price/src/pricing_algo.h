#ifndef PRICING_ALGO_H
#define PRICING_ALGO_H

#include "gurobi_c++.h"
#include <vector>
#include <numeric> 
#include <algorithm>

/**
 * @brief Result structure for optimization with pooling
 */
struct PoolResult {
    double obj_val;                        ///< Objective value (-1 if no solution found)
    int num_solutions;                     ///< Number of solutions found
    std::vector<std::vector<double>> solutions;  ///< Solutions (x_v values for each node)
};

/**
 * @brief Initialize the LP model for the pricing problem
 * 
 * Creates a Gurobi model with:
 * - Variables x_v for each node (0 <= x_v <= 1, continuous)
 * - Variables y_e for each edge (0 <= y_e, continuous)
 * - Constraints: y_e <= x_u and y_e <= x_v for each edge e=(u,v)
 * 
 * The model is set to maximization mode with output suppressed.
 * Initial objective coefficients: x_v have coefficient -1.0, y_e have coefficient 0.0
 * 
 * @param nedges Number of edges in the graph
 * @param nnodes Number of nodes in the graph
 * @param tail Array of tail nodes for each edge
 * @param head Array of head nodes for each edge
 * @return Pointer to the initialized GRBModel, or nullptr on error
 */
GRBModel* init_LP_pricer(int nedges, int nnodes, int* tail, int* head);

/**
 * @brief Set objective coefficients for edge variables (y_e)
 * 
 * Updates the objective coefficients of all y_e variables based on the 
 * dual values from the master problem. This function should be called
 * at each pricing iteration with updated dual values.
 * 
 * @param model Pointer to the Gurobi model
 * @param coeffs Array of coefficients for edge variables y_e
 * @param nedges Number of edges (length of coeffs array)
 */
void set_coeffs(GRBModel* model, double pi_val, double* edge_coeffs, double* node_coeffs, int nedges, int nnodes);

/**
 * @brief Optimize the model and handle different solution scenarios
 * 
 * This function:
 * 1. First tries to solve the model as is
 * 2. If objective > 1e-6, returns the solution
 * 3. If objective ≈ 0, tries fixing each x_v coefficient to 0 and re-solving
 * 4. Returns the first solution found with positive objective
 * 
 * @param model Pointer to the Gurobi model
 * @param nnodes Number of nodes in the graph
 * @param nedges Number of edges in the graph
 * @return PoolResult containing objective value, number of solutions, and solution vectors
 */
PoolResult solve_pricing(GRBModel* model, int nnodes, int nedges, double pi_val);

/**
 * @brief Free memory allocated for the LP model
 * 
 * Properly deallocates the Gurobi model. The associated environment
 * is handled automatically by Gurobi.
 * 
 * @param model Pointer to the Gurobi model to be freed
 */
void free_LP_pricer(GRBModel* model);

#endif // PRICING_ALGO_H