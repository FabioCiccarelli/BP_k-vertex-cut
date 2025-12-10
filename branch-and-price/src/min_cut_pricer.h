// min_cut_pricer.h
#pragma once

#include <vector>      
#include <map>        
#include <utility>    

// LEMON libraries
#include <lemon/list_graph.h>  
#include <lemon/preflow.h>    

using namespace std;

/**
 * @brief Result structure for optimization with pooling
 */
struct PoolResult {
    std::vector<double> obj_vals;              ///< Objective values for each solution
    int num_solutions;                         ///< Number of solutions found
    std::vector<std::vector<int>> solutions;   ///< Binary solutions (0/1 for each node)
    std::vector<int*> subsets;                 ///< Array of subset pointers (for SCIP)
    std::vector<int> subset_sizes;             ///< Size of each subset
};

class MinCutPricer {
public:
    // Builder
    MinCutPricer(int nnodes, bool** adj_matr, int n_cliques, const std::vector<std::vector<bool>>* cliques);
    
    // Destructor
    ~MinCutPricer();

    // Main pricing function
    PoolResult solve(
        int nnodes,
        int nedges,
        int n_cliques,
        double mu_val,
        double* node_duals,
        double* clique_duals,
        int** adj_lists, 
        int* adj_sizes,
        vector<int> down_branched_vertices,
        vector<int> up_branched_vertices,
        bool onlyone
    );

private:
    // --- LEMON GRAPH ---
    lemon::ListDigraph g;                           // static graph
    lemon::ListDigraph::Node s, t;                  // source and sink nodes
    lemon::ListDigraph::ArcMap<double> capacity;    // capacity map (the only one that changes dynamically)


    // 1. vertices of the network (variables x_v, y_e, gamma_c)
    vector<lemon::ListDigraph::Node> x_nodes;  // Map index v -> Node object of x_v
    vector<lemon::ListDigraph::Node> cl_nodes; // Map index c -> Node object of gamma_c

    // 2. For the cost arcs (those that change most often)
    vector<lemon::ListDigraph::Arc> s_arcs_x;  // Arcs s -> x_v
    vector<lemon::ListDigraph::Arc> t_arcs_x;  // Arcs x_v -> t
    vector<lemon::ListDigraph::Arc> s_arcs_cl; // Arcs s -> gamma_c
    vector<lemon::ListDigraph::Arc> t_arcs_cl; // Arcs gamma_c -> t

    // 3. for the down-branching arcs (x_u <= x_v)
    // the key is the pair of indices (u, v),
    // the value is the corresponding LEMON arc handle.
    map<pair<int, int>, lemon::ListDigraph::Arc> branching_arcs;
    
    // 4. Preflow solver
    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<double>>* preflow_solver;
};