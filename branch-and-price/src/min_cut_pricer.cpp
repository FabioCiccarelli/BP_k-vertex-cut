#include "min_cut_pricer.h"

#include <iostream>
#include <limits> 
#include <set>

// Define LEMON tolerance static member to resolve linking issue
namespace lemon {
    double Tolerance<double>::def_epsilon = 1e-10;
}

MinCutPricer::MinCutPricer(int nnodes, bool** adj_matr, int n_cliques, const std::vector<std::vector<bool>>* cliques)
    : g(), capacity(g) {

    cout << "Building the static graph for MinCutPricer..." << endl;

    const double INF = std::numeric_limits<double>::infinity();

    // --- 1. Create source and sink ---
    s = g.addNode();
    t = g.addNode();

    // --- 2. Create nodes for variables and arcs for costs ---

    // A. Variables x_v 
    x_nodes.resize(nnodes);
    s_arcs_x.resize(nnodes);
    t_arcs_x.resize(nnodes);
    for (int i = 0; i < nnodes; ++i) {
        x_nodes[i] = g.addNode();
        s_arcs_x[i] = g.addArc(s, x_nodes[i]);
        t_arcs_x[i] = g.addArc(x_nodes[i], t);
    }


    // C. Variables cl_c (clique variables)
    cl_nodes.resize(n_cliques);
    s_arcs_cl.resize(n_cliques);
    t_arcs_cl.resize(n_cliques);
    for (int i = 0; i < n_cliques; ++i) {
        cl_nodes[i] = g.addNode();
        s_arcs_cl[i] = g.addArc(s, cl_nodes[i]);
        t_arcs_cl[i] = g.addArc(cl_nodes[i], t);
    }

    // --- 3. ARCS FOR THE STRUCTURAL "PRECEDENCE" CONSTRAINTS (INFINITE CAPACITY) ---
    
    // Arcs for constraints of type x_v <= cl_c
    for (int c = 0; c < n_cliques; ++c)
    {
        for (int v = 0; v < nnodes; ++v)
        {
            if ((*cliques)[c][v] == 1) // Node v is in clique c
            {
                lemon::ListDigraph::Arc arc = g.addArc(x_nodes[v], cl_nodes[c]);
                capacity[arc] = INF;
            }
        }
    }

    // --- 4. ARCS FOR POSSIBLE DOWN-BRANCHING CONSTRAINTS (CAPACITY 0) ---
    for (int u = 0; u < nnodes; ++u) {
        for (int v = u+1; v < nnodes; ++v) {
            if (u == v) continue;
            // check if u and v are adjacent
            if (adj_matr[u][v]) {
                lemon::ListDigraph::Arc arc1 = g.addArc(x_nodes[u], x_nodes[v]);
                branching_arcs[{u, v}] = arc1;
                capacity[arc1] = 0;

                lemon::ListDigraph::Arc arc2 = g.addArc(x_nodes[v], x_nodes[u]);
                branching_arcs[{v, u}] = arc2;
                capacity[arc2] = 0;
            }
        }
    }

    // --- 5. INITIALIZE PREFLOW SOLVER ---
    preflow_solver = new lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<double>>(g, capacity, s, t);

    cout << "Initial graph successfully built." << endl;
}

// Destructor
MinCutPricer::~MinCutPricer() {
    delete preflow_solver;
}


PoolResult MinCutPricer::solve(
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
    bool optpricing) {

    const double INF = std::numeric_limits<double>::infinity();
    PoolResult result;
    result.num_solutions = 0;

    int duplicates = 0;
    
    // Helper set to check for duplicate solutions
    set<vector<int>> seen_solutions;

    // --- 1. SET BASE CAPACITIES ---
    double base_positive_coeffs_sum = 0.0;
    
    // Set capacities for cliques
    for (int i = 0; i < n_cliques; ++i) {
        double coeff = clique_duals[i];
        if (coeff > 0) {
            capacity[s_arcs_cl[i]] = coeff;
            capacity[t_arcs_cl[i]] = 0;
            base_positive_coeffs_sum += coeff;
        } else {
            capacity[s_arcs_cl[i]] = 0;
            capacity[t_arcs_cl[i]] = -coeff;
        }
    }

    // --- 2. ACTIVATE DOWN-BRANCHING CONSTRAINTS ---
    for (int v : down_branched_vertices) {
        for (int i = 0; i < adj_sizes[v]; ++i) {
            int u = adj_lists[v][i];
            capacity[branching_arcs.at({u, v})] = INF;
        }
    }

    // --- 2b. HANDLE UP-BRANCHING CONSTRAINTS ---
    // Create a boolean array to mark up-branched vertices for fast lookup
    vector<bool> is_up_branched(nnodes, false);
    for (int v : up_branched_vertices) {
        is_up_branched[v] = true;
        // Force up-branched vertices to be fixed to 0
        capacity[s_arcs_x[v]] = 0;
        capacity[t_arcs_x[v]] = INF;
    }

    // --- 3. ITERATE OVER ALL NODES, MODIFYING EACH ONE AT A TIME ---
    for (int modified_node = 0; modified_node < nnodes; ++modified_node) {
        // Skip up-branched vertices and preFixed vertices
        if (is_up_branched[modified_node]) {
            continue;
        }
        

        double positive_coeffs_sum = base_positive_coeffs_sum;
        
        // Set base capacities for all nodes
        for (int i = 0; i < nnodes; ++i) {
            // Skip up-branched vertices and preFixed vertices (already set above)
            if (is_up_branched[i]) {
                continue;
            }

            double coeff = node_duals[i];
            if (coeff > 0) {
                capacity[s_arcs_x[i]] = coeff;
                capacity[t_arcs_x[i]] = 0;
                positive_coeffs_sum += coeff;
            } else {
                capacity[s_arcs_x[i]] = 0;
                capacity[t_arcs_x[i]] = -coeff;
            }
        }
        
        // Modify the coefficient of the current node (subtract pi_val)
        double original_coeff = node_duals[modified_node];
        double modified_coeff = node_duals[modified_node] + (mu_val > 0 ? mu_val : 0);
        
        // Update positive_coeffs_sum accounting for the change
        if (original_coeff > 0) {
            positive_coeffs_sum -= original_coeff;
        }
        if (modified_coeff > 0) {
            positive_coeffs_sum += modified_coeff;
        }
        
        // Update capacities for the modified node
        if (modified_coeff > 0) {
            capacity[s_arcs_x[modified_node]] = modified_coeff;
            capacity[t_arcs_x[modified_node]] = 0;
        } else {
            capacity[s_arcs_x[modified_node]] = 0;
            capacity[t_arcs_x[modified_node]] = -modified_coeff;
        }

        // --- 4. SOLVE THE MAX FLOW PROBLEM ---
        preflow_solver->init();
        preflow_solver->run();
        double flow_value = preflow_solver->flowValue();

        double reduced_cost = positive_coeffs_sum - flow_value;
        double threshold = 1e-6 - (mu_val < 0 ? mu_val : 0);

        // --- 5. EXTRACT COLUMN IF REDUCED COST IS POSITIVE ---
        if (reduced_cost > threshold){
            vector<int> new_col_vars;  
            vector<int> new_col_binary(nnodes, 0);  

            // Extract nodes in the source partition
            int num_vertices = 0;
            for (int i = 0; i < nnodes; ++i) {
                // minCut returns true if the node is in the source side of the cut
                if (preflow_solver->minCut(x_nodes[i])) {
                    new_col_vars.push_back(i);  
                    new_col_binary[i] = 1;      
                    num_vertices++;
                }
            }

            // Only add solutions with more than 1 vertex and that are new
            if (num_vertices > 1 && !new_col_vars.empty()) {
                if (seen_solutions.find(new_col_vars) == seen_solutions.end()) {
                    seen_solutions.insert(new_col_vars);
                    result.solutions.push_back(new_col_binary);
                    result.obj_vals.push_back(reduced_cost);
                    
                    // Allocate and populate subset array
                    int* subset = new int[nnodes];
                    for (int i = 0; i < nnodes; ++i) {
                        subset[i] = new_col_binary[i];
                    }
                    result.subsets.push_back(subset);
                    result.subset_sizes.push_back(num_vertices);
                    
                    // cout << "Extracting column with reduced cost: " << reduced_cost << endl;
                    result.num_solutions++;
                }
                else
                {
                    duplicates++;
                }
            }
        }
        
        // Stop if we found enough columns (or duplicates). This is a parameter which can be tuned.
        if(!optpricing)
        {
            if (result.num_solutions >= 10 || duplicates >= 10) {
                break;
            }
        
        }
    }
    
    // --- 6. RESET: DEACTIVATE DOWN-BRANCHING CONSTRAINTS ---
    for (int v : down_branched_vertices) {
        for (int i = 0; i < adj_sizes[v]; ++i) {
            int u = adj_lists[v][i];
            capacity[branching_arcs.at({u, v})] = 0;
        }
    }

    // Reset up-branched vertices capacities
    for (int v : up_branched_vertices) {
        capacity[s_arcs_x[v]] = 0;
        capacity[t_arcs_x[v]] = 0;
    }

      
    return result;
}