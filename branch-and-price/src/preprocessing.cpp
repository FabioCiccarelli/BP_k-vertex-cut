/**@file   preprocessing.cpp
 * @brief  Preprocessing routines for k-vertex cut problem
 * @author Fabio Ciccarelli
 */

#include "preprocessing.h"
#include "gurobi_c++.h"
#include <iostream>
#include <cstring>

using namespace std;

/**
 * @brief Computes maximum stable set in the subgraph induced by neighbors of node v
 * 
 * @param env Gurobi environment
 * @param graph Graph structure
 * @param v Node whose neighborhood we consider
 * @return Size of maximum stable set in neighborhood of v
 */
static
int maxStableSetInNeighborhood(GRBEnv& env, const Graph& graph, int v)
{
    int output = 0;
    int nnodes = graph.nnodes;
    int nedges = graph.nedges;
    
    try {
        // Create model
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_LogToConsole, 0);
        model.set(GRB_IntParam_Threads, 1);
        
        // Get neighbors of v
        int num_neighbors;
        const int* neighbors = graph.getNeighbors(v, num_neighbors);
        
        // Create binary variable for each neighbor
        GRBVar* x = new GRBVar[nnodes];
        for (int i = 0; i < nnodes; i++) {
            x[i] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
        }
        
        // Set bounds: neighbors cannot be selected
        for (int i = 0; i < num_neighbors; i++) {
            x[neighbors[i]].set(GRB_DoubleAttr_UB, 0.0);
        }

        x[v].set(GRB_DoubleAttr_UB, 0.0);
        
        // Add edge constraints: for each edge, at most one endpoint can be selected
        for (int e = 0; e < nedges; e++) {
            int u = graph.tail[e];
            int w = graph.head[e];
            
            model.addConstr(x[u] + x[w] <= 1);
        }
        
        // Objective: maximize sum of selected nodes
        GRBLinExpr obj = 0;
        for (int i = 0; i < nnodes; i++) {
            obj += x[i];
        }
        model.setObjective(obj, GRB_MAXIMIZE);
        
        // Optimize
        model.optimize();
        
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            output = (int)(model.get(GRB_DoubleAttr_ObjVal) + 0.5);
        }
        
        delete[] x;
        
    } catch(GRBException& e) {
        cerr << "Gurobi exception in maxStableSetInNeighborhood: " << e.getMessage() << endl;
    } catch(...) {
        cerr << "Unknown exception in maxStableSetInNeighborhood" << endl;
    }
    
    return output;
}

int performPreprocessing(const Graph& graph, int k, bool* preFixed)
{
    int nnodes = graph.nnodes;
    
    cout << "\n=== Preprocessing ===" << endl;
    cout << "Computing bounds for " << nnodes << " nodes with k=" << k << endl;
    
    // Initialize preFixed array
    memset(preFixed, false, nnodes * sizeof(bool));
    
    int num_fixed = 0;
    
    // Create Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);
    
    // Compute bound for each node
    int* bound = new int[nnodes];
    for (int i = 0; i < nnodes; i++) {
        bound[i] = maxStableSetInNeighborhood(env, graph, i);
    }
    
    // Fix nodes whose bound is less than k
    int v = 0;
    while (v < nnodes) {
        if (k >= bound[v]+2 && !preFixed[v]) {
            preFixed[v] = true;
            num_fixed++;
            for(int u=0; u<nnodes; u++){
				if(!preFixed[u])
					bound[u]=maxStableSetInNeighborhood(env, graph, u);
			}
            v = 0;
        }
        else v++;
    }
    
    delete[] bound;
    
    return num_fixed;
}

int computeIterativeSeparatorSolution(const Graph& graph, int k, bool* solution)
{
    int nnodes = graph.nnodes;
    int overall_sol_cost = 0;
    
    cout << "\n=== Computing heuristic solution ===" << endl;
    cout << "Target: k = " << k << " connected components" << endl;
    
    // Initialize solution array (no nodes removed initially)
    memset(solution, 0, nnodes * sizeof(bool));
    int solution_size = 0;
    
    // Track which nodes are still in the current graph
    bool* removed = new bool[nnodes];
    memset(removed, 0, nnodes * sizeof(bool));
    
    // Count initial connected components
    int num_components = countConnectedComponents(graph, removed);
    cout << "Initial graph has " << num_components << " connected component(s)" << endl;
    
    if (num_components >= k) {
        cout << "Graph already has >= k components, no solution needed!" << endl;
        delete[] removed;
        cout << "==============================================\n" << endl;
        return 0;
    }
    
    // Iteratively compute separators
    int iteration = 0;
    Graph* current_graph = nullptr;
    vector<int> old_to_new, new_to_old;
    
    while (num_components < k) {
        iteration++;
        cout << "\n--- Iteration " << iteration << " ---" << endl;
        cout << "Current components: " << num_components << ", target: " << k << endl;
        
        // Build residual graph with nodes not yet removed
        if (current_graph != nullptr) {
            delete current_graph;
        }
        current_graph = buildResidualGraph(graph, removed, old_to_new, new_to_old);
        
        // Compute minimum separator on current graph
        bool* separator_residual = new bool[current_graph->nnodes];
        int connectivity = computeVertexConnectivity(*current_graph, true, separator_residual);
        
        cout << "Connectivity of residual graph: " << connectivity << endl;

        overall_sol_cost += connectivity;
        
        // Count separator size
        int separator_size = 0;
        for (int i = 0; i < current_graph->nnodes; i++) {
            if (separator_residual[i]) {
                separator_size++;
            }
        }
        cout << "Separator size: " << separator_size << endl;
        
        if (separator_size == 0) {
            cout << "WARNING: Empty separator found! Graph may be complete or disconnected." << endl;
            delete[] separator_residual;
            break;
        }
        
        // Map separator back to original graph indices and add to solution
        cout << "Adding nodes to solution: {";
        bool first = true;
        for (int i = 0; i < current_graph->nnodes; i++) {
            if (separator_residual[i]) {
                int original_idx = new_to_old[i];
                if (!solution[original_idx]) {
                    solution[original_idx] = true;
                    removed[original_idx] = true;
                    solution_size++;
                    if (!first) cout << ", ";
                    cout << original_idx;
                    first = false;
                }
            }
        }
        cout << "}" << endl;
        
        delete[] separator_residual;
        
        // Count components after removing separator
        num_components = countConnectedComponents(graph, removed);
        cout << "After removing separator: " << num_components << " component(s)" << endl;
        cout << "Overall solution cost so far: " << overall_sol_cost << endl;
        
        // Safety check: if we removed all nodes or separator didn't increase components
        int remaining_nodes = 0;
        for (int i = 0; i < nnodes; i++) {
            if (!removed[i]) remaining_nodes++;
        }
        
        if (remaining_nodes == 0) {
            cout << "WARNING: All nodes removed!" << endl;
            break;
        }
        
        if (num_components >= k) {
            cout << "Target achieved! Graph now has " << num_components << " >= " << k << " components" << endl;
            break;
        }
        
        // Safety limit on iterations
        if (iteration >= nnodes) {
            cout << "WARNING: Too many iterations (" << iteration << "), stopping." << endl;
            break;
        }
    }
    
    // Clean up
    if (current_graph != nullptr) {
        delete current_graph;
    }
    delete[] removed;
    
    cout << "\n=== Solution Summary ===" << endl;
    cout << "Total iterations: " << iteration << endl;
    cout << "Final number of components: " << num_components << endl;
    cout << "Solution size: " << solution_size << " nodes" << endl;
    cout << "Nodes in solution: {";
    bool first = true;
    for (int i = 0; i < nnodes; i++) {
        if (solution[i]) {
            if (!first) cout << ", ";
            cout << i;
            first = false;
        }
    }
    cout << "}" << "\n\n";
    
    
    return solution_size;
}



int ILPwarmStart(const Graph& graph, int k, bool* solution)
{
    int nnodes = graph.nnodes;
    int nedges = graph.nedges;
    int cut_cost = 0;
    
    // Initialize solution array (no nodes in cut initially)
    memset(solution, false, nnodes * sizeof(bool));
    
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);
    
    try {
        // Create model
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_LogToConsole, 0);
        model.set(GRB_IntParam_Threads, 1);
        model.set(GRB_DoubleParam_TimeLimit, 20.0); 
        
        // Create variables y_vi (v ∈ V, i ∈ K)
        // y[v][i] = 1 if vertex v is assigned to class i, 0 otherwise
        vector<vector<GRBVar>> y(nnodes, vector<GRBVar>(k));
        
        for (int v = 0; v < nnodes; v++) {
            for (int i = 0; i < k; i++) {
                char varname[100];
                sprintf(varname, "y_%d_%d", v, i);
                // Binary variable with coefficient 1.0 in objective (maximize assignment)
                y[v][i] = model.addVar(0.0, 1.0, graph.vertex_weights[v], GRB_BINARY, varname);
            }
        }
        
        // Constraint 1: Σ_{i∈K} y_vi ≤ 1 for each v ∈ V
        // Each vertex can be assigned to at most one class
        for (int v = 0; v < nnodes; v++) {
            GRBLinExpr expr = 0;
            for (int i = 0; i < k; i++) {
                expr += y[v][i];
            }
            char consname[100];
            sprintf(consname, "vertex_assignment_%d", v);
            model.addConstr(expr <= 1, consname);
        }
        
        // Constraint 2: y_ui + Σ_{j∈K,j≠i} y_vj ≤ 1 for each i ∈ K, each edge e = {u,v} ∈ E
        // If vertex u is assigned to class i, then vertex v cannot be assigned to any other class j ≠ i
        for (int e = 0; e < nedges; e++) {
            int u = graph.tail[e];
            int v = graph.head[e];
            
            for (int i = 0; i < k; i++) {
                GRBLinExpr expr = 0;
                
                // Add y_ui
                expr += y[u][i];
                
                // Add Σ_{j≠i} y_vj
                for (int j = 0; j < k; j++) {
                    if (j != i) {
                        expr += y[v][j];
                    }
                }
                
                char consname[100];
                sprintf(consname, "edge_compatibility_%d_%d_%d_%d", e, i, u, v);
                model.addConstr(expr <= 1, consname);
            }
        }
        
        // Constraint 3: Σ_{v∈V} y_vi ≥ 1 for each i ∈ K
        // Each class must contain at least one vertex (non-empty classes)
        for (int i = 0; i < k; i++) {
            GRBLinExpr expr = 0;
            for (int v = 0; v < nnodes; v++) {
                expr += y[v][i];
            }
            char consname[100];
            sprintf(consname, "component_nonempty_%d", i);
            model.addConstr(expr >= 1, consname);
        }
        
        // Objective: maximize the sum of all y_vi (equivalently minimize cut vertices)
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        
        // Optimize
        model.optimize();
        
        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            // Extract solution: vertices with no assignment (all y[v][i] = 0) are in the cut
            for (int v = 0; v < nnodes; v++) {
                bool assigned = false;
                for (int i = 0; i < k; i++) {
                    if (y[v][i].get(GRB_DoubleAttr_X) > 0.5) {
                        assigned = true;
                        break;
                    }
                }
                if (!assigned) {
                    solution[v] = true;
                    cut_cost += graph.vertex_weights[v];
                }
            }
            
            cout << "\n\nILP Warm Start: found solution with cut cost = " << cut_cost << " in " << model.get(GRB_DoubleAttr_Runtime) << " seconds\n\n";
        } else {
            cout << "ILP Warm Start: no feasible solution found (status = " << status << ")" << endl;
        }
        
    } catch(GRBException& e) {
        cerr << "Gurobi exception in ILPwarmStart: " << e.getMessage() << endl;
    } catch(...) {
        cerr << "Unknown exception in ILPwarmStart" << endl;
    }

    return cut_cost;
}
