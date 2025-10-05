/**
 * @file main.cpp
 * @brief K-vertex cut problem solver using SCIP optimization
 * @author Fabio Ciccarelli
 * @date 2025
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

using namespace std;

/**
 * @brief Structure representing a graph
 * 
 * This structure contains all the necessary information to represent
 * an undirected graph including vertices, edges, and endpoint arrays.
 */
struct Graph {
    int n_vertices;                        ///< Number of vertices in the graph
    int n_edges;                           ///< Number of edges in the graph
    vector<pair<int, int>> edges;          ///< Vector of edge pairs (u,v)
    vector<int> tail;                      ///< Array of edge tails (first endpoint)
    vector<int> head;                      ///< Array of edge heads (second endpoint)
    
    /**
     * @brief Add an edge to the graph
     * @param u First vertex of the edge (0-based indexing)
     * @param v Second vertex of the edge (0-based indexing)
     */
    void addEdge(int u, int v) {
        edges.push_back({u, v});
        tail.push_back(u);
        head.push_back(v);
    }
};

/**
 * @brief Read a graph from DIMACS format file
 * 
 * Parses a graph file in DIMACS format and populates the Graph structure.
 * The DIMACS format consists of:
 * - Comment lines starting with 'c'
 * - Problem line: "p edge <n_vertices> <n_edges>"
 * - Edge lines: "e <u> <v>" where u and v are 1-indexed vertices
 * 
 * @param filename Path to the DIMACS format file
 * @param graph Reference to Graph structure to be populated
 * @return true if file was read successfully, false otherwise
 */
bool readDIMACSGraph(const string& filename, Graph& graph) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return false;
    }
    
    string line;
    bool foundHeader = false;
    
    while (getline(file, line)) {
        if (line.empty() || line[0] == 'c') {
            continue; // Skip comments and empty lines
        }
        
        if (line[0] == 'p') {
            // Header: p edge n_vertices n_edges
            stringstream ss(line);
            string p, edge;
            ss >> p >> edge >> graph.n_vertices >> graph.n_edges;
            
            if (p != "p" || edge != "edge") {
                cerr << "Error: Invalid DIMACS header" << endl;
                return false;
            }
            
            foundHeader = true;
            graph.edges.reserve(graph.n_edges);
            graph.tail.reserve(graph.n_edges);
            graph.head.reserve(graph.n_edges);
        }
        else if (line[0] == 'e') {
            // Edge: e u v
            if (!foundHeader) {
                cerr << "Error: Header not found before edges" << endl;
                return false;
            }
            
            stringstream ss(line);
            string e;
            int u, v;
            ss >> e >> u >> v;
            
            if (e != "e") {
                cerr << "Error: Invalid edge format" << endl;
                return false;
            }
            
            // DIMACS nodes are 1-indexed, convert to 0-indexed
            graph.addEdge(u - 1, v - 1);
        }
    }
    
    file.close();
    
    if (!foundHeader) {
        cerr << "Error: DIMACS header not found" << endl;
        return false;
    }
    
    if (graph.edges.size() != static_cast<size_t>(graph.n_edges)) {
        cerr << "Warning: Number of edges read (" << graph.edges.size() 
             << ") differs from declared (" << graph.n_edges << ")" << endl;
    }
    
    return true;
}

/**
 * @brief Create and solve the k-vertex cut SCIP model
 * 
 * Mathematical model:
 * - Variables: y_vi ∈ {0,1} for each vertex v ∈ V and class i ∈ K = {1,...,k}
 * - Objective: max Σ_{i∈K} Σ_{v∈V} y_vi (maximize total assignments)
 * - Constraints:
 *   1. Σ_{i∈K} y_vi ≤ 1 for each v ∈ V (each vertex assigned to at most one connected component)
 *   2. y_ui + Σ_{j∈K,j≠i} y_vj ≤ 1 for each edge {u,v} ∈ E, each i ∈ K (compatibility)
 *   3. Σ_{v∈V} y_iv ≥ 1 for each i ∈ K (each connected component must contain at least one vertex)
 * 
 * @param graph Input graph structure
 * @param k Number of vertex classes (must be positive)
 * @param time_limit Time limit in seconds (if > 0, otherwise no limit)
 * @param output_file Output file path for results (empty string if no output)
 * @param input_filename Original input filename for logging
 * @return SCIP return code (SCIP_OKAY if successful)
 */
SCIP_RETCODE solveKVertexCutModel(const Graph& graph, int k, double time_limit, 
                                  const string& output_file, const string& input_filename) {
    SCIP* scip = nullptr;
    
    // Initialize SCIP
    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    
    // Create the problem
    SCIP_CALL(SCIPcreateProbBasic(scip, "k-vertex-cut"));
    
    // Set optimization direction (maximization)
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE));
    
    // Set time limit if specified
    if (time_limit > 0) {
        SCIP_CALL(SCIPsetRealParam(scip, "limits/time", time_limit));
        cout << "Time limit set to " << time_limit << " seconds" << endl;
    }
    
    cout << "Creating model for graph with " << graph.n_vertices << " vertices, " 
         << graph.n_edges << " edges and k = " << k << endl;
    
    // Create variables y_vi (v ∈ V, i ∈ K)
    // y[v][i] = 1 if vertex v is assigned to class i, 0 otherwise
    vector<vector<SCIP_VAR*>> y(graph.n_vertices, vector<SCIP_VAR*>(k));
    
    for (int v = 0; v < graph.n_vertices; v++) {
        for (int i = 0; i < k; i++) {
            char varname[100];
            sprintf(varname, "y_%d_%d", v, i);
            
            SCIP_VAR* var;
            // Binary variable with coefficient 1.0 in objective function
            SCIP_CALL(SCIPcreateVarBasic(scip, &var, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
            SCIP_CALL(SCIPaddVar(scip, var));
            
            y[v][i] = var;
        }
    }
    
    // Constraint 1: Σ_{i∈K} y_vi ≤ 1 for each v ∈ V
    // Each vertex can be assigned to at most one class
    for (int v = 0; v < graph.n_vertices; v++) {
        SCIP_CONS* cons;
        char consname[100];
        sprintf(consname, "vertex_assignment_%d", v);
        
        SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons, consname, 0, nullptr, nullptr, -SCIPinfinity(scip), 1.0));
        
        for (int i = 0; i < k; i++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, cons, y[v][i], 1.0));
        }
        
        SCIP_CALL(SCIPaddCons(scip, cons));
        SCIP_CALL(SCIPreleaseCons(scip, &cons));
    }
    
    // Constraint 2: y_ui + Σ_{j∈K,j≠i} y_vj ≤ 1 for each i ∈ K, each edge e = {u,v} ∈ E
    // If vertex u is assigned to class i, then vertex v cannot be assigned to any other class j ≠ i
    // This ensures compatibility between adjacent vertices
    for (int e = 0; e < graph.n_edges; e++) {
        int u = graph.tail[e];
        int v = graph.head[e];
        
        for (int i = 0; i < k; i++) {
            SCIP_CONS* cons;
            char consname[100];
            sprintf(consname, "edge_compatibility_%d_%d_%d_%d", e, i, u, v);
            
            SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons, consname, 0, nullptr, nullptr, -SCIPinfinity(scip), 1.0));
            
            // Add y_ui
            SCIP_CALL(SCIPaddCoefLinear(scip, cons, y[u][i], 1.0));
            
            // Add Σ_{j≠i} y_vj
            for (int j = 0; j < k; j++) {
                if (j != i) {
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, y[v][j], 1.0));
                }
            }
            
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
        }
    }
    
    // Constraint 3: Σ_{v∈V} y_iv ≥ 1 for each i ∈ K
    // Each class must contain at least one vertex (non-empty classes)
    for (int i = 0; i < k; i++) {
        SCIP_CONS* cons;
        char consname[100];
        sprintf(consname, "component_nonempty_%d", i);
        
        SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons, consname, 0, nullptr, nullptr, 1.0, SCIPinfinity(scip)));
        
        for (int v = 0; v < graph.n_vertices; v++) {
            SCIP_CALL(SCIPaddCoefLinear(scip, cons, y[v][i], 1.0));
        }
        
        SCIP_CALL(SCIPaddCons(scip, cons));
        SCIP_CALL(SCIPreleaseCons(scip, &cons));
    }

    // save an LP file for debugging
    // SCIP_CALL( SCIPwriteOrigProblem(scip, "k-vertex-cut.lp", "lp", FALSE) );

    // Set number of threads to 1 for reproducibility
    SCIP_CALL(SCIPsetIntParam(scip, "parallel/maxnthreads", 1));
    
    cout << "Model created. Starting optimization..." << endl;
    
    // Solve the problem
    SCIP_CALL(SCIPsolve(scip));
    
    // Collect solving statistics
    SCIP_STATUS status = SCIPgetStatus(scip);
    double solving_time = SCIPgetSolvingTime(scip);
    SCIP_Longint nnodes_solved = SCIPgetNNodes(scip);
    int nvars = SCIPgetNVars(scip);
    int ncons = SCIPgetNConss(scip);
    
    // Convert status to string
    string status_str;
    switch (status) {
        case SCIP_STATUS_OPTIMAL:
            status_str = "Optimal";
            break;
        case SCIP_STATUS_INFEASIBLE:
            status_str = "Infeasible";
            break;
        case SCIP_STATUS_UNBOUNDED:
            status_str = "Unbounded";
            break;
        case SCIP_STATUS_TIMELIMIT:
            status_str = "TimeLimit";
            break;
        case SCIP_STATUS_BESTSOLLIMIT:
            status_str = "BestSolLimit";
            break;
        default:
            status_str = "OTHER_" + to_string(static_cast<int>(status));
            break;
    }
    
    // Print results to console
    cout << "Solution status: " << status_str << endl;
    cout << "Solving time: " << solving_time << " seconds" << endl;
    cout << "Nodes solved: " << nnodes_solved << endl;
    
    double objval = -1.0;  // Default value if no solution found
    int min_vertex_cut = -1;
    
    if (status == SCIP_STATUS_OPTIMAL || status == SCIP_STATUS_BESTSOLLIMIT || status == SCIP_STATUS_TIMELIMIT) {
        SCIP_SOL* sol = SCIPgetBestSol(scip);
        if (sol != nullptr) {
            objval = SCIPgetSolOrigObj(scip, sol);
            min_vertex_cut = graph.n_vertices - static_cast<int>(objval);
            
            cout << "Objective value: " << objval << endl;
            cout << "Minimum cardinality cut: " << min_vertex_cut << endl;
            cout << "\nVertex-class assignments:" << endl;
            
            for (int v = 0; v < graph.n_vertices; v++) {
                for (int i = 0; i < k; i++) {
                    double val = SCIPgetSolVal(scip, sol, y[v][i]);
                    if (val > 0.5) {  // Binary variable is 1
                        cout << "Vertex " << v + 1 << " -> Class " << i + 1 << endl;
                    }
                }
            }
        }
    }
    
    // Write results to output file if specified
    if (!output_file.empty()) {
        ofstream outfile(output_file, ios::app);  // Append mode
        if (outfile.is_open()) {
            // Check if file is empty (to write header)
            outfile.seekp(0, ios::end);
            bool file_empty = (outfile.tellp() == 0);
            
            if (file_empty) {
                // Write header
                outfile << "input_file\tk\ttime_limit\tnnodes\tnedges\ttime\tstatus\t"
                       << "best_incumbent_value\tmin_vertex_cut\tnvariables\tnconstraints\t"
                       << "branching_tree_nodes" << endl;
            }
            
            // Write data row
            outfile << input_filename << "\t"
                   << k << "\t"
                   << (time_limit > 0 ? to_string(time_limit) : "inf") << "\t"
                   << graph.n_vertices << "\t"
                   << graph.n_edges << "\t"
                   << solving_time << "\t"
                   << status_str << "\t"
                   << SCIPgetDualbound(scip) << "\t"
                   << (objval >= 0 ? to_string(objval) : "NA") << "\t"
                   << (min_vertex_cut >= 0 ? to_string(min_vertex_cut) : "NA") << "\t"
                   << nvars << "\t"
                   << ncons << "\t"
                   << nnodes_solved << endl;
            
            outfile.close();
            cout << "\nResults appended to: " << output_file << endl;
        } else {
            cerr << "Warning: Could not open output file " << output_file << " for writing" << endl;
        }
    }
    
    // Release variables
    for (int v = 0; v < graph.n_vertices; v++) {
        for (int i = 0; i < k; i++) {
            SCIP_CALL(SCIPreleaseVar(scip, &y[v][i]));
        }
    }
    
    // Free SCIP
    SCIP_CALL(SCIPfree(&scip));
    
    return SCIP_OKAY;
}

/**
 * @brief Main function - entry point of the k-vertex cut solver
 * 
 * Parses command line arguments, reads the graph from DIMACS file,
 * and solves the k-vertex cut optimization problem.
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 *             argv[1]: path to DIMACS graph file
 *             argv[2]: number of classes k (positive integer)
 *             argv[3]: time limit in seconds (optional, default: no limit)
 *             argv[4]: output file for results (optional, default: no output file)
 * @return 0 on success, 1 on error
 */
int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        cout << "Usage: " << argv[0] << " <graph_file.dimacs> <k> [time_limit] [output_file]" << endl;
        cout << "  graph_file.dimacs: Path to input graph in DIMACS format" << endl;
        cout << "  k: Number of connected components (positive integer)" << endl;
        cout << "  time_limit: Maximum solving time in seconds (optional, default: no limit)" << endl;
        cout << "  output_file: File to append results (optional, default: no output)" << endl;
        cout << "Example: " << argv[0] << " data/10th_DIMACS/karate.graph.dimacs 3 300 results.txt" << endl;
        return 1;
    }
    
    string filename = argv[1];
    int k = atoi(argv[2]);
    double time_limit = -1.0;  // No time limit by default
    string output_file = "";
    
    if (argc >= 4) {
        time_limit = atof(argv[3]);
        if (time_limit <= 0) {
            cerr << "Error: time_limit must be a positive number" << endl;
            return 1;
        }
    }
    
    if (argc >= 5) {
        output_file = argv[4];
    }
    
    if (k <= 0) {
        cerr << "Error: k must be a positive integer" << endl;
        return 1;
    }
    
    // Read the graph
    Graph graph;
    if (!readDIMACSGraph(filename, graph)) {
        return 1;
    }
    
    cout << "Graph successfully loaded:" << endl;
    cout << "  Vertices: " << graph.n_vertices << endl;
    cout << "  Edges: " << graph.n_edges << endl;
    cout << "  k: " << k << endl;
    if (time_limit > 0) {
        cout << "  Time limit: " << time_limit << " seconds" << endl;
    }
    if (!output_file.empty()) {
        cout << "  Output file: " << output_file << endl;
    }
    cout << endl;
    
    // Solve the model
    SCIP_RETCODE retcode = solveKVertexCutModel(graph, k, time_limit, output_file, filename);
    
    if (retcode != SCIP_OKAY) {
        cerr << "Error in SCIP model resolution" << endl;
        return 1;
    }
    
    return 0;
}
