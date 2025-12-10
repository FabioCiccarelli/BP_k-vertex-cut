/**@file   graph.cpp
 * @brief  Graph data structure and related functions for k-vertex cut problem
 * @author Fabio Ciccarelli
 */

#include "graph.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <stack>
#include <algorithm>
#include <limits>
#include <climits>

// LEMON library
#include <lemon/list_graph.h>
#include <lemon/preflow.h>

using namespace std;

// ============================================================================
// Graph Constructor/Destructor and Memory Management
// ============================================================================

Graph::Graph() 
    : nnodes(0), nedges(0), tail(nullptr), head(nullptr), 
      adj_matrix(nullptr), adj_matrix_compl(nullptr), adj_lists(nullptr), node_degree(nullptr), 
      vertex_weights(nullptr), cliques(nullptr)
{
}

Graph::Graph(int n_nodes, int n_edges, int* edge_tails, int* edge_heads, int* v_weights, bool copy_edges)
    : nnodes(n_nodes), nedges(n_edges)
{
    // Handle edge arrays
    if (copy_edges) {
        tail = new int[nedges];
        head = new int[nedges];
        memcpy(tail, edge_tails, nedges * sizeof(int));
        memcpy(head, edge_heads, nedges * sizeof(int));
    } else {
        tail = edge_tails;
        head = edge_heads;
    }
    
    if (copy_edges) {
        this->vertex_weights = new int[nnodes];
        memcpy(this->vertex_weights, v_weights, nnodes * sizeof(int));
    } else {
        this->vertex_weights = v_weights;
    }

    // Initialize other arrays
    adj_matrix = nullptr;
    adj_matrix_compl = nullptr;
    adj_lists = nullptr;
    node_degree = nullptr;
    cliques = nullptr;
    
    // Build data structures
    buildAdjacencyMatrix();
    buildAdjacencyLists();
}

Graph::~Graph()
{
    free();
}

Graph::Graph(const Graph& other)
    : nnodes(0), nedges(0), tail(nullptr), head(nullptr),
      adj_matrix(nullptr), adj_matrix_compl(nullptr), adj_lists(nullptr), node_degree(nullptr),
      vertex_weights(nullptr), cliques(nullptr)
{
    copyFrom(other);
}

Graph& Graph::operator=(const Graph& other)
{
    if (this != &other) {
        free();
        copyFrom(other);
    }
    return *this;
}

void Graph::free()
{
    if (tail != nullptr) {
        delete[] tail;
        tail = nullptr;
    }
    
    if (head != nullptr) {
        delete[] head;
        head = nullptr;
    }
    
    if (adj_matrix != nullptr) {
        for (int i = 0; i < nnodes; i++) {
            delete[] adj_matrix[i];
        }
        delete[] adj_matrix;
        adj_matrix = nullptr;
    }
    
    if (adj_matrix_compl != nullptr) {
        for (int i = 0; i < nnodes; i++) {
            delete[] adj_matrix_compl[i];
        }
        delete[] adj_matrix_compl;
        adj_matrix_compl = nullptr;
    }
    
    if (adj_lists != nullptr) {
        for (int i = 0; i < nnodes; i++) {
            if (adj_lists[i] != nullptr) {
                delete[] adj_lists[i];
            }
        }
        delete[] adj_lists;
        adj_lists = nullptr;
    }
    
    if (node_degree != nullptr) {
        delete[] node_degree;
        node_degree = nullptr;
    }

    if (vertex_weights != nullptr) {
        delete[] vertex_weights;
        vertex_weights = nullptr;
    }
    
    if (cliques != nullptr) {
        delete cliques;
        cliques = nullptr;
    }
}

void Graph::copyFrom(const Graph& other)
{
    nnodes = other.nnodes;
    nedges = other.nedges;
    
    // Copy edge arrays
    if (other.tail != nullptr) {
        tail = new int[nedges];
        memcpy(tail, other.tail, nedges * sizeof(int));
    }
    
    if (other.head != nullptr) {
        head = new int[nedges];
        memcpy(head, other.head, nedges * sizeof(int));
    }
    
    // Copy adjacency matrix
    if (other.adj_matrix != nullptr) {
        adj_matrix = new bool*[nnodes];
        for (int i = 0; i < nnodes; i++) {
            adj_matrix[i] = new bool[nnodes];
            memcpy(adj_matrix[i], other.adj_matrix[i], nnodes * sizeof(bool));
        }
    }
    
    // Copy complementary adjacency matrix
    if (other.adj_matrix_compl != nullptr) {
        adj_matrix_compl = new bool*[nnodes];
        for (int i = 0; i < nnodes; i++) {
            adj_matrix_compl[i] = new bool[nnodes];
            memcpy(adj_matrix_compl[i], other.adj_matrix_compl[i], nnodes * sizeof(bool));
        }
    }
    
    // Copy node degrees
    if (other.node_degree != nullptr) {
        node_degree = new int[nnodes];
        memcpy(node_degree, other.node_degree, nnodes * sizeof(int));
    }

    // Copy vertex weights
    if (other.vertex_weights != nullptr) {
        vertex_weights = new int[nnodes];
        memcpy(vertex_weights, other.vertex_weights, nnodes * sizeof(int));
    }
    
    // Copy adjacency lists
    if (other.adj_lists != nullptr && other.node_degree != nullptr) {
        adj_lists = new int*[nnodes];
        for (int i = 0; i < nnodes; i++) {
            if (node_degree[i] > 0) {
                adj_lists[i] = new int[node_degree[i]];
                memcpy(adj_lists[i], other.adj_lists[i], node_degree[i] * sizeof(int));
            } else {
                adj_lists[i] = nullptr;
            }
        }
    }
    
    // Copy cliques if computed
    if (other.cliques != nullptr) {
        cliques = new vector<vector<bool>>(*other.cliques);
    }
}

void Graph::buildAdjacencyMatrix()
{
    // cout << "Building adjacency matrix for " << nnodes << " nodes from " << nedges << " edges..." << endl;
    
    // Allocate adjacency matrix
    adj_matrix = new bool*[nnodes];
    for (int i = 0; i < nnodes; i++) {
        adj_matrix[i] = new bool[nnodes];
        memset(adj_matrix[i], 0, nnodes * sizeof(bool));
    }
    
    // Allocate complementary adjacency matrix (initialized to true)
    adj_matrix_compl = new bool*[nnodes];
    for (int i = 0; i < nnodes; i++) {
        adj_matrix_compl[i] = new bool[nnodes];
        for (int j = 0; j < nnodes; j++) {
            adj_matrix_compl[i][j] = (i != j);  
        }
    }
    
    // Fill adjacency matrix from edge list
    for (int e = 0; e < nedges; e++) {
        int u = tail[e];
        int v = head[e];
        
        adj_matrix[u][v] = true;
        adj_matrix[v][u] = true;
        
        // Update complementary matrix
        adj_matrix_compl[u][v] = false;
        adj_matrix_compl[v][u] = false;
    }
    
    // cout << "Adjacency matrix built successfully." << endl;
}

void Graph::buildAdjacencyLists()
{
    // cout << "Building adjacency lists..." << endl;
    
    // Initialize node degrees
    node_degree = new int[nnodes];
    memset(node_degree, 0, nnodes * sizeof(int));
    
    // Count degrees
    for (int e = 0; e < nedges; e++) {
        int u = tail[e];
        int v = head[e];
        
        if (u >= 0 && u < nnodes && v >= 0 && v < nnodes) {
            node_degree[u]++;
            node_degree[v]++;
        }
    }
    
    // Print degree statistics
    int max_degree = 0;
    int min_degree = nnodes;
    int total_degree = 0;
    for (int i = 0; i < nnodes; i++) {
        total_degree += node_degree[i];
        if (node_degree[i] > max_degree) max_degree = node_degree[i];
        if (node_degree[i] < min_degree) min_degree = node_degree[i];
    }

    // cout << "Degree statistics: min=" << min_degree << ", max=" << max_degree
    //      << ", avg=" << (total_degree / (double)nnodes) << endl;

    // Allocate adjacency lists
    adj_lists = new int*[nnodes];
    vector<int> current_sizes(nnodes, 0);  
    for (int i = 0; i < nnodes; i++) {
        if (node_degree[i] > 0) {
            adj_lists[i] = new int[node_degree[i]];
        } else {
            adj_lists[i] = nullptr;
        }
    }
    
    // Fill adjacency lists
    for (int e = 0; e < nedges; e++) {
        int u = tail[e];
        int v = head[e];
        
        if (u >= 0 && u < nnodes && v >= 0 && v < nnodes) {
            adj_lists[u][current_sizes[u]++] = v;
            adj_lists[v][current_sizes[v]++] = u;
        }
    }
    
    // cout << "Adjacency lists built successfully." << endl;
}

// ============================================================================
// Graph Query Methods
// ============================================================================

bool Graph::hasEdge(int u, int v) const
{
    if (u < 0 || u >= nnodes || v < 0 || v >= nnodes) {
        return false;
    }
    return adj_matrix[u][v];
}

const int* Graph::getNeighbors(int node, int& num_neighbors) const
{
    if (node < 0 || node >= nnodes) {
        num_neighbors = 0;
        return nullptr;
    }
    num_neighbors = node_degree[node];
    return adj_lists[node];
}

int Graph::getDegree(int node) const
{
    if (node < 0 || node >= nnodes) {
        return 0;
    }
    return node_degree[node];
}

// ============================================================================
// Graph Algorithms
// ============================================================================

int countConnectedComponents(const Graph& graph)
{
    if (graph.nnodes == 0) {
        return 0;
    }
    
    vector<bool> visited(graph.nnodes, false);
    int num_components = 0;
    
    // DFS implementation using stack
    for (int start = 0; start < graph.nnodes; start++) {
        if (visited[start]) {
            continue;
        }
        
        // Found a new component
        num_components++;
        
        // DFS from start node
        stack<int> st;
        st.push(start);
        visited[start] = true;
        
        while (!st.empty()) {
            int node = st.top();
            st.pop();
            
            // Visit all neighbors
            int num_neighbors;
            const int* neighbors = graph.getNeighbors(node, num_neighbors);
            for (int i = 0; i < num_neighbors; i++) {
                int neighbor = neighbors[i];
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    st.push(neighbor);
                }
            }
        }
    }
    
    return num_components;
}

int countConnectedComponents(const Graph& graph, const bool* removed)
{
    if (graph.nnodes == 0) {
        return 0;
    }
    
    vector<bool> visited(graph.nnodes, false);
    int num_components = 0;
    
    // DFS implementation using stack (skipping removed nodes)
    for (int start = 0; start < graph.nnodes; start++) {
        // Skip removed nodes
        if (removed[start] || visited[start]) {
            continue;
        }
        
        // Found a new component
        num_components++;
        
        // DFS from start node
        stack<int> st;
        st.push(start);
        visited[start] = true;
        
        while (!st.empty()) {
            int node = st.top();
            st.pop();
            
            // Visit all neighbors (that are not removed)
            int num_neighbors;
            const int* neighbors = graph.getNeighbors(node, num_neighbors);
            for (int i = 0; i < num_neighbors; i++) {
                int neighbor = neighbors[i];
                if (!removed[neighbor] && !visited[neighbor]) {
                    visited[neighbor] = true;
                    st.push(neighbor);
                }
            }
        }
    }
    
    return num_components;
}

void findMaximalClique(const Graph& graph, int u, int v, vector<bool>& clique)
{
    clique.assign(graph.nnodes, false);
    clique[u] = true;
    clique[v] = true;
        
    // Check if u and v are actually adjacent
    if (!graph.hasEdge(u, v)) {
        cerr << "    ERROR: Nodes " << u << " and " << v << " are NOT adjacent!" << endl;
        return;
    }
    
    // Try to extend the clique: add vertices that are adjacent to all current clique members
    int extensions = 0;
    for (int w = 0; w < graph.nnodes; w++) {
        if (clique[w]) {
            continue;
        }
        
        // Check if w is adjacent to all vertices in current clique
        bool adjacent_to_all = true;
        for (int i = 0; i < graph.nnodes; i++) {
            if (clique[i] && !graph.hasEdge(w, i)) {
                adjacent_to_all = false;
                break;
            }
        }
        
        if (adjacent_to_all) {
            clique[w] = true;
            extensions++;
        }
    }
}

int buildEdgeCoveringCliques(const Graph& graph, vector<vector<bool>>& cliques, int clique_option)
{
    cliques.clear();

    if(clique_option == 0)
    {
        // clique_option==0 -> each edge is its own clique
        for (int e = 0; e < graph.nedges; e++) {
            int u = graph.tail[e];
            int v = graph.head[e];
            vector<bool> clique(graph.nnodes, false);
            clique[u] = true;
            clique[v] = true;
            cliques.push_back(clique);
        }
        return (int)cliques.size();
    }

    
    if (graph.nedges == 0) {
        cout << "WARNING: No edges to cover!" << endl;
        return 0;
    }
    
    // Track which edges are covered
    vector<bool> edge_covered(graph.nedges, false);
    int covered_edges = 0;
    
    // Greedy algorithm: for each uncovered edge, find a maximal clique containing it
    int clique_count = 0;
    while (covered_edges < graph.nedges) {
        // Find first uncovered edge
        int edge_idx = -1;
        for (int e = 0; e < graph.nedges; e++) {
            if (!edge_covered[e]) {
                edge_idx = e;
                break;
            }
        }
        
        if (edge_idx == -1) {
            cout << "WARNING: No uncovered edge found but covered_edges=" << covered_edges 
                 << " < nedges=" << graph.nedges << endl;
            break;
        }
        
        int u = graph.tail[edge_idx];
        int v = graph.head[edge_idx];
        
      
        // Validation check
        if (u < 0 || u >= graph.nnodes || v < 0 || v >= graph.nnodes) {
            cerr << "ERROR: Invalid edge indices! u=" << u << ", v=" << v 
                 << " (valid range: 0-" << (graph.nnodes-1) << ")" << endl;
            edge_covered[edge_idx] = true;
            covered_edges++;
            continue;
        }
        
        // Find maximal clique containing edge (u,v)
        vector<bool> clique;
        
        if (clique_option == 1) {
            // Edge partition mode: only add vertices whose edges to the clique are all uncovered
            clique.assign(graph.nnodes, false);
            clique[u] = true;
            clique[v] = true;
            
            // Try to extend the clique with vertices that:
            // 1. Are adjacent to all current clique members
            // 2. Have all edges to clique members uncovered
            for (int w = 0; w < graph.nnodes; w++) {
                if (clique[w]) continue;
                
                // Check if w is adjacent to all vertices in current clique
                bool adjacent_to_all = true;
                bool all_edges_uncovered = true;
                
                for (int i = 0; i < graph.nnodes; i++) {
                    if (clique[i]) {
                        if (!graph.hasEdge(w, i)) {
                            adjacent_to_all = false;
                            break;
                        }
                        // Check if edge (w, i) is already covered
                        for (int e = 0; e < graph.nedges; e++) {
                            int t = graph.tail[e];
                            int h = graph.head[e];
                            if ((t == w && h == i) || (t == i && h == w)) {
                                if (edge_covered[e]) {
                                    all_edges_uncovered = false;
                                }
                                break;
                            }
                        }
                        if (!all_edges_uncovered) break;
                    }
                }
                
                if (adjacent_to_all && all_edges_uncovered) {
                    clique[w] = true;
                }
            }
        } else {
            // clique_option == 2: Edge cover mode with maximal cliques (edges can be covered multiple times)
            findMaximalClique(graph, u, v, clique);
        }
        
        // Count clique size
        int clique_size = 0;
        for (int i = 0; i < graph.nnodes; i++) {
            if (clique[i]) clique_size++;
        }

        cliques.push_back(clique);
        clique_count++;
        
        // Mark all edges covered by this clique
        for (int e = 0; e < graph.nedges; e++) {
            if (!edge_covered[e]) {
                int t = graph.tail[e];
                int h = graph.head[e];
                
                if (t >= 0 && t < graph.nnodes && h >= 0 && h < graph.nnodes && 
                    clique[t] && clique[h]) {
                    edge_covered[e] = true;
                    covered_edges++;
                }
            }
        }
    }
    
    return (int)cliques.size();
}

Graph* buildResidualGraph(const Graph& graph, const bool* removed,
                          vector<int>& old_to_new, vector<int>& new_to_old)
{
    // cout << "\n=== Building Residual Graph ===" << endl;
    
    // Count number of remaining nodes
    int num_remaining = 0;
    for (int i = 0; i < graph.nnodes; i++) {
        if (!removed[i]) {
            num_remaining++;
        }
    }
    
    // cout << "Original nodes: " << graph.nnodes << ", Removed: " << (graph.nnodes - num_remaining) 
    //      << ", Remaining: " << num_remaining << endl;
    
    // Build node mapping
    old_to_new.assign(graph.nnodes, -1);
    new_to_old.clear();
    new_to_old.reserve(num_remaining);
    
    int new_idx = 0;
    for (int old_idx = 0; old_idx < graph.nnodes; old_idx++) {
        if (!removed[old_idx]) {
            old_to_new[old_idx] = new_idx;
            new_to_old.push_back(old_idx);                
            new_idx++;
        }
    }

    // Map vertex weights for residual graph
    int* vertex_weights_residual = new int[num_remaining];
    for (int i = 0; i < num_remaining; i++) {
        int old_idx = new_to_old[i];
        vertex_weights_residual[i] = graph.vertex_weights[old_idx];
    }
    
    // Count edges in residual graph
    int num_edges_residual = 0;
    for (int e = 0; e < graph.nedges; e++) {
        int u = graph.tail[e];
        int v = graph.head[e];
        
        // Edge remains only if both endpoints are not removed
        if (!removed[u] && !removed[v]) {
            num_edges_residual++;
        }
    }
    
    // cout << "Original edges: " << graph.nedges << ", Residual edges: " << num_edges_residual << endl;
    
    if (num_edges_residual == 0) {
        cout << "WARNING: Residual graph has no edges!" << endl;
        exit(-1);
    }
    
    // Build edge arrays for residual graph
    int* tail_residual = new int[num_edges_residual];
    int* head_residual = new int[num_edges_residual];
    
    int edge_count = 0;
    for (int e = 0; e < graph.nedges; e++) {
        int u = graph.tail[e];
        int v = graph.head[e];
        
        if (!removed[u] && !removed[v]) {
            // Map old indices to new indices
            tail_residual[edge_count] = old_to_new[u];
            head_residual[edge_count] = old_to_new[v];
            edge_count++;
        }
    }
    
    // Create residual graph (it will copy the arrays internally)
    Graph* residual = new Graph(num_remaining, num_edges_residual, 
                                 tail_residual, head_residual, vertex_weights_residual, true);
    
    // Clean up temporary arrays
    delete[] tail_residual;
    delete[] head_residual;
    delete[] vertex_weights_residual;
    
    // cout << "Residual graph created successfully." << endl;
    // cout << "================================\n" << endl;
    
    return residual;
}



int computeVertexConnectivity(const Graph& graph, bool per_component, bool* separator)
{
    // cout << "\n=== Computing Vertex Connectivity ===" << endl;
    
    int n = graph.nnodes;
    
    // Initialize separator if provided (binary array of length nnodes)
    if (separator != nullptr) {
        memset(separator, 0, n * sizeof(bool));
    }
    
    // Special cases
    if (n <= 1) {
        // cout << "Graph has <= 1 node, connectivity = 0" << endl;
        return 0;
    }
    
    // Check if graph is connected
    int num_components = countConnectedComponents(graph);
    if (num_components > 1) {
        // cout << "Graph is disconnected (" << num_components << " components)" << endl;
        if (!per_component) {
            // cout << "Returning connectivity = 0 (graph is already disconnected)" << endl;
            return 0;
        }
        
        // Compute connectivity for each component and return minimum
        // cout << "Computing minimum connectivity among all components..." << endl;
        
        // Identify components using DFS
        vector<int> component_id(n, -1);
        vector<bool> visited(n, false);
        int comp_count = 0;
        
        for (int start = 0; start < n; start++) {
            if (visited[start]) continue;
            
            // DFS to mark component
            stack<int> st;
            st.push(start);
            visited[start] = true;
            component_id[start] = comp_count;
            
            while (!st.empty()) {
                int node = st.top();
                st.pop();
                
                int num_neighbors;
                const int* neighbors = graph.getNeighbors(node, num_neighbors);
                for (int i = 0; i < num_neighbors; i++) {
                    int neighbor = neighbors[i];
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        component_id[neighbor] = comp_count;
                        st.push(neighbor);
                    }
                }
            }
            comp_count++;
        }
        
        // Build residual graph for each component and compute its connectivity
        int min_comp_connectivity = n; // Initialize to large value
        bool* best_comp_separator = nullptr;
        vector<int> best_old_to_new, best_new_to_old;
        
        for (int comp = 0; comp < comp_count; comp++) {
            // Mark nodes not in this component as removed
            bool* removed = new bool[n];
            int comp_size = 0;
            for (int i = 0; i < n; i++) {
                removed[i] = (component_id[i] != comp);
                if (!removed[i]) comp_size++;
            }

            // cout << "\n  Component " << comp << " has " << comp_size << " nodes" << endl;

            if (comp_size <= 2) {
                // cout << "    Component too small, connectivity = " << (comp_size - 1) << endl;
                delete[] removed;
                continue;
            }
            
            // Build residual graph for this component
            vector<int> old_to_new, new_to_old;
            Graph* comp_graph = buildResidualGraph(graph, removed, old_to_new, new_to_old);
            
            // Allocate separator for this component (in component coordinates)
            bool* comp_separator = new bool[comp_graph->nnodes];
            
            // Recursively compute connectivity of this component with separator
            int comp_connectivity = computeVertexConnectivity(*comp_graph, false, comp_separator);
            // cout << "    Component " << comp << " connectivity = " << comp_connectivity << endl;
            
            if (comp_connectivity < min_comp_connectivity) {
                min_comp_connectivity = comp_connectivity;
                
                // Save this separator and mappings for later conversion
                if (best_comp_separator != nullptr) {
                    delete[] best_comp_separator;
                }
                best_comp_separator = comp_separator;
                best_old_to_new = old_to_new;
                best_new_to_old = new_to_old;
            } else {
                delete[] comp_separator;
            }
            
            delete comp_graph;
            delete[] removed;
        }
        
        // Map the best separator back to original graph coordinates
        if (separator != nullptr && best_comp_separator != nullptr) {
            memset(separator, 0, n * sizeof(bool));
            for (int i = 0; i < (int)best_new_to_old.size(); i++) {
                if (best_comp_separator[i]) {
                    int original_idx = best_new_to_old[i];
                    separator[original_idx] = true;
                }
            }
            
            // Print the separator in original coordinates
            // cout << "\nMinimum separator (in original graph coordinates): {";
            // bool first = true;
            // for (int v = 0; v < n; v++) {
            //     if (separator[v]) {
            //         if (!first) cout << ", ";
            //         cout << v;
            //         first = false;
            //     }
            // }
            // cout << "}" << endl;

        }
        
        // Clean up
        if (best_comp_separator != nullptr) {
            delete[] best_comp_separator;
        }
        
        // cout << "\nMinimum connectivity among (non singleton) components = " << min_comp_connectivity << endl;
        // cout << "=====================================\n" << endl;

        if (min_comp_connectivity == n) {
            // All components were singletons or too small
            return 0;
        }
        return min_comp_connectivity;
    }
    
    // cout << "Graph is connected" << endl;
    
    // Check if graph is complete
    bool is_complete = true;
    for (int i = 0; i < n && is_complete; i++) {
        if (graph.node_degree[i] != n - 1) {
            is_complete = false;
            break;
        }
    }
    
    if (is_complete) {
        // For complete graph, connectivity is sum of weights of all nodes except one
        int total_weight = 0;
        for (int i = 0; i < n; i++) {
            total_weight += graph.vertex_weights[i];
        }
        // Find minimum node weight
        int min_weight = graph.vertex_weights[0];
        for (int i = 1; i < n; i++) {
            if (graph.vertex_weights[i] < min_weight) {
                min_weight = graph.vertex_weights[i];
            }
        }
        // cout << "Graph is complete K" << n << ", weighted connectivity = " << (total_weight - min_weight) << endl;
        return total_weight - min_weight;
    }
    
    // Find node with minimum weighted degree (upper bound on connectivity)
    int min_weighted_degree = INT_MAX;
    int min_degree_node = -1;
    for (int i = 0; i < n; i++) {
        // Compute sum of weights of neighbors
        int weighted_degree = 0;
        int num_neighbors;
        const int* neighbors = graph.getNeighbors(i, num_neighbors);
        for (int j = 0; j < num_neighbors; j++) {
            weighted_degree += graph.vertex_weights[neighbors[j]];
        }
        
        if (weighted_degree < min_weighted_degree) {
            min_weighted_degree = weighted_degree;
            min_degree_node = i;
        }
    }
    
    // cout << "Graph has " << n << " nodes, " << graph.nedges << " edges" << endl;
    // cout << "Minimum weighted degree: " << min_weighted_degree << " (upper bound on weighted connectivity)" << endl;
    
    // Initialize minimum connectivity to min_weighted_degree
    int min_connectivity = min_weighted_degree;
    
    // Keep track of best separator found (binary array)
    bool* best_separator = new bool[n];
    memset(best_separator, 0, n * sizeof(bool));
    int best_source = -1;
    // int best_sink = -1;
    
    // Initialize best_separator with neighbors of minimum weighted degree node
    // (this is a valid separator that disconnects that node from the rest)
    if (min_degree_node != -1) {
        int num_neighbors;
        const int* neighbors = graph.getNeighbors(min_degree_node, num_neighbors);
        for (int i = 0; i < num_neighbors; i++) {
            best_separator[neighbors[i]] = true;
        }
        best_source = min_degree_node;
        // cout << "Initial separator: neighbors of node " << min_degree_node 
        //      << " (weighted size " << min_weighted_degree << ")" << endl;
    }
    
    // We need to check max-flow between non-adjacent pairs
   
    // For each pair of non-adjacent nodes, compute max-flow with node capacities
    for (int source = 0; source < n && min_connectivity > 1; source++) {
        for (int sink = source + 1; sink < n && min_connectivity > 1; sink++) {
            // Skip adjacent pairs
            if (graph.hasEdge(source, sink)) {
                continue;
            }
            
            // Build auxiliary digraph with node splitting
            // Each node v becomes v_in and v_out with arc (v_in, v_out) of capacity w_v
            lemon::ListDigraph digraph;
            lemon::ListDigraph::ArcMap<int> capacity(digraph);
            
            // Create nodes: for each original node, create in-node and out-node
            vector<lemon::ListDigraph::Node> node_in(n);
            vector<lemon::ListDigraph::Node> node_out(n);
            
            for (int v = 0; v < n; v++) {
                node_in[v] = digraph.addNode();
                node_out[v] = digraph.addNode();
                
                // Connect in to out with capacity equal to vertex weight (except for source and sink)
                if (v != source && v != sink) {
                    lemon::ListDigraph::Arc arc = digraph.addArc(node_in[v], node_out[v]);
                    capacity[arc] = graph.vertex_weights[v];
                } else {
                    // Source and sink have infinite internal capacity
                    lemon::ListDigraph::Arc arc = digraph.addArc(node_in[v], node_out[v]);
                    // Use sum of all weights as "infinite" capacity
                    int total_weight = 0;
                    for (int i = 0; i < n; i++) {
                        total_weight += graph.vertex_weights[i];
                    }
                    capacity[arc] = total_weight;
                }
            }
            
            // Add edges: for each edge (u,v) in original graph, add arcs in both directions
            // Compute total weight for "infinite" capacity on edges
            int total_weight = 0;
            for (int i = 0; i < n; i++) {
                total_weight += graph.vertex_weights[i];
            }
            
            for (int e = 0; e < graph.nedges; e++) {
                int u = graph.tail[e];
                int v = graph.head[e];
                
                // Arc from u_out to v_in with infinite capacity
                lemon::ListDigraph::Arc arc1 = digraph.addArc(node_out[u], node_in[v]);
                capacity[arc1] = total_weight; // large enough
                
                // Arc from v_out to u_in with infinite capacity
                lemon::ListDigraph::Arc arc2 = digraph.addArc(node_out[v], node_in[u]);
                capacity[arc2] = total_weight; // large enough
            }
            
            // Compute max flow from source_out to sink_in
            lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<int>> 
                preflow(digraph, capacity, node_out[source], node_in[sink]);
            
            preflow.runMinCut();
            int flow_value = preflow.flowValue();
            
            // flow_value is already the total weight of the separator
            if (flow_value < min_connectivity) {
                min_connectivity = flow_value;
                best_source = source;
                // best_sink = sink;
                
                // Extract minimum vertex separator from min-cut
                // A node v is in the separator if the arc (v_in, v_out) is in the cut
                memset(best_separator, 0, n * sizeof(bool));
                for (int v = 0; v < n; v++) {
                    // Skip source and sink (they can't be in the separator for this pair)
                    if (v == source || v == sink) continue;
                    
                    // Check if v_in is on source side and v_out is on sink side
                    // This means the arc (v_in, v_out) is in the cut
                    if (preflow.minCut(node_in[v]) && !preflow.minCut(node_out[v])) {
                        best_separator[v] = true;
                    }
                }
                
            }
        }
    }
    
    // cout << "Weighted vertex connectivity = " << min_connectivity << endl;
    
    // Store the separator if requested
    if (separator != nullptr && best_source != -1) {
        memcpy(separator, best_separator, n * sizeof(bool));
        // Compute separator weight and size for debugging
        // int sep_count = 0;
        // int sep_weight = 0;
        // for (int v = 0; v < n; v++) {
        //     if (best_separator[v]) {
        //         sep_count++;
        //         sep_weight += graph.vertex_weights[v];
        //     }
        // }
        // cout << "Minimum vertex separator:" << endl;
        // cout << "Separates node " << best_source;
        // if (best_sink != -1) {
        //     cout << " and node " << best_sink;
        // } else {
        //     cout << " from the rest of the graph";
        // }
        // cout << endl;
        // cout << "Separator size: " << sep_count << ", weight: " << sep_weight << endl;
    }
    
    // Clean up
    delete[] best_separator;

    // cout << "=====================================\n" << endl;

    return min_connectivity;
}

// ============================================================================
// Optional Computed Properties
// ============================================================================

int Graph::computeEdgeCoveringCliques(int clique_option)
{
    // cout << "\n=== Computing Edge-Covering Cliques ===" << endl;
    
    // Free previous cliques if any
    if (cliques != nullptr) {
        delete cliques;
    }
    
    // Allocate new vector
    cliques = new vector<vector<bool>>();
    
    // Use the standalone function to compute cliques
    int n_cliques = buildEdgeCoveringCliques(*this, *cliques, clique_option);
    
    // cout << "Computed " << n_cliques << " cliques covering all edges" << endl;
    // cout << "========================================\n" << endl;
    
    return n_cliques;
}

const vector<vector<bool>>* Graph::getCliques() const
{
    return cliques;
}

int Graph::getNCliques() const
{
    return (cliques != nullptr) ? cliques->size() : 0;
}

// ============================================================================
// Isolated Cliques Detection
// ============================================================================

int identifyIsolatedCliques(const Graph& graph, const bool* preFixed, bool* is_in_isolated_clique)
{
    int n = graph.nnodes;
    
    // Initialize output array
    memset(is_in_isolated_clique, 0, n * sizeof(bool));
    
    // Track visited nodes for DFS
    vector<bool> visited(n, false);
    vector<int> component_nodes; 
    
    int total_nodes_in_cliques = 0;
    int num_cliques = 0;
    
    // Explore all connected components
    for (int start = 0; start < n; start++) {
        // Skip already visited or pre-fixed nodes
        if (visited[start] || (preFixed != nullptr && preFixed[start])) {
            continue;
        }
        
        // Find all nodes in this component using DFS
        component_nodes.clear();
        stack<int> st;
        st.push(start);
        visited[start] = true;
        
        while (!st.empty()) {
            int node = st.top();
            st.pop();
            component_nodes.push_back(node);
            
            // Visit all neighbors
            int num_neighbors;
            const int* neighbors = graph.getNeighbors(node, num_neighbors);
            for (int i = 0; i < num_neighbors; i++) {
                int neighbor = neighbors[i];
                if (!visited[neighbor] && (preFixed == nullptr || !preFixed[neighbor])) {
                    visited[neighbor] = true;
                    st.push(neighbor);
                }
            }
        }
        
        // Check if this component is a clique
        int comp_size = component_nodes.size();
        bool is_clique = true;
        
        // Singleton is always a clique!
        if (comp_size == 1) {
            is_clique = true;
        }
        else {
            // Check all pairs of nodes in the component
            for (int i = 0; i < comp_size && is_clique; i++) {
                for (int j = i + 1; j < comp_size && is_clique; j++) {
                    if (!graph.hasEdge(component_nodes[i], component_nodes[j])) {
                        is_clique = false;
                    }
                }
            }
        }
        
        // If it's a clique, mark all its nodes
        if (is_clique) {
            num_cliques++;
            for (int node : component_nodes) {
                is_in_isolated_clique[node] = true;
                total_nodes_in_cliques++;
            }
            
            // cout << "  Found isolated clique of size " << comp_size << ": {";
            // for (size_t i = 0; i < component_nodes.size(); i++) {
            //     cout << component_nodes[i];
            //     if (i < component_nodes.size() - 1) cout << ", ";
            // }
            // cout << "}" << endl;
        }
    }
    
    // cout << "Total isolated cliques found: " << num_cliques << endl;
    // cout << "Total nodes in isolated cliques: " << total_nodes_in_cliques << endl;
    // cout << "====================================\n" << endl;
    
    return num_cliques;
}



