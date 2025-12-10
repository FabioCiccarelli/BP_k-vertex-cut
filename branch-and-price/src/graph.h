/**@file   graph.h
 * @brief  Graph data structure and related functions for k-vertex cut problem
 * @author Fabio Ciccarelli
 *
 * This file defines a graph structure that contains all graph-related data
 * and operations. It is designed to be independent from SCIP-specific types
 * to allow for better future reusability.
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

/**
 * @brief Graph structure for undirected graphs
 * 
 * This structure contains all the essential data for representing an undirected graph:
 * - Number of nodes and edges
 * - Adjacency matrix 
 * - Adjacency lists for fast neighbor iteration
 * - Node degrees
 * - Edge representation (tail and head arrays)
 * - Optional computed properties (cliques)
 */
struct Graph {
    int nnodes;                        ///< Number of nodes in the graph
    int nedges;                        ///< Number of edges in the graph
    
    int* tail;                         ///< Array of edge tail nodes (0-indexed)
    int* head;                         ///< Array of edge head nodes (0-indexed)
    
    bool** adj_matrix;                 ///< Adjacency matrix [i][j] = true if edge (i,j) exists
    bool** adj_matrix_compl;           ///< Complementary adjacency matrix [i][j] = true if edge (i,j) does NOT exist
    int** adj_lists;                   ///< Adjacency lists: adj_lists[i] contains neighbors of node i
    int* node_degree;                  ///< Degree of each node (also size of adj_lists[i])
    int* vertex_weights;               ///< Weights of each vertex

    /* Optional computed properties */
    std::vector<std::vector<bool>>* cliques;  ///< Edge-covering family of cliques (NULL if not computed)
    
    /**
     * @brief Constructor - creates an empty graph
     */
    Graph();
    
    /**
     * @brief Constructor - creates a graph from edge list
     * @param n_nodes number of nodes
     * @param n_edges number of edges
     * @param edge_tails array of edge tail nodes (0-indexed)
     * @param edge_heads array of edge head nodes (0-indexed)
     * @param vertex_weights array of vertex weights
     * @param copy_edges if true, copies the edge arrays; if false, takes ownership
     */
    Graph(int n_nodes, int n_edges, int* edge_tails, int* edge_heads, int* vertex_weights, bool copy_edges = true);
    
    /**
     * @brief Destructor - frees all allocated memory
     */
    ~Graph();
    
    /**
     * @brief Copy constructor
     */
    Graph(const Graph& other);
    
    /**
     * @brief Assignment operator
     */
    Graph& operator=(const Graph& other);
    
    /**
     * @brief Checks if an edge exists between two nodes
     * @param u first node (0-indexed)
     * @param v second node (0-indexed)
     * @return true if edge exists, false otherwise
     */
    bool hasEdge(int u, int v) const;
    
    /**
     * @brief Gets the neighbors of a node
     * @param node node index (0-indexed)
     * @param num_neighbors output parameter for number of neighbors
     * @return pointer to array of neighbor indices
     */
    const int* getNeighbors(int node, int& num_neighbors) const;
    
    /**
     * @brief Gets the degree of a node
     * @param node node index (0-indexed)
     * @return degree of the node
     */
    int getDegree(int node) const;
    
    /**
     * @brief Computes and stores edge-covering family of cliques
     * 
     * Finds a family of cliques such that every edge is covered by at least one clique.
     * Results are stored internally and can be accessed via getCliques().
     * 
     * @return number of cliques computed
     */
    int computeEdgeCoveringCliques(int clique_option);
    
    /**
     * @brief Returns the computed cliques
     * @return pointer to cliques vector (NULL if not computed)
     */
    const std::vector<std::vector<bool>>* getCliques() const;
    
    /**
     * @brief Returns number of computed cliques
     * @return number of cliques (0 if not computed)
     */
    int getNCliques() const;
    
private:
    /**
     * @brief Builds the adjacency matrix from edge list
     */
    void buildAdjacencyMatrix();
    
    /**
     * @brief Builds adjacency lists and computes node degrees
     */
    void buildAdjacencyLists();
    
    /**
     * @brief Frees all allocated memory
     */
    void free();
    
    /**
     * @brief Deep copy from another graph
     */
    void copyFrom(const Graph& other);
};

/**
 * @brief Counts the number of connected components in a graph
 * 
 * Uses Depth-First Search (DFS) to find all connected components.
 * 
 * @param graph the graph to analyze
 * @return number of connected components
 */
int countConnectedComponents(const Graph& graph);

/**
 * @brief Counts the number of connected components in a graph with some nodes removed
 * 
 * This function considers a "residual graph" where some nodes are marked as removed.
 * It counts the connected components in the graph induced by the non-removed nodes.
 * 
 * @param graph the original graph
 * @param removed array indicating which nodes are removed (true = removed)
 * @return number of connected components in the residual graph
 */
int countConnectedComponents(const Graph& graph, const bool* removed);

/**
 * @brief Builds an edge-covering family of cliques using a greedy approach
 * 
 * Given a graph, finds a family of cliques such that every edge is covered
 * by at least one clique. Uses a greedy algorithm: repeatedly finds maximal
 * cliques until all edges are covered.
 * 
 * @param graph the graph to analyze
 * @param cliques output vector of cliques (each clique is a bool vector)
 * @return number of cliques in the covering
 */
int buildEdgeCoveringCliques(const Graph& graph, std::vector<std::vector<bool>>& cliques, bool only_edges);

/**
 * @brief Finds a maximal clique containing a given edge
 * 
 * Starting from an edge (u,v), extends it to a maximal clique by greedily
 * adding vertices that are adjacent to all current clique members.
 * 
 * @param graph the graph to analyze
 * @param u first endpoint of the edge (0-indexed)
 * @param v second endpoint of the edge (0-indexed)
 * @param clique output vector indicating which nodes are in the clique
 */
void findMaximalClique(const Graph& graph, int u, int v, std::vector<bool>& clique);

/**
 * @brief Builds a residual graph by removing specified nodes and their incident edges
 * 
 * Creates a new graph that is induced by the non-removed nodes. The returned graph
 * has nodes re-indexed from 0 to (num_remaining_nodes - 1). A mapping from new indices
 * to original indices is provided.
 * 
 * 
 * @param graph the original graph
 * @param removed array indicating which nodes to remove (true = remove)
 * @param old_to_new output: mapping from original node indices to new indices (-1 if removed)
 * @param new_to_old output: mapping from new node indices to original indices
 * @return pointer to newly allocated residual Graph (caller is responsible for deletion)
 */
Graph* buildResidualGraph(const Graph& graph, const bool* removed, 
                          std::vector<int>& old_to_new, std::vector<int>& new_to_old);

/**
 * @brief Computes the vertex connectivity of a graph
 * 
 * The vertex connectivity (or node connectivity) is the minimum weight subset of nodes
 * that must be removed to disconnect the graph or reduce it to a single node.
 * Uses Menger's theorem and max-flow computations with LEMON library.
 *
 * Algorithm:
 * - If graph is disconnected, returns 0 (or min connectivity of components if per_component=true)
 * - If graph is complete, returns the sum of the vertex weights minus the minimum vertex weight
 * - Otherwise, computes max-flow between non-adjacent pairs using node splitting
 *   (each node v is split into v_in and v_out with capacity w_v)
 * 
 * Time complexity: O(n^2 * max_flow_time) where max_flow_time depends on the algorithm
 * 
 * @param graph the graph to analyze
 * @param per_component if true and graph is disconnected, compute minimum connectivity
 *                      among all connected components; if false, return 0 for disconnected graphs
 * @param separator output: bool array of length nnodes where separator[v] = true if
 *                  vertex v is in the minimum vertex separator (if NULL, separator is not computed)
 *                  Note: caller must allocate array of size graph.nnodes
 * @return the vertex connectivity of the graph
 */
int computeVertexConnectivity(const Graph& graph, bool per_component = false, 
                               bool* separator = nullptr);

/**
 * @brief Identifies all isolated cliques in the graph (including singletons)
 * 
 * An isolated clique is a connected component where all nodes are pairwise adjacent.
 * This includes singleton nodes (which are cliques of size 1).
 * 
 * The function performs the following steps:
 * 1. Uses DFS to identify all connected components (excluding preFixed nodes)
 * 2. For each component, checks if it forms a clique
 * 3. Marks all nodes in isolated cliques in the output array
 * 
 * A component is a clique if:
 * - It has size 1 (singleton - always a clique)
 * - All pairs of nodes in the component are adjacent
 * 
 * @param graph the graph to analyze
 * @param preFixed array indicating which nodes are fixed by preprocessing (can be NULL)
 * @param is_in_isolated_clique output array of size graph.nnodes; 
 *                              is_in_isolated_clique[v] = true if node v is part of an isolated clique
 *                              (caller must allocate this array)
 * @return total number of nodes that are in isolated cliques
 */
int identifyIsolatedCliques(const Graph& graph, const bool* preFixed, bool* is_in_isolated_clique);


#endif // GRAPH_H
