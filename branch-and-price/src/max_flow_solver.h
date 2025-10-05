#ifndef MAX_FLOW_SOLVER_H
#define MAX_FLOW_SOLVER_H

#include <vector>
#include <queue>
#include <deque>
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;

struct PushRelabel {
  using T = double;

  struct Edge {
    int to, rev;
    T   res;
    T   cap0;
  };

  // Graph structure
  int n;
  vector<vector<Edge>> g;
  vector<int> len_adj;
  vector<vector<int>> adj_ind;

  // Runtime state
  int s = -1, t = -1;
  vector<T> excess;
  vector<int> height;
  vector<int> cur;
  vector<int> cnt;

  // Buckets per height
  vector<int> headBucket;
  vector<int> nextActive;
  vector<char> active;
  int H = 0;

  // Global relabel trigger
  long long work_counter = 0;
  long long GR_THRESHOLD = 0;

  // Numerical tolerance
  T epsRel = 1e-12;
  T epsAbs = 1e-12;

  // METHODS
  explicit PushRelabel(int n_);
  
  void set_relative_epsilon(T eps);
  void add_edge(int u, int v, T c);
  void update_edge_capacity(int u, int v, T new_c);
  PushRelabel create_pricing_copy() const;
  
  T max_flow(int s_, int t_);
  void mincut_source_mask(vector<unsigned char> &Smask) const;
  T mincut_value_from_mask(const vector<unsigned char> &Smask) const;

private:
  inline void buckets_clear();
  inline void activate(int v);
  inline int pop_highest();
  inline void clamp_zero(T &x);
  inline void push_edge(int u, Edge &e);
  inline void relabel_node(int u);
  inline void gap_heuristic(int k);
  void global_relabel();
  void discharge(int u);
};

#endif 