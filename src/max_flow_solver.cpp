#include "max_flow_solver.h"

/*
  Push–Relabel (highest-label) with Gap Heuristic and Global Relabeling
  - Floating capacities (double) with a numerical tolerance (EPS)
  - Residual-only representation
  - Buckets implemented as singly-linked lists (head/next) per height level
  - Extracts source-side min-cut as a 0/1 mask

  Notes:
    - epsRel controls the relative tolerance; an absolute epsAbs is computed
      per run as epsAbs = epsRel * (1 + max_cap0) to scale with instance.
    - All "positivity" checks use (> epsAbs); small residuals/excess are clamped to zero.
*/



PushRelabel::PushRelabel(int n_) : n(n_), 
                                g(n_), 
                                len_adj(n_), 
                                adj_ind(n_, vector<int>(n_, -1)) {}

void PushRelabel::set_relative_epsilon(T eps) { epsRel = eps; }

void PushRelabel::add_edge(int u, int v, T c) {
    Edge a{v, (int)g[v].size(), c, c};
    Edge b{u, (int)g[u].size(), 0, 0};
    g[u].push_back(a);
    g[v].push_back(b);

    adj_ind[u][v] = len_adj[u];
    len_adj[u]++;
    len_adj[v]++;
    }

void PushRelabel::update_edge_capacity(int u, int v, T new_c) {
    // Note: this does not preserve flow feasibility if called mid-run (to do before solve)
    int ind = adj_ind[u][v];
    if (ind == -1) {
        throw runtime_error("Trying to update a non-existing edge");
    }
    Edge &e = g[u][ind];
    e.cap0 = new_c;
    e.res = new_c; // reset residual to new capacity
}

// Create a temporary copy of the graph (for pricing)
PushRelabel PushRelabel::create_pricing_copy() const {
    PushRelabel copy(n);
    // copy the graph structure
    copy.g = g; 
    copy.len_adj = len_adj;
    copy.adj_ind = adj_ind;
    copy.epsRel = epsRel;
    return copy;
}

// ---- Bucket management (linked lists) ----
inline void PushRelabel::buckets_clear() {
    int L = 2 * n + 2;
    headBucket.assign(L, -1);
    nextActive.assign(n, -1);
    active.assign(n, 0);
    H = 0;
}

inline void PushRelabel::activate(int v) {
    if (v == s || v == t) return;
    if (!active[v] && excess[v] > epsAbs) {
        int hv = height[v];
        nextActive[v] = headBucket[hv];
        headBucket[hv] = v;
        active[v] = 1;
        if (hv > H) H = hv;
    }
}

inline int PushRelabel::pop_highest() {
    while (H >= 0 && (H >= (int)headBucket.size() || headBucket[H] == -1)) --H;
    if (H < 0) return -1;
    int v = headBucket[H];
    headBucket[H] = nextActive[v];
    nextActive[v] = -1;
    active[v] = 0; // lazy removal: caller re-activates if still has excess
    return v;
}

// ---- Core primitives ----
inline void PushRelabel::clamp_zero(T &x) {
    if (abs(x) <= epsAbs) x = 0;
}

inline void PushRelabel::push_edge(int u, Edge &e) {
    if (height[u] != height[e.to] + 1) return;     // admissibility
    if (e.res <= epsAbs || excess[u] <= epsAbs) return;
    T f = min(excess[u], e.res);
    e.res -= f;
    g[e.to][e.rev].res += f;
    excess[u]    -= f;
    excess[e.to] += f;

    // clamp tiny numerical noise
    clamp_zero(e.res);
    clamp_zero(g[e.to][e.rev].res);
    clamp_zero(excess[u]);
    clamp_zero(excess[e.to]);

    work_counter += 1;
    activate(e.to);
}

inline void PushRelabel::relabel_node(int u) {
    int old = height[u];
    int newh = 2 * n;  // local "infinity"
    for (const auto &e : g[u]) if (e.res > epsAbs) {
        newh = min(newh, height[e.to] + 1);
    }
    cnt[old]--;
    height[u] = newh;
    if (newh >= (int)cnt.size()) cnt.resize(newh + 1, 0);
    cnt[newh]++;
    cur[u] = 0; // reset current arc after relabel
    if (newh > H) H = newh;
    work_counter += 1;
}

inline void PushRelabel::gap_heuristic(int k) {
    // If cnt[k] == 0 and k < n, nodes with h > k cannot reach t.
    for (int v = 0; v < n; ++v) {
        if (v == s || v == t) continue;
        int hv = height[v];
        if (hv > k && hv < n + 1) {
        cnt[hv]--;
        height[v] = n + 1;
        if (n + 1 >= (int)cnt.size()) cnt.resize(n + 2, 0);
        cnt[n + 1]++;
        active[v] = 0;       // lazy removal from old bucket
        nextActive[v] = -1;
        cur[v] = 0;
        }
    }
    if (n + 1 > H) H = n + 1;
}

void PushRelabel::global_relabel() {
    // BFS from t on the reverse residual graph: traverse y <- x if residual(y->x) > 0
    vector<int> dist(n, n + 1);
    deque<int> dq;
    dist[t] = 0;
    dq.push_back(t);
    while (!dq.empty()) {
        int x = dq.front(); dq.pop_front();
        for (const auto &a : g[x]) {
        const Edge &rev = g[a.to][a.rev];
        if (rev.res > epsAbs && dist[a.to] == n + 1) {
            dist[a.to] = dist[x] + 1;
            dq.push_back(a.to);
        }
        }
    }

    // rebuild heights and counts
    int L = 2 * n + 2;
    cnt.assign(L, 0);
    for (int v = 0; v < n; ++v) {
        height[v] = dist[v];
        if (height[v] > 2 * n) height[v] = 2 * n; // guard
    }
    height[s] = n;
    height[t] = 0;
    for (int v = 0; v < n; ++v) ++cnt[height[v]];

    // rebuild buckets with current excess
    buckets_clear();
    for (int v = 0; v < n; ++v) {
        cur[v] = 0;
        if (v != s && v != t && excess[v] > epsAbs) activate(v);
    }
    work_counter = 0;
    }

    void PushRelabel::discharge(int u) {
    while (excess[u] > epsAbs) {
        // scan admissible edges starting from current arc
        auto &adj = g[u];
        for (; cur[u] < (int)adj.size(); ++cur[u]) {
        Edge &e = adj[cur[u]];
        if (e.res > epsAbs && height[u] == height[e.to] + 1) {
            push_edge(u, e);
            if (excess[u] <= epsAbs) break;
        }
        }
        if (excess[u] <= epsAbs) break;

        // no admissible edges -> relabel
        int old = height[u];
        relabel_node(u);
        if (old < n && cnt[old] == 0) {
        gap_heuristic(old);
        }

        if (work_counter >= GR_THRESHOLD) {
        global_relabel();
        }
    }
}

// ---- Max-flow driver ----
PushRelabel::T PushRelabel::max_flow(int s_, int t_) {
    s = s_; t = t_;

    // Initialize state
    excess.assign(n, 0);
    height.assign(n, 0);
    cur.assign(n, 0);
    cnt.assign(2 * n + 2, 0);
    buckets_clear();

    // Reset residuals to original capacities (cap0), compute max capacity scale
    T maxCap = 0;
    for (int u = 0; u < n; ++u) {
        for (auto &e : g[u]) {
        e.res = e.cap0;
        if (e.cap0 > maxCap) maxCap = e.cap0;
        }
    }
    // Set absolute epsilon scaled to instance magnitude
    epsAbs = epsRel * (1.0 + maxCap);

    // Preflow: saturate all arcs out of s
    for (auto &e : g[s]) {
        if (e.res > epsAbs) {
        T f = e.res;
        e.res -= f;
        g[e.to][e.rev].res += f;
        excess[e.to] += f;
        }
    }
    // Clamp tiny noise
    for (int v = 0; v < n; ++v) clamp_zero(excess[v]);

    // Set global relabel trigger proportional to |E|
    long long m = 0;
    for (int i = 0; i < n; ++i) m += (long long)g[i].size();
    GR_THRESHOLD = max(1LL, 4LL * m);
    work_counter = 0;

    // Initial global relabel for good labels
    global_relabel();

    // Main loop: highest-label selection via buckets
    while (true) {
        int u = pop_highest();
        if (u == -1) break;                    // no active vertices
        if (excess[u] <= epsAbs) continue;     // stale
        if (height[u] != H) {                  // stale (level changed)
        activate(u);
        continue;
        }
        discharge(u);
        if (excess[u] > epsAbs) activate(u);
    }

    // Flow value: sum original capacity minus residual on arcs out of s
    T flow = 0;
    for (const auto &e : g[s]) {
        T contrib = e.cap0 - e.res;
        if (contrib > epsAbs) flow += contrib;
    }
    return flow;
}

// ---- Extract source-side min-cut as 0/1 mask ----
// S[v] = 1 iff v is reachable from s in the residual graph
void PushRelabel::mincut_source_mask(vector<unsigned char> &Smask) const {
    Smask.assign(n, 0);
    queue<int> q;
    Smask[s] = 1;
    q.push(s);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (const auto &e : g[u]) {
        if (e.res > epsAbs && !Smask[e.to]) {
            Smask[e.to] = 1;
            q.push(e.to);
        }
        }
    }
}

// Optional: compute min-cut value from a given mask (uses original capacities)
PushRelabel::T PushRelabel::mincut_value_from_mask(const vector<unsigned char> &Smask) const {
    T cut = 0;
    for (int u = 0; u < n; ++u) if (Smask[u]) {
        for (const auto &e : g[u]) if (!Smask[e.to]) {
        if (e.cap0 > epsAbs) cut += e.cap0; // sum original capacities crossing S -> ~S
        }
    }
    return cut;
}

