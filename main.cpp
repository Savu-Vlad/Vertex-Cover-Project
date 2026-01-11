// main.cpp
// Vertex Cover: 4 algoritmi + benchmark CSV
// C++17
//
// Algs:
//  - bb   : Branch & Bound exact (cu timeout)
//  - fpt  : FPT (Buss kernel + bounded search) exact (cu timeout)
//  - match: 2-approx via maximal matching
//  - lp   : LP relax + rounding (2-approx) via GLPK (optional)
//
// Build with LP enabled (Linux):
//   sudo apt-get install -y libglpk-dev
//   g++ -O2 -std=c++17 -DUSE_GLPK main.cpp -o vc -lglpk
//
// Build without LP (lp falls back to match):
//   g++ -O2 -std=c++17 main.cpp -o vc

#include <bits/stdc++.h>
#include <filesystem>
using namespace std;

#ifdef USE_GLPK
#include <glpk.h>
#endif

struct DynBitset {
    int nbits = 0;
    vector<uint64_t> b;

    DynBitset() = default;
    explicit DynBitset(int n) { init(n); }

    void init(int n) {
        nbits = n;
        b.assign((n + 63) / 64, 0ULL);
    }

    inline void trim_last() {
        if (b.empty()) return;
        int r = nbits & 63;
        if (r == 0) return;
        uint64_t mask = (r == 64) ? ~0ULL : ((1ULL << r) - 1ULL);
        b.back() &= mask;
    }

    inline void set(int i) { b[i >> 6] |= (1ULL << (i & 63)); }
    inline void reset(int i) { b[i >> 6] &= ~(1ULL << (i & 63)); }
    inline bool test(int i) const { return (b[i >> 6] >> (i & 63)) & 1ULL; }

    inline bool empty() const {
        for (auto x : b) if (x) return false;
        return true;
    }

    inline int count() const {
        int c = 0;
        for (auto x : b) c += __builtin_popcountll(x);
        return c;
    }

    // first set bit index, -1 if none
    inline int first_set() const {
        for (int bi = 0; bi < (int)b.size(); bi++) {
            uint64_t x = b[bi];
            if (!x) continue;
            int t = __builtin_ctzll(x);
            int idx = (bi << 6) + t;
            if (idx < nbits) return idx;
        }
        return -1;
    }

    inline DynBitset and_not(const DynBitset& mask) const {
        DynBitset res;
        res.nbits = nbits;
        res.b.resize(b.size());
        for (size_t i = 0; i < b.size(); i++) res.b[i] = b[i] & ~mask.b[i];
        res.trim_last();
        return res;
    }

    inline DynBitset bit_or(const DynBitset& other) const {
        DynBitset res;
        res.nbits = nbits;
        res.b.resize(b.size());
        for (size_t i = 0; i < b.size(); i++) res.b[i] = b[i] | other.b[i];
        return res;
    }

    template <class F>
    inline void for_each_setbit(F&& f) const {
        for (int bi = 0; bi < (int)b.size(); bi++) {
            uint64_t x = b[bi];
            while (x) {
                uint64_t lsb = x & -x;
                int t = __builtin_ctzll(x);
                int idx = (bi << 6) + t;
                if (idx >= nbits) break;
                f(idx);
                x ^= lsb;
            }
        }
    }
};

static inline int intersection_count(const DynBitset& a, const DynBitset& c) {
    int s = 0;
    for (size_t i = 0; i < a.b.size(); i++) {
        s += __builtin_popcountll(a.b[i] & c.b[i]);
    }
    return s;
}

struct Graph {
    int n = 0;
    vector<pair<int,int>> edges;      // index -> (u,v)
    vector<vector<int>> adjEdgeIdx;   // v -> list of edge indices
    vector<DynBitset> incMask;        // v -> bitset over edges indices incident to v
    DynBitset allEdges;               // bitset with all edges set

    void build(int N, const vector<pair<int,int>>& inEdges) {
        n = N;
        edges = inEdges;
        int m = (int)edges.size();

        adjEdgeIdx.assign(n, {});
        incMask.assign(n, DynBitset(m));
        for (int i = 0; i < n; i++) incMask[i].init(m);

        for (int ei = 0; ei < m; ei++) {
            auto [u,v] = edges[ei];
            adjEdgeIdx[u].push_back(ei);
            adjEdgeIdx[v].push_back(ei);
            incMask[u].set(ei);
            incMask[v].set(ei);
        }

        allEdges.init(m);
        for (int ei = 0; ei < m; ei++) allEdges.set(ei);
    }
};

static inline bool is_vertex_cover(const Graph& g, const vector<int>& cover) {
    vector<char> in(g.n, false);
    for (int v : cover) if (0 <= v && v < g.n) in[v] = true;
    for (auto [u,v] : g.edges) {
        if (!in[u] && !in[v]) return false;
    }
    return true;
}

static inline vector<int> cleanup_cover(const Graph& g, const vector<int>& coverIn) {
    vector<char> in(g.n, false);
    vector<int> cover = coverIn;
    sort(cover.begin(), cover.end());
    cover.erase(unique(cover.begin(), cover.end()), cover.end());
    for (int v : cover) if (0 <= v && v < g.n) in[v] = true;

    bool changed = true;
    while (changed) {
        changed = false;

        // Important: iterate all vertices; order matters less with repetition
        for (int v = 0; v < g.n; v++) {
            if (!in[v]) continue;

            in[v] = false;
            bool ok = true;

            for (int ei : g.adjEdgeIdx[v]) {
                auto [a, b] = g.edges[ei];
                if (!in[a] && !in[b]) { ok = false; break; }
            }

            if (!ok) {
                in[v] = true;
            } else {
                changed = true; // we successfully removed v
            }
        }
    }

    vector<int> out;
    for (int v = 0; v < g.n; v++) if (in[v]) out.push_back(v);
    return out;
}


static inline int greedy_maximal_matching_count(const Graph& g, const DynBitset& rem) {
    vector<char> used(g.n, false);
    int cnt = 0;
    rem.for_each_setbit([&](int ei){
        auto [u,v] = g.edges[ei];
        if (!used[u] && !used[v]) {
            used[u] = used[v] = true;
            cnt++;
        }
    });
    return cnt;
}

// -------------------- Approx 2-approx via maximal matching (edges) --------------------
static vector<int> approx_matching_2(const Graph& g, bool do_cleanup = true) {
    DynBitset rem = g.allEdges;
    vector<char> in(g.n, false);

    while (!rem.empty()) {
        int ei = rem.first_set();
        auto [u,v] = g.edges[ei];
        in[u] = true; in[v] = true;

        DynBitset mask = g.incMask[u].bit_or(g.incMask[v]);
        rem = rem.and_not(mask);
    }

    vector<int> cover;
    cover.reserve(g.n);
    for (int i = 0; i < g.n; i++) if (in[i]) cover.push_back(i);

    if (do_cleanup) return cleanup_cover(g, cover);
    return cover; // already unique
}


// -------------------- LP relaxation + rounding (2-approx) --------------------
static vector<int> approx_lp_rounding_2(const Graph& g, long long timeout_ms, bool do_cleanup = true) {
#ifdef USE_GLPK
    int n = g.n;
    int m = (int)g.edges.size();
    if (m == 0) return {};

    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "vc_lp");
    glp_set_obj_dir(lp, GLP_MIN);

    // variables x_v
    glp_add_cols(lp, n);
    for (int v = 1; v <= n; v++) {
        glp_set_col_bnds(lp, v, GLP_DB, 0.0, 1.0); // 0 <= x_v <= 1
        glp_set_obj_coef(lp, v, 1.0);              // min sum x_v
    }

    // constraints x_u + x_v >= 1 for each edge
    glp_add_rows(lp, m);
    for (int i = 1; i <= m; i++) {
        glp_set_row_bnds(lp, i, GLP_LO, 1.0, 0.0);
    }

    // matrix (2 nonzeros per row)
    vector<int> ia(1 + 2*m), ja(1 + 2*m);
    vector<double> ar(1 + 2*m);

    for (int i = 0; i < m; i++) {
        auto [u, v] = g.edges[i];
        int row = i + 1;
        int k1 = 2*i + 1;
        int k2 = 2*i + 2;

        ia[k1] = row; ja[k1] = u + 1; ar[k1] = 1.0;
        ia[k2] = row; ja[k2] = v + 1; ar[k2] = 1.0;
    }
    glp_load_matrix(lp, 2*m, ia.data(), ja.data(), ar.data());

    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.presolve = GLP_ON;

    if (timeout_ms > 0) {
        parm.tm_lim = (timeout_ms > (long long)INT_MAX) ? INT_MAX : (int)timeout_ms;
    }

    int ret = glp_simplex(lp, &parm);
    if (ret != 0) {
        glp_delete_prob(lp);
        return approx_matching_2(g, do_cleanup);
    }

    vector<int> cover;
    cover.reserve(n);
    const double EPS = 1e-9;

    // rounding: x_v >= 0.5 => v in cover
    for (int v = 0; v < n; v++) {
        double x = glp_get_col_prim(lp, v + 1);
        if (x >= 0.5 - EPS) cover.push_back(v);
    }

    glp_delete_prob(lp);

    if (do_cleanup) return cleanup_cover(g, cover);
    return cover;
#else
    (void)timeout_ms;
    return approx_matching_2(g, do_cleanup);
#endif
}


// -------------------- Quick exact solver: two cliques connected by ONE bridge --------------------
static vector<int> bfs_component_skip_edge(const Graph& g, int start, int banned_ei, vector<char>& vis) {
    vector<int> comp;
    queue<int> q;
    vis[start] = true;
    q.push(start);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        comp.push_back(u);
        for (int ei : g.adjEdgeIdx[u]) {
            if (ei == banned_ei) continue;
            auto [a,b] = g.edges[ei];
            int v = (a == u) ? b : a;
            if (!vis[v]) { vis[v] = true; q.push(v); }
        }
    }
    return comp;
}

static bool is_clique_component(const Graph& g, const vector<int>& comp, const vector<char>& inComp) {
    long long s = (long long)comp.size();
    if (s <= 1) return true;
    long long required = s * (s - 1) / 2;
    long long cnt = 0;
    for (auto [u,v] : g.edges) {
        if (inComp[u] && inComp[v]) cnt++;
    }
    return cnt == required;
}

static optional<vector<int>> solve_bridge_of_two_cliques(const Graph& g) {
    int n = g.n;
    int m = (int)g.edges.size();
    if (m == 0) return vector<int>{}; // empty graph

    // Tarjan bridges
    vector<int> disc(n, -1), low(n, -1), parentEdge(n, -1);
    int timer = 0;
    vector<int> bridges;

    function<void(int)> dfs = [&](int u) {
        disc[u] = low[u] = ++timer;
        for (int ei : g.adjEdgeIdx[u]) {
            auto [a,b] = g.edges[ei];
            int v = (a == u) ? b : a;
            if (disc[v] == -1) {
                parentEdge[v] = ei;
                dfs(v);
                low[u] = min(low[u], low[v]);
                if (low[v] > disc[u]) bridges.push_back(ei);
            } else if (ei != parentEdge[u]) {
                low[u] = min(low[u], disc[v]);
            }
        }
    };

    for (int i = 0; i < n; i++) {
        if (disc[i] == -1) dfs(i);
    }
    if (bridges.empty()) return nullopt;

    for (int be : bridges) {
        auto [u0, v0] = g.edges[be];

        vector<char> vis(n, false);
        auto A = bfs_component_skip_edge(g, u0, be, vis);
        if (vis[v0]) continue; // not a true split
        auto B = bfs_component_skip_edge(g, v0, be, vis);

        vector<char> inA(n, false), inB(n, false);
        for (int x : A) inA[x] = true;
        for (int x : B) inB[x] = true;

        bool okEdges = true;
        for (int ei = 0; ei < m; ei++) {
            auto [x,y] = g.edges[ei];
            bool xa = inA[x], ya = inA[y];
            bool xb = inB[x], yb = inB[y];

            if ((xa && ya) || (xb && yb)) continue;
            if (ei == be) continue;

            okEdges = false;
            break;
        }
        if (!okEdges) continue;

        if (!is_clique_component(g, A, inA)) continue;
        if (!is_clique_component(g, B, inB)) continue;

        // optimal cover: (|A|-1) + (|B|-1), avoid omitting both bridge endpoints
        int omitA = -1;
        for (int x : A) if (x != u0) { omitA = x; break; }
        if (omitA == -1) omitA = u0;

        int omitB = -1;
        for (int x : B) if (x != v0) { omitB = x; break; }
        if (omitB == -1) omitB = v0;

        vector<int> cover;
        cover.reserve(n);

        for (int x : A) if (x != omitA) cover.push_back(x);
        for (int x : B) if (x != omitB) cover.push_back(x);

        if (!is_vertex_cover(g, cover)) cover.push_back(u0);

        cover = cleanup_cover(g, cover);
        if (is_vertex_cover(g, cover)) return cover;
    }

    return nullopt;
}

// -------------------- Exact Branch & Bound (with timeout) --------------------
struct BranchBoundSolver {
    const Graph& g;
    vector<int> bestCover;
    int bestSize;

    bool use_timeout = false;
    bool timed_out = false;
    chrono::high_resolution_clock::time_point deadline;
    uint64_t calls = 0;

    explicit BranchBoundSolver(const Graph& graph, long long timeout_ms) : g(graph) {
        bestSize = g.n + 1;
        if (timeout_ms > 0) {
            use_timeout = true;
            deadline = chrono::high_resolution_clock::now() + chrono::milliseconds(timeout_ms);
        }
    }

    inline void check_timeout() {
        if (!use_timeout || timed_out) return;
        if ((++calls & 0x3FFFULL) == 0) {
            if (chrono::high_resolution_clock::now() >= deadline) timed_out = true;
        }
    }

    void dfs(const DynBitset& rem, vector<int>& cur) {
        check_timeout();
        if (timed_out) return;

        if ((int)cur.size() >= bestSize) return;

        int lb = greedy_maximal_matching_count(g, rem);
        if ((int)cur.size() + lb >= bestSize) return;

        if (rem.empty()) {
            bestSize = (int)cur.size();
            bestCover = cur;
            return;
        }

        int ei = rem.first_set();
        auto [u,v] = g.edges[ei];

        int degU = intersection_count(rem, g.incMask[u]);
        int degV = intersection_count(rem, g.incMask[v]);

        int first = (degU >= degV) ? u : v;
        int second = (first == u) ? v : u;

        {
            cur.push_back(first);
            DynBitset rem2 = rem.and_not(g.incMask[first]);
            dfs(rem2, cur);
            cur.pop_back();
            if (timed_out) return;
        }
        {
            cur.push_back(second);
            DynBitset rem2 = rem.and_not(g.incMask[second]);
            dfs(rem2, cur);
            cur.pop_back();
        }
    }

    vector<int> solve() {
        // Upper bound from matching (fast, good enough for pruning)
        auto ub = approx_matching_2(g);
        bestCover = ub;
        bestSize = (int)ub.size();

        vector<int> cur;
        dfs(g.allEdges, cur);

        return cleanup_cover(g, bestCover);
    }
};

// -------------------- Exact FPT: Buss kernel + bounded search, iterative deepening (with timeout) --------------------
struct FPTSolver {
    const Graph& g;

    bool use_timeout = false;
    bool timed_out = false;
    chrono::high_resolution_clock::time_point deadline;
    uint64_t calls = 0;

    explicit FPTSolver(const Graph& graph, long long timeout_ms) : g(graph) {
        if (timeout_ms > 0) {
            use_timeout = true;
            deadline = chrono::high_resolution_clock::now() + chrono::milliseconds(timeout_ms);
        }
    }

    inline void check_timeout() {
        if (!use_timeout || timed_out) return;
        if ((++calls & 0x3FFFULL) == 0) {
            if (chrono::high_resolution_clock::now() >= deadline) timed_out = true;
        }
    }

    bool apply_reductions(DynBitset& rem, int& k, vector<int>& forced) {
        while (true) {
            check_timeout();
            if (timed_out) return false;

            if (k < 0) return false;
            if (rem.empty()) return true;

            vector<int> deg(g.n, 0);
            for (int v = 0; v < g.n; v++) deg[v] = intersection_count(rem, g.incMask[v]);

            bool changed = false;

            // rule: deg(v) > k => v forced
            for (int v = 0; v < g.n; v++) {
                if (deg[v] > k) {
                    forced.push_back(v);
                    rem = rem.and_not(g.incMask[v]);
                    k--;
                    changed = true;
                    break;
                }
            }
            if (changed) continue;

            // rule: deg(v) == 1 => include its unique neighbor
            for (int v = 0; v < g.n; v++) {
                if (deg[v] == 1) {
                    int eidx = -1;
                    for (int ei : g.adjEdgeIdx[v]) {
                        if (rem.test(ei)) { eidx = ei; break; }
                    }
                    if (eidx == -1) continue;
                    auto [a,b] = g.edges[eidx];
                    int u = (a == v) ? b : a;

                    forced.push_back(u);
                    rem = rem.and_not(g.incMask[u]);
                    k--;
                    changed = true;
                    break;
                }
            }
            if (changed) continue;

            // Buss kernel edge bound
            long long mrem = rem.count();
            if (mrem > 1LL * k * k) return false;

            return true;
        }
    }

    bool vc_decision(DynBitset rem, int k, vector<int>& cur, vector<int>& sol) {
        check_timeout();
        if (timed_out) return false;

        vector<int> forced;
        int k2 = k;

        if (!apply_reductions(rem, k2, forced)) return false;
        if (timed_out) return false;

        int forcedCnt = (int)forced.size();
        for (int v : forced) cur.push_back(v);

        if (rem.empty()) {
            sol = cur;
            while (forcedCnt--) cur.pop_back();
            return true;
        }

        if (k2 == 0) {
            while (forcedCnt--) cur.pop_back();
            return false;
        }

        int ei = rem.first_set();
        auto [u,v] = g.edges[ei];

        // branch include u
        cur.push_back(u);
        if (vc_decision(rem.and_not(g.incMask[u]), k2 - 1, cur, sol)) {
            cur.pop_back();
            while (forcedCnt--) cur.pop_back();
            return true;
        }
        cur.pop_back();

        check_timeout();
        if (timed_out) { while (forcedCnt--) cur.pop_back(); return false; }

        // branch include v
        cur.push_back(v);
        if (vc_decision(rem.and_not(g.incMask[v]), k2 - 1, cur, sol)) {
            cur.pop_back();
            while (forcedCnt--) cur.pop_back();
            return true;
        }
        cur.pop_back();

        while (forcedCnt--) cur.pop_back();
        return false;
    }

    vector<int> solve_minimum() {
        // Upper bound from matching
        auto ub = approx_matching_2(g);
        int UB = (int)ub.size();

        int LB = greedy_maximal_matching_count(g, g.allEdges);
        if (LB == UB) return ub;

        vector<int> cur, sol;
        for (int k = LB; k <= UB; k++) {
            cur.clear(); sol.clear();
            if (vc_decision(g.allEdges, k, cur, sol)) {
                return cleanup_cover(g, sol);
            }
            if (timed_out) break;
        }

        return ub; // valid cover
    }
};

// -------------------- IO helpers --------------------
static Graph read_graph_from_stream(istream& in) {
    int N, M;
    in >> N >> M;
    vector<pair<int,int>> ed;
    ed.reserve(M);

    unordered_set<uint64_t> seen;
    seen.reserve((size_t)M * 2);

    for (int i = 0; i < M; i++) {
        int x, y;
        in >> x >> y;
        if (x == y) continue;
        int a = min(x,y), b = max(x,y);
        uint64_t key = (uint64_t(a) << 32) ^ uint64_t(b);
        if (seen.insert(key).second) ed.push_back({a,b});
    }

    Graph g;
    g.build(N, ed);
    return g;
}

static void print_solution(const vector<int>& cover) {
    cout << cover.size() << "\n";
    for (size_t i = 0; i < cover.size(); i++) {
        if (i) cout << ' ';
        cout << cover[i];
    }
    cout << "\n";
}

static string pad2(int x) {
    ostringstream oss;
    oss << setw(2) << setfill('0') << x;
    return oss.str();
}

// -------------------- Benchmark --------------------
struct BenchRow {
    string name;
    int n=0, m=0;

    int opt=-1;

    double t_bb_ms=0;

    int fpt=-1;
    double t_fpt_ms=0;

    int match_raw=-1;
    int match=-1;            // CLEAN
    double t_match_ms=0;
    double ratio_match_raw=-1;
    double ratio_match=-1;   // CLEAN ratio

    int lp_raw=-1;
    int lp=-1;               // CLEAN
    double t_lp_ms=0;
    double ratio_lp_raw=-1;
    double ratio_lp=-1;      // CLEAN ratio
};


static double ms_since(const chrono::high_resolution_clock::time_point& t0) {
    auto t1 = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::duration<double, std::milli>>(t1 - t0).count();
}

static vector<BenchRow> run_bench(const string& folder, long long timeout_ms) {
    vector<string> files;

    for (int i = 1; i <= 99; i++) {
        string fn = folder + "/" + pad2(i) + ".in";
        if (filesystem::exists(fn)) files.push_back(fn);
    }
    if (files.empty()) {
        for (auto& p : filesystem::directory_iterator(folder)) {
            if (!p.is_regular_file()) continue;
            if (p.path().extension() == ".in") files.push_back(p.path().string());
        }
        sort(files.begin(), files.end());
    }

    vector<BenchRow> rows;

    for (auto& path : files) {
        ifstream fin(path);
        if (!fin) continue;
        Graph g = read_graph_from_stream(fin);

        BenchRow r;
        r.name = filesystem::path(path).filename().string();
        r.n = g.n;
        r.m = (int)g.edges.size();

        // Quick exact special-case: two cliques + one bridge
        optional<vector<int>> special;
        {
            auto t0 = chrono::high_resolution_clock::now();
            special = solve_bridge_of_two_cliques(g);
            (void)ms_since(t0);
        }

        // BB (timeout)
        {
            auto t0 = chrono::high_resolution_clock::now();
            BranchBoundSolver bb(g, timeout_ms);
            auto cov = bb.solve();
            r.t_bb_ms = ms_since(t0);
            if (!is_vertex_cover(g, cov)) cerr << "[ERR] BB invalid cover on " << r.name << "\n";
            if (bb.timed_out) cerr << "[INFO] BB timeout on " << r.name << "\n";
        }

        // FPT exact baseline (timeout) OR special-case exact
        bool opt_known = false;
        if (special.has_value()) {
            auto cov = cleanup_cover(g, *special);
            r.fpt = (int)cov.size();
            r.t_fpt_ms = 0.0;
            r.opt = r.fpt;
            opt_known = true;
        } else {
            auto t0 = chrono::high_resolution_clock::now();
            FPTSolver fpt(g, timeout_ms);
            auto cov = fpt.solve_minimum();
            r.t_fpt_ms = ms_since(t0);
            r.fpt = (int)cov.size();
            if (!is_vertex_cover(g, cov)) cerr << "[ERR] FPT invalid cover on " << r.name << "\n";
            if (fpt.timed_out) {
                cerr << "[INFO] FPT timeout on " << r.name << " (OPT unknown)\n";
                r.opt = -1;
                opt_known = false;
            } else {
                r.opt = r.fpt; // OPT exact from FPT
                opt_known = true;
            }
        }

        // Matching 2-approx (RAW + CLEAN)
        {
            auto t0 = chrono::high_resolution_clock::now();
            auto cov_raw = approx_matching_2(g, false);
            r.t_match_ms = ms_since(t0);

            r.match_raw = (int)cov_raw.size();
            auto cov = cleanup_cover(g, cov_raw);
            r.match = (int)cov.size();

            if (!is_vertex_cover(g, cov)) cerr << "[ERR] MATCH invalid cover on " << r.name << "\n";
        }


        // LP rounding 2-approx (RAW + CLEAN)
        {
            auto t0 = chrono::high_resolution_clock::now();
            auto cov_raw = approx_lp_rounding_2(g, timeout_ms, false);
            r.t_lp_ms = ms_since(t0);

            r.lp_raw = (int)cov_raw.size();
            auto cov = cleanup_cover(g, cov_raw);
            r.lp = (int)cov.size();

            if (!is_vertex_cover(g, cov)) cerr << "[ERR] LP invalid cover on " << r.name << "\n";
        }


        if (opt_known && r.opt > 0) {
            r.ratio_match_raw = double(r.match_raw) / double(r.opt);
            r.ratio_match     = double(r.match)     / double(r.opt);

            r.ratio_lp_raw    = double(r.lp_raw)    / double(r.opt);
            r.ratio_lp        = double(r.lp)        / double(r.opt);
        } else {
            r.ratio_match_raw = -1.0;
            r.ratio_match     = -1.0;
            r.ratio_lp_raw    = -1.0;
            r.ratio_lp        = -1.0;
        }


        rows.push_back(r);
    }

    return rows;
}

static void write_csv(const string& outPath, const vector<BenchRow>& rows) {
    ofstream out(outPath);
    out << "test,n,m,opt,bb_ms,fpt,fpt_ms,"
           "match_raw,match,match_ms,match_raw_ratio,match_ratio,"
           "lp_raw,lp,lp_ms,lp_raw_ratio,lp_ratio\n";
    out.setf(std::ios::fixed);
    out << setprecision(6);

    for (auto& r : rows) {
        out << r.name << ","
            << r.n << "," << r.m << ","
            << r.opt << "," << r.t_bb_ms << ","
            << r.fpt << "," << r.t_fpt_ms << ","
            << r.match_raw << "," << r.match << "," << r.t_match_ms << ","
            << r.ratio_match_raw << "," << r.ratio_match << ","
            << r.lp_raw << "," << r.lp << "," << r.t_lp_ms << ","
            << r.ratio_lp_raw << "," << r.ratio_lp
            << "\n";
    }
}


// -------------------- main --------------------
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // args:
    // --alg=bb|fpt|match|lp
    // --bench <folder> --out <csv>
    // --timeout_ms <int>   (default 2000)
    string alg = "bb";
    bool bench = false;
    string benchFolder = "tests";
    string outCsv = "results.csv";
    long long timeout_ms = 2000; // default 2 seconds

    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s.rfind("--alg=", 0) == 0) {
            alg = s.substr(6);
        } else if (s == "--bench") {
            bench = true;
            if (i + 1 < argc) benchFolder = argv[++i];
        } else if (s == "--out") {
            if (i + 1 < argc) outCsv = argv[++i];
        } else if (s == "--timeout_ms") {
            if (i + 1 < argc) timeout_ms = stoll(argv[++i]);
        }
    }

    if (bench) {
        auto rows = run_bench(benchFolder, timeout_ms);
        write_csv(outCsv, rows);
        cerr << "Wrote CSV: " << outCsv << " (" << rows.size() << " tests)\n";
        return 0;
    }

    Graph g = read_graph_from_stream(cin);

    // If this is the bridge-of-two-cliques pattern, solve exactly instantly
    if (auto special = solve_bridge_of_two_cliques(g); special.has_value()) {
        auto cover = cleanup_cover(g, *special);
        print_solution(cover);
        return 0;
    }

    vector<int> cover;

    if (alg == "bb") {
        BranchBoundSolver bb(g, timeout_ms);
        cover = bb.solve();
    } else if (alg == "fpt") {
        FPTSolver fpt(g, timeout_ms);
        cover = fpt.solve_minimum();
    } else if (alg == "match") {
        cover = approx_matching_2(g, false);
    } else if (alg == "lp") {
        cover = approx_lp_rounding_2(g, timeout_ms, false);
    } else {
        BranchBoundSolver bb(g, timeout_ms);
        cover = bb.solve();
    }

    cover = cleanup_cover(g, cover);

    if (!is_vertex_cover(g, cover)) {
        cerr << "[ERR] produced invalid cover\n";
    }

    print_solution(cover);
    return 0;
}
