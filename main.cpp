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
// Build with LP enabled:
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

// -------------------- Bitset --------------------
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

// -------------------- Graph --------------------
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
    for (int v : cover) in[v] = true;

    // remove redundant vertices
    for (int v : cover) {
        if (!in[v]) continue;
        in[v] = false;

        bool ok = true;
        for (int ei : g.adjEdgeIdx[v]) {
            auto [a,b] = g.edges[ei];
            if (!in[a] && !in[b]) { ok = false; break; }
        }
        if (!ok) in[v] = true;
    }

    vector<int> out;
    out.reserve(cover.size());
    for (int v = 0; v < g.n; v++) if (in[v]) out.push_back(v);
    return out;
}

// lower bound: size of a (greedy) maximal matching (matching <= OPT)
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

// -------------------- MATCH 2-approx (raw) --------------------
static vector<int> approx_matching_2_raw(const Graph& g) {
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
    for (int i = 0; i < g.n; i++) if (in[i]) cover.push_back(i);
    return cover;
}

// -------------------- LP relaxation + rounding (raw) --------------------
struct LPRawResult {
    vector<int> cover_raw;
    bool ok_opt = false;       // LP solved to optimality
    double obj = NAN;          // objective if ok_opt
};

static LPRawResult approx_lp_rounding_2_raw(const Graph& g, long long timeout_ms) {
    LPRawResult res;
#ifdef USE_GLPK
    int n = g.n;
    int m = (int)g.edges.size();
    if (m == 0) { res.cover_raw = {}; return res; }

    glp_term_out(GLP_OFF);

    glp_prob* lp = glp_create_prob();
    glp_set_prob_name(lp, "vc_lp");
    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_cols(lp, n);
    for (int v = 1; v <= n; v++) {
        glp_set_col_bnds(lp, v, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, v, 1.0);
    }

    glp_add_rows(lp, m);
    for (int i = 1; i <= m; i++) {
        glp_set_row_bnds(lp, i, GLP_LO, 1.0, 0.0);
    }

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
    if (timeout_ms > 0) parm.tm_lim = (timeout_ms > (long long)INT_MAX) ? INT_MAX : (int)timeout_ms;

    int ret = glp_simplex(lp, &parm);
    if (ret != 0) {
        glp_delete_prob(lp);
        res.cover_raw = approx_matching_2_raw(g);
        return res;
    }

    int st = glp_get_status(lp);
    if (st != GLP_OPT) {
        // daca nu e OPT, nu folosim (nici pt LB, nici pt rounding garantat)
        glp_delete_prob(lp);
        res.cover_raw = approx_matching_2_raw(g);
        return res;
    }

    res.ok_opt = true;
    res.obj = glp_get_obj_val(lp);

    vector<int> cover;
    cover.reserve(n);
    const double EPS = 1e-9;

    for (int v = 0; v < n; v++) {
        double x = glp_get_col_prim(lp, v + 1);
        if (x >= 0.5 - EPS) cover.push_back(v);
    }

    glp_delete_prob(lp);
    res.cover_raw = std::move(cover);
    return res;
#else
    (void)timeout_ms;
    res.cover_raw = approx_matching_2_raw(g);
    return res;
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
    if (m == 0) return vector<int>{};

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

    for (int i = 0; i < n; i++) if (disc[i] == -1) dfs(i);
    if (bridges.empty()) return nullopt;

    for (int be : bridges) {
        auto [u0, v0] = g.edges[be];

        vector<char> vis(n, false);
        auto A = bfs_component_skip_edge(g, u0, be, vis);
        if (vis[v0]) continue;
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

    vector<int> solve_raw() {
        // Upper bound (curat) pt pruning
        auto ub = cleanup_cover(g, approx_matching_2_raw(g));
        bestCover = ub;
        bestSize = (int)ub.size();

        vector<int> cur;
        dfs(g.allEdges, cur);
        return bestCover; // poate fi deja curat, dar ok
    }
};

// -------------------- Exact FPT (with timeout) --------------------
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

    vector<int> solve_minimum_raw() {
        auto ub_clean = cleanup_cover(g, approx_matching_2_raw(g));
        int UB = (int)ub_clean.size();

        int LB = greedy_maximal_matching_count(g, g.allEdges);
        if (LB == UB) return ub_clean;

        vector<int> cur, sol;
        for (int k = LB; k <= UB; k++) {
            cur.clear(); sol.clear();
            if (vc_decision(g.allEdges, k, cur, sol)) {
                return sol; // raw (poate avea duplicate), curatam extern
            }
            if (timed_out) break;
        }

        return ub_clean; // valid cover
    }
};

// -------------------- IO --------------------
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

// -------------------- Timing helpers --------------------
using Clock = chrono::high_resolution_clock;

static double ms_since(const Clock::time_point& t0) {
    auto t1 = Clock::now();
    return chrono::duration_cast<chrono::duration<double, std::milli>>(t1 - t0).count();
}

struct RunResult {
    vector<int> raw;
    vector<int> clean;
    double solve_ms = 0;
    double cleanup_ms = 0;
    double total_ms = 0;

    bool timed_out = false;

    // LP info (optional)
    bool lp_opt = false;
    double lp_obj = NAN;
};

template <class SolveFn>
static RunResult run_once_pipeline(const Graph& g, SolveFn&& fn, bool do_cleanup) {
    RunResult rr;

    auto t0 = Clock::now();
    rr = fn(); // must fill raw + timed_out + (lp info optional)
    rr.solve_ms = ms_since(t0);

    if (do_cleanup) {
        auto t1 = Clock::now();
        rr.clean = cleanup_cover(g, rr.raw);
        rr.cleanup_ms = ms_since(t1);
    } else {
        rr.clean = rr.raw;
        rr.cleanup_ms = 0.0;
    }

    rr.total_ms = rr.solve_ms + rr.cleanup_ms;
    return rr;
}

static RunResult median_by_total(vector<RunResult>& runs) {
    vector<int> idx(runs.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b){
        return runs[a].total_ms < runs[b].total_ms;
    });
    return std::move(runs[idx[idx.size()/2]]);
}

// -------------------- Benchmark --------------------
struct BenchRow {
    string name;
    int n=0, m=0;

    int opt=-1;

    int lb_match = 0;
    int lb_lpceil = -1;
    int lb = 0;
    string ratio_base = "na"; // opt | lb | na

    // BB
    int bb=-1; // clean size
    int bb_raw=-1;
    int bb_timedout=0;
    double bb_solve_ms=0, bb_cleanup_ms=0, bb_ms=0;

    // FPT
    int fpt=-1;
    int fpt_raw=-1;
    int fpt_timedout=0;
    double fpt_solve_ms=0, fpt_cleanup_ms=0, fpt_ms=0;

    // MATCH
    int match=-1;
    int match_raw=-1;
    double match_solve_ms=0, match_cleanup_ms=0, match_ms=0;

    // LP
    int lp=-1;
    int lp_raw=-1;
    double lp_solve_ms=0, lp_cleanup_ms=0, lp_ms=0;
    bool lp_opt=false;
    double lp_obj=NAN;

    double match_ratio=-1, match_raw_ratio=-1;
    double lp_ratio=-1, lp_raw_ratio=-1;
};

static void write_csv(const string& outPath, const vector<BenchRow>& rows) {
    ofstream out(outPath);
    out << "test,n,m,opt,lb_match,lb_lpceil,lb,ratio_base,"
           "bb_raw,bb,bb_solve_ms,bb_cleanup_ms,bb_ms,bb_timedout,"
           "fpt_raw,fpt,fpt_solve_ms,fpt_cleanup_ms,fpt_ms,fpt_timedout,"
           "match_raw,match,match_solve_ms,match_cleanup_ms,match_ms,match_raw_ratio,match_ratio,"
           "lp_raw,lp,lp_solve_ms,lp_cleanup_ms,lp_ms,lp_raw_ratio,lp_ratio\n";
    out.setf(std::ios::fixed);
    out << setprecision(6);

    for (auto& r : rows) {
        out << r.name << ","
            << r.n << "," << r.m << ","
            << r.opt << ","
            << r.lb_match << "," << r.lb_lpceil << "," << r.lb << "," << r.ratio_base << ","

            << r.bb_raw << "," << r.bb << ","
            << r.bb_solve_ms << "," << r.bb_cleanup_ms << "," << r.bb_ms << "," << r.bb_timedout << ","

            << r.fpt_raw << "," << r.fpt << ","
            << r.fpt_solve_ms << "," << r.fpt_cleanup_ms << "," << r.fpt_ms << "," << r.fpt_timedout << ","

            << r.match_raw << "," << r.match << ","
            << r.match_solve_ms << "," << r.match_cleanup_ms << "," << r.match_ms << ","
            << r.match_raw_ratio << "," << r.match_ratio << ","

            << r.lp_raw << "," << r.lp << ","
            << r.lp_solve_ms << "," << r.lp_cleanup_ms << "," << r.lp_ms << ","
            << r.lp_raw_ratio << "," << r.lp_ratio
            << "\n";
    }
}

static vector<BenchRow> run_bench(const string& folder, long long timeout_ms, int reps, bool do_cleanup) {
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
    rows.reserve(files.size());

    for (auto& path : files) {
        ifstream fin(path);
        if (!fin) continue;
        Graph g = read_graph_from_stream(fin);

        BenchRow r;
        r.name = filesystem::path(path).filename().string();
        r.n = g.n;
        r.m = (int)g.edges.size();

        // quick special-case
        optional<vector<int>> special = solve_bridge_of_two_cliques(g);

        // LB from greedy maximal matching
        r.lb_match = greedy_maximal_matching_count(g, g.allEdges);

        // --- BB (median of total time) ---
        {
            vector<RunResult> runs;
            runs.reserve(reps);
            for (int it = 0; it < reps; it++) {
                auto rr = run_once_pipeline(g, [&]() -> RunResult {
                    RunResult x;
                    BranchBoundSolver bb(g, timeout_ms);
                    x.raw = bb.solve_raw();
                    x.timed_out = bb.timed_out;
                    return x;
                }, do_cleanup);
                runs.push_back(std::move(rr));
            }
            RunResult med = median_by_total(runs);
            r.bb_raw = (int)med.raw.size();
            r.bb = (int)med.clean.size();
            r.bb_solve_ms = med.solve_ms;
            r.bb_cleanup_ms = med.cleanup_ms;
            r.bb_ms = med.total_ms;
            r.bb_timedout = med.timed_out ? 1 : 0;
        }

        // --- FPT (median of total time); OPT known if special or at least one rep finished ---
        {
            if (special.has_value()) {
                auto raw = *special;
                auto clean = do_cleanup ? cleanup_cover(g, raw) : raw;
                r.fpt_raw = (int)raw.size();
                r.fpt = (int)clean.size();
                r.fpt_solve_ms = 0.0;
                r.fpt_cleanup_ms = 0.0;
                r.fpt_ms = 0.0;
                r.fpt_timedout = 0;
                r.opt = r.fpt;
            } else {
                vector<RunResult> runs;
                runs.reserve(reps);

                int best_exact_opt = INT_MAX;
                bool any_exact = false;

                for (int it = 0; it < reps; it++) {
                    auto rr = run_once_pipeline(g, [&]() -> RunResult {
                        RunResult x;
                        FPTSolver fpt(g, timeout_ms);
                        x.raw = fpt.solve_minimum_raw();
                        x.timed_out = fpt.timed_out;
                        return x;
                    }, do_cleanup);

                    if (!rr.timed_out) {
                        any_exact = true;
                        best_exact_opt = min(best_exact_opt, (int)rr.clean.size());
                    }

                    runs.push_back(std::move(rr));
                }

                RunResult med = median_by_total(runs);
                r.fpt_raw = (int)med.raw.size();
                r.fpt = (int)med.clean.size();
                r.fpt_solve_ms = med.solve_ms;
                r.fpt_cleanup_ms = med.cleanup_ms;
                r.fpt_ms = med.total_ms;
                r.fpt_timedout = med.timed_out ? 1 : 0;

                if (any_exact) r.opt = best_exact_opt;
            }
        }

        // --- MATCH (median) ---
        {
            vector<RunResult> runs;
            runs.reserve(reps);
            for (int it = 0; it < reps; it++) {
                auto rr = run_once_pipeline(g, [&]() -> RunResult {
                    RunResult x;
                    x.raw = approx_matching_2_raw(g);
                    return x;
                }, do_cleanup);
                runs.push_back(std::move(rr));
            }
            RunResult med = median_by_total(runs);
            r.match_raw = (int)med.raw.size();
            r.match = (int)med.clean.size();
            r.match_solve_ms = med.solve_ms;
            r.match_cleanup_ms = med.cleanup_ms;
            r.match_ms = med.total_ms;
        }

        // --- LP (median) + LP lower bound if OPT ---
        {
            vector<RunResult> runs;
            runs.reserve(reps);

            bool saw_lp_opt = false;
            double lp_obj_opt = NAN;

            for (int it = 0; it < reps; it++) {
                auto rr = run_once_pipeline(g, [&]() -> RunResult {
                    RunResult x;
                    auto lp = approx_lp_rounding_2_raw(g, timeout_ms);
                    x.raw = std::move(lp.cover_raw);
                    x.lp_opt = lp.ok_opt;
                    x.lp_obj = lp.obj;
                    return x;
                }, do_cleanup);
                if (rr.lp_opt) { saw_lp_opt = true; lp_obj_opt = rr.lp_obj; }
                runs.push_back(std::move(rr));
            }

            RunResult med = median_by_total(runs);
            r.lp_raw = (int)med.raw.size();
            r.lp = (int)med.clean.size();
            r.lp_solve_ms = med.solve_ms;
            r.lp_cleanup_ms = med.cleanup_ms;
            r.lp_ms = med.total_ms;

            r.lp_opt = saw_lp_opt;
            r.lp_obj = lp_obj_opt;

            if (saw_lp_opt && std::isfinite(lp_obj_opt)) {
                // ceil(lp_obj) as integer LB (LP optimum is a LB on OPT)
                int lb_lp = (int)ceil(lp_obj_opt - 1e-9);
                r.lb_lpceil = lb_lp;
            } else {
                r.lb_lpceil = -1;
            }
        }

        // finalize LB
        r.lb = r.lb_match;
        if (r.lb_lpceil >= 0) r.lb = max(r.lb, r.lb_lpceil);

        // ratios: prefer OPT if known (>0), else LB (>0), else -1
        int denom = -1;
        if (r.opt > 0) {
            denom = r.opt;
            r.ratio_base = "opt";
        } else if (r.lb > 0) {
            denom = r.lb;
            r.ratio_base = "lb";
        } else {
            r.ratio_base = "na";
        }

        if (denom > 0) {
            r.match_raw_ratio = (double)r.match_raw / (double)denom;
            r.match_ratio = (double)r.match / (double)denom;
            r.lp_raw_ratio = (double)r.lp_raw / (double)denom;
            r.lp_ratio = (double)r.lp / (double)denom;
        } else {
            r.match_raw_ratio = r.match_ratio = -1.0;
            r.lp_raw_ratio = r.lp_ratio = -1.0;
        }

        // sanity
        if (!is_vertex_cover(g, do_cleanup ? cleanup_cover(g, { }) : vector<int>{})) {
            // ignore
        }
        if (!is_vertex_cover(g, cleanup_cover(g, approx_matching_2_raw(g)))) {
            cerr << "[ERR] internal sanity failed on " << r.name << "\n";
        }

        rows.push_back(std::move(r));
    }

    return rows;
}

// -------------------- main --------------------
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // args:
    // --alg=bb|fpt|match|lp
    // --bench <folder> --out <csv>
    // --timeout_ms <int>   (default 2000)
    // --reps <int>         (default 5)
    // --no_cleanup         (skip cleanup_cover; raport raw==clean)
    string alg = "bb";
    bool bench = false;
    string benchFolder = "tests";
    string outCsv = "results.csv";
    long long timeout_ms = 2000;
    int reps = 5;
    bool do_cleanup = true;

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
        } else if (s == "--reps") {
            if (i + 1 < argc) reps = max(1, stoi(argv[++i]));
        } else if (s == "--no_cleanup") {
            do_cleanup = false;
        }
    }

    if (bench) {
        auto rows = run_bench(benchFolder, timeout_ms, reps, do_cleanup);
        write_csv(outCsv, rows);
        cerr << "Wrote CSV: " << outCsv << " (" << rows.size() << " tests)\n";
        return 0;
    }

    Graph g = read_graph_from_stream(cin);

    // special-case exact
    if (auto special = solve_bridge_of_two_cliques(g); special.has_value()) {
        auto cover = do_cleanup ? cleanup_cover(g, *special) : *special;
        print_solution(cover);
        return 0;
    }

    vector<int> raw, clean;

    if (alg == "bb") {
        BranchBoundSolver bb(g, timeout_ms);
        raw = bb.solve_raw();
    } else if (alg == "fpt") {
        FPTSolver fpt(g, timeout_ms);
        raw = fpt.solve_minimum_raw();
    } else if (alg == "match") {
        raw = approx_matching_2_raw(g);
    } else if (alg == "lp") {
        raw = approx_lp_rounding_2_raw(g, timeout_ms).cover_raw;
    } else {
        BranchBoundSolver bb(g, timeout_ms);
        raw = bb.solve_raw();
    }

    clean = do_cleanup ? cleanup_cover(g, raw) : raw;

    if (!is_vertex_cover(g, clean)) {
        cerr << "[ERR] produced invalid cover\n";
    }

    print_solution(clean);
    return 0;
}
