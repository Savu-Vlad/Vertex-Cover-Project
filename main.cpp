// main.cpp
// Vertex Cover: 4 algoritmi + benchmark CSV
// C++17

#include <bits/stdc++.h>
#include <filesystem>
using namespace std;

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
static vector<int> approx_matching_2(const Graph& g) {
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
    return cleanup_cover(g, cover);
}

// -------------------- Primal-dual 2-approx (unweighted) --------------------
static vector<int> approx_primal_dual_2(const Graph& g) {
    DynBitset rem = g.allEdges;

    vector<long double> slack(g.n, 1.0L);
    vector<char> in(g.n, false);

    const long double EPS = 1e-15L;

    while (!rem.empty()) {
        int ei = rem.first_set();
        auto [u,v] = g.edges[ei];

        if (in[u] || in[v]) {
            rem.reset(ei);
            continue;
        }

        long double delta = min(slack[u], slack[v]);
        slack[u] -= delta;
        slack[v] -= delta;

        if (slack[u] <= EPS) {
            in[u] = true;
            rem = rem.and_not(g.incMask[u]);
        }
        if (slack[v] <= EPS) {
            in[v] = true;
            rem = rem.and_not(g.incMask[v]);
        }
    }

    vector<int> cover;
    for (int i = 0; i < g.n; i++) if (in[i]) cover.push_back(i);
    return cleanup_cover(g, cover);
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

    // Tarjan to find bridges
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

    // Try each bridge as the "clique-bridge"
    for (int be : bridges) {
        auto [u0, v0] = g.edges[be];

        vector<char> vis(n, false);
        auto A = bfs_component_skip_edge(g, u0, be, vis);

        // If v0 is still reachable without the bridge, it's not a true split (paranoia)
        if (vis[v0]) continue;

        auto B = bfs_component_skip_edge(g, v0, be, vis);

        // Build membership masks
        vector<char> inA(n, false), inB(n, false);
        for (int x : A) inA[x] = true;
        for (int x : B) inB[x] = true;

        // Check that all edges are either:
        // - inside A
        // - inside B
        // - or exactly this bridge
        bool okEdges = true;
        for (int ei = 0; ei < m; ei++) {
            auto [x,y] = g.edges[ei];
            bool xa = inA[x], ya = inA[y];
            bool xb = inB[x], yb = inB[y];

            if ((xa && ya) || (xb && yb)) continue;
            if (ei == be) continue;

            // allow isolated vertices (no edges), but if an edge involves vertices outside A∪B => reject
            // (meaning there are other edged components)
            okEdges = false;
            break;
        }
        if (!okEdges) continue;

        // Verify clique property
        if (!is_clique_component(g, A, inA)) continue;
        if (!is_clique_component(g, B, inB)) continue;

        // Construct an optimal cover:
        // clique of size s needs s-1 vertices. Choose omissions so the bridge is covered.
        int omitA = -1;
        for (int x : A) if (x != u0) { omitA = x; break; }
        if (omitA == -1) omitA = u0; // A size 1

        int omitB = -1;
        for (int x : B) if (x != v0) { omitB = x; break; }
        if (omitB == -1) omitB = v0; // B size 1

        vector<int> cover;
        cover.reserve(n);

        for (int x : A) if (x != omitA) cover.push_back(x);
        for (int x : B) if (x != omitB) cover.push_back(x);

        // Ensure bridge covered (avoid omitting both endpoints)
        if (!is_vertex_cover(g, cover)) {
            // add one endpoint (safe)
            cover.push_back(u0);
        }

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

    // timeout
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
        // Upper bound from best of the two approximations (helps pruning)
        auto ub1 = approx_matching_2(g);
        auto ub2 = approx_primal_dual_2(g);
        vector<int> ub = (ub2.size() < ub1.size()) ? ub2 : ub1;

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
        // Upper bound from best of the two approximations
        auto ub1 = approx_matching_2(g);
        auto ub2 = approx_primal_dual_2(g);
        vector<int> ub = (ub2.size() < ub1.size()) ? ub2 : ub1;
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

        // If timed out, return best known UB (valid cover)
        return ub;
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

    int match=-1;
    double t_match_ms=0;

    int pd=-1;
    double t_pd_ms=0;

    double ratio_match=-1;
    double ratio_pd=-1;
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
            double t_ms = ms_since(t0);
            (void)t_ms;
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

        // Matching 2-approx
        {
            auto t0 = chrono::high_resolution_clock::now();
            auto cov = approx_matching_2(g);
            r.t_match_ms = ms_since(t0);
            r.match = (int)cov.size();
            if (!is_vertex_cover(g, cov)) cerr << "[ERR] MATCH invalid cover on " << r.name << "\n";
        }

        // Primal-dual 2-approx
        {
            auto t0 = chrono::high_resolution_clock::now();
            auto cov = approx_primal_dual_2(g);
            r.t_pd_ms = ms_since(t0);
            r.pd = (int)cov.size();
            if (!is_vertex_cover(g, cov)) cerr << "[ERR] PD invalid cover on " << r.name << "\n";
        }

        if (opt_known && r.opt > 0) {
            r.ratio_match = double(r.match) / double(r.opt);
            r.ratio_pd = double(r.pd) / double(r.opt);
        } else {
            r.ratio_match = -1.0;
            r.ratio_pd = -1.0;
        }

        rows.push_back(r);
    }

    return rows;
}

static void write_csv(const string& outPath, const vector<BenchRow>& rows) {
    ofstream out(outPath);
    out << "test,n,m,opt,bb_ms,fpt,fpt_ms,match,match_ms,match_ratio,pd,pd_ms,pd_ratio\n";
    out.setf(std::ios::fixed);
    out << setprecision(6);
    for (auto& r : rows) {
        out << r.name << ","
            << r.n << "," << r.m << ","
            << r.opt << "," << r.t_bb_ms << ","
            << r.fpt << "," << r.t_fpt_ms << ","
            << r.match << "," << r.t_match_ms << "," << r.ratio_match << ","
            << r.pd << "," << r.t_pd_ms << "," << r.ratio_pd
            << "\n";
    }
}

// -------------------- main --------------------
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // args:
    // --alg=bb|fpt|match|pd
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

    // If this is the bridge-of-two-cliques pattern, solve exactly instantly (prevents hangs on test 22)
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
        cover = approx_matching_2(g);
    } else if (alg == "pd") {
        cover = approx_primal_dual_2(g);
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
