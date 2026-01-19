#!/usr/bin/env python3
"""
A small script I wrote to spit out 25 test inputs for the vertex-cover
problem. Each file goes into `tests/` and is named `01.in`..`25.in`.

File format: first line `N M`, then `M` lines `u v` with 0-based node ids.

I included a mix of cases: tiny examples, paths, cycles, cliques, stars,
some trees/grids, bipartite examples, a couple of hub-like graphs and a
few larger instances so you can check performance. It's deterministic so
it's easy to reproduce the same set later.

Run: `python3 scripts/generate_tests.py`
"""
from pathlib import Path
from itertools import combinations
import random

OUTDIR = Path(__file__).resolve().parent.parent / "tests"
OUTDIR.mkdir(exist_ok=True)

def write_test(name, n, edges):
    path = OUTDIR / name
    with path.open("w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u,v in edges:
            f.write(f"{u} {v}\n")
    print(f"Wrote {path} (n={n}, m={len(edges)})")

# 01: one node, no edges
write_test("01.in", 1, [])

# 02: two nodes with one edge
write_test("02.in", 2, [(0,1)])

# 03: tiny path (3 nodes)
write_test("03.in", 3, [(0,1),(1,2)])

# 04: small cycle (4 nodes)
write_test("04.in", 4, [(0,1),(1,2),(2,3),(3,0)])

# 05: complete graph K5 (all pairs)
n=5
edges=[(i,j) for i in range(n) for j in range(i+1,n)]
write_test("05.in", n, edges)

# 06: small star (center 0, 9 leaves)
n=10
edges=[(0,i) for i in range(1,n)]
write_test("06.in", n, edges)

# 07: path of length 20 (simple long-ish chain)
n=20
edges=[(i,i+1) for i in range(n-1)]
write_test("07.in", n, edges)

# 08: 20-node cycle
n=20
edges=[(i,(i+1)%n) for i in range(n)]
write_test("08.in", n, edges)

# 09: complete bipartite K10,10 (left 0..9, right 10..19)
n=20
edges=[(i,j) for i in range(0,10) for j in range(10,20)]
write_test("09.in", n, edges)

# 10: mostly a path (N=100) with a few short chords added
n=100
edges=[(i,i+1) for i in range(n-1)]
# add a few (i, i+2) shortcuts to make it slightly less trivial
for i in range(0,21):
    edges.append((i,i+2))
write_test("10.in", n, edges)

# 11: denser-ish graph on 100 nodes — first 2000 pairs in lexicographic order
n=100
edges=[]
count=2000
for i in range(n):
    for j in range(i+1,n):
        edges.append((i,j))
        if len(edges)>=count:
            break
    if len(edges)>=count:
        break
write_test("11.in", n, edges)

# 12: three disconnected components in one file:
#   - path of 5 nodes (0..4)
#   - cycle of 7 nodes (5..11)
#   - star with center 12 and leaves 13..19
edges=[]
edges += [(i,i+1) for i in range(0,4)]
edges += [(i,i+1) for i in range(5,11)]
edges.append((11,5))
edges += [(12,i) for i in range(13,20)]
write_test("12.in", 20, edges)

# 13: a binary tree layout on 50 nodes (parent = (i-1)//2)
n=50
edges=[((i-1)//2, i) for i in range(1,n)]
write_test("13.in", n, edges)

# 14: 10x10 grid (nice for locality tests)
rows,cols = 10,10
n = rows*cols
edges=[]
for r in range(rows):
    for c in range(cols):
        u = r*cols + c
        if c+1<cols:
            edges.append((u, r*cols + (c+1)))
        if r+1<rows:
            edges.append((u, (r+1)*cols + c))
write_test("14.in", n, edges)

# 15: clique K15 (all pairs)
n=15
edges=[(i,j) for i in range(n) for j in range(i+1,n)]
write_test("15.in", n, edges)

# 16: sparse bipartite-ish graph (N=200, ~150 edges) — deterministic spread
n=200
edges=[]
for idx in range(150):
    left = idx % 100
    right = 100 + ((idx*7) % 100)
    edges.append((left, right))
# remove duplicates but keep the order
seen=set()
uniq=[]
for u,v in edges:
    if (u,v) not in seen:
        uniq.append((u,v)); seen.add((u,v))
write_test("16.in", n, uniq)

# 17: hub-heavy graph: a few hubs (0,1,2) connected to most nodes
n=100
hubs=[0,1,2]
edges=[]
for h in hubs:
    for v in range(3,n):
        edges.append((h,v))
# link the hubs together too
for a,b in combinations(hubs,2):
    edges.append((a,b))
write_test("17.in", n, edges)

# 18: big path with 1000 nodes (tests linear behavior)
n=1000
edges=[(i,i+1) for i in range(n-1)]
write_test("18.in", n, edges)

# 19: big star (center 0, 999 leaves)
n=1000
edges=[(0,i) for i in range(1,n)]
write_test("19.in", n, edges)

# 20: deterministic "random-like" chunk (first 2500 pairs) on 500 nodes
n=500
edges=[]
count=2500
for i in range(n):
    for j in range(i+1,n):
        edges.append((i,j))
        if len(edges)>=count:
            break
    if len(edges)>=count:
        break
write_test("20.in", n, edges)

# 21: 2x250 grid (thin rectangle)
rows,cols=2,250
n=rows*cols
edges=[]
for r in range(rows):
    for c in range(cols):
        u = r*cols + c
        if c+1<cols:
            edges.append((u, r*cols + (c+1)))
        if r+1<rows:
            edges.append((u, (r+1)*cols + c))
write_test("21.in", n, edges)

# 22: two big cliques (0..99 and 100..199) with a single bridge (0-100)
n=200
edges=[]
for i in range(0,100):
    for j in range(i+1,100):
        edges.append((i,j))
for i in range(100,200):
    for j in range(i+1,200):
        edges.append((i,j))
edges.append((0,100))
write_test("22.in", n, edges)

# 23: empty graph (no edges) — useful to check edge cases
write_test("23.in", 30, [])

# 24: moderate deterministic chunk (1000 pairs) on 250 nodes
n=250
edges=[]
count=1000
for i in range(n):
    for j in range(i+1,n):
        edges.append((i,j))
        if len(edges)>=count:
            break
    if len(edges)>=count:
        break
write_test("24.in", n, edges)

# 25: small forest built from four long paths plus a couple of cross links
n=200
edges=[]
start=0
for k in range(4):
    length = 50
    for i in range(start, start+length-1):
        edges.append((i, i+1))
    start += length
edges.append((0,50))
edges.append((100,150))
write_test("25.in", n, edges)
# 26: Windmill / friendship graph F_k (k triangles share the same hub 0)
# Edge order is "bad" for greedy matching: leaf-leaf edges come first.
# This can push matching-based 2-approx close to its worst-case.
k = 60
n = 1 + 2*k
edges = []
for i in range(k):
    a = 1 + 2*i
    b = a + 1
    edges.append((a, b))          # leaf-leaf first (hurts greedy)
for i in range(k):
    a = 1 + 2*i
    b = a + 1
    edges.append((0, a))
    edges.append((0, b))
write_test("26.in", n, edges)

# 27: Same windmill F_k, but edge order is "good":
# hub edges first (greedy tends to include hub early).
k = 60
n = 1 + 2*k
edges = []
for i in range(k):
    a = 1 + 2*i
    b = a + 1
    edges.append((0, a))          # hub edges first (helps greedy)
    edges.append((0, b))
for i in range(k):
    a = 1 + 2*i
    b = a + 1
    edges.append((a, b))
write_test("27.in", n, edges)

# 28: Circulant graph C(n; 1,2,7) on 30 nodes (6-regular, many triangles)
# Good "non-bipartite / many cycles" stress test, still small enough for exact.
n = 30
jumps = [1, 2, 7]
E = set()
for i in range(n):
    for d in jumps:
        u = i
        v = (i + d) % n
        a, b = sorted((u, v))
        E.add((a, b))
edges = sorted(E)
write_test("28.in", n, edges)

# 29: Lollipop graph: clique K20 plus a path of 20 nodes attached to node 0
# Hybrid structure: very dense part + sparse tail.
clique = 20
tail = 20
n = clique + tail
edges = []
# clique on 0..19
for i in range(0, clique):
    for j in range(i+1, clique):
        edges.append((i, j))
# tail path on 20..39
for i in range(clique, clique + tail - 1):
    edges.append((i, i+1))
# connect tail to clique node 0
edges.append((0, clique))
write_test("29.in", n, edges)

# 30: "Almost bipartite": K12,12 plus an odd cycle inside the left part
# Breaks bipartite exactness, but still structured.
L = 12
R = 12
n = L + R
E = set()
# K12,12 edges: left 0..11, right 12..23
for i in range(L):
    for j in range(L, L+R):
        E.add((i, j))
# add an odd cycle on left part: 0-1-2-3-4-0 (5-cycle)
cycle = [0, 1, 2, 3, 4]
for i in range(len(cycle)):
    u = cycle[i]
    v = cycle[(i+1) % len(cycle)]
    a, b = sorted((u, v))
    E.add((a, b))
edges = sorted(E)
write_test("30.in", n, edges)

# 31: Random-ish Erdos-Renyi G(n,p) (deterministic seed)
# Useful to see ratios > 1, and runtime scaling.
import random
rnd = random.Random(12345)
n = 120
p = 0.08
edges = []
for i in range(n):
    for j in range(i+1, n):
        if rnd.random() < p:
            edges.append((i,j))
write_test("31.in", n, edges)

# 32: Same n, higher density (harder)
rnd = random.Random(54321)
n = 120
p = 0.18
edges = []
for i in range(n):
    for j in range(i+1, n):
        if rnd.random() < p:
            edges.append((i,j))
write_test("32.in", n, edges)

# 33: Planted vertex cover (you know OPT exactly = k)
# Construct: choose cover set C size k, connect every other vertex to C, no edges among non-cover.
# Then OPT is exactly k (because every edge touches C and non-cover has no internal edges).
n = 200
k = 30
C = set(range(k))
edges = []
for u in range(k, n):
    # connect to 3 cover vertices deterministically
    edges.append((u, (u*3) % k))
    edges.append((u, (u*7) % k))
    edges.append((u, (u*11) % k))
# add a little density inside C so it's not trivial
for i in range(k):
    for j in range(i+1, k):
        if ((i*17 + j*31) % 7) == 0:
            edges.append((i,j))
# dedup
edges = sorted(set(tuple(sorted(e)) for e in edges))
write_test("33.in", n, edges)

# 34: Barbell graph: two cliques joined by a path (not a single bridge only)
# More complex than test 22; stresses solvers but not as extreme.
clique = 30
pathlen = 40
n = 2*clique + pathlen
edges = []
# left clique 0..29
for i in range(0, clique):
    for j in range(i+1, clique):
        edges.append((i,j))
# right clique (clique+pathlen)..(2*clique+pathlen-1)
offset = clique + pathlen
for i in range(offset, offset+clique):
    for j in range(i+1, offset+clique):
        edges.append((i,j))
# path in the middle: 30..(30+pathlen-1)
for i in range(clique, clique+pathlen-1):
    edges.append((i, i+1))
# connect left clique node 0 to path start, and path end to right clique node offset
edges.append((0, clique))
edges.append((clique+pathlen-1, offset))
write_test("34.in", n, edges)

# 35: 3-partite dense-ish (lots of triangles but not a clique)
# A,B,C sizes: 40 each, connect A-B and B-C and A-C with patterned sparsity
n = 120
A = range(0,40)
B = range(40,80)
C = range(80,120)
edges = []
for u in A:
    for v in B:
        if ((u*13 + v*7) % 5) != 0:
            edges.append((u,v))
for u in B:
    for v in C:
        if ((u*11 + v*17) % 6) != 0:
            edges.append((u,v))
for u in A:
    for v in C:
        if ((u*19 + v*23) % 7) != 0:
            edges.append((u,v))
write_test("35.in", n, edges)

# 36: "Hard-for-greedy" matching order on a bipartite graph (still bipartite)
# This isn't worst-case theoretically, but it tends to increase cover sizes for the greedy order.
L = 60
R = 60
n = L + R
edges = []
# build layered connections
for i in range(L):
    for d in (0,1,2,3):
        j = (i + d) % R
        edges.append((i, L+j))
# reorder edges to be "bad": group by right side first
edges.sort(key=lambda e: (e[1], e[0]))
write_test("36.in", n, edges)

# 37/38: same random-ish graph, different edge order
def make_gnp_edges(n, p, seed):
    rnd = random.Random(seed)
    edges=[]
    for i in range(n):
        for j in range(i+1,n):
            if rnd.random() < p:
                edges.append((i,j))
    return edges

n = 160
edges = make_gnp_edges(n, 0.10, 12345)
write_test("37.in", n, sorted(edges))              # lexicographic order
write_test("38.in", n, list(reversed(sorted(edges))))  # reversed order

# 39/40: same structured graph (dense core + sparse tail), different edge order
n = 200
core = 40
edges=[]
# clique core
for i in range(core):
    for j in range(i+1, core):
        edges.append((i,j))
# tail path
for i in range(core, n-1):
    edges.append((i, i+1))
# cross edges deterministic
for i in range(core, n, 3):
    edges.append((0, i))
    edges.append((1, i))

edges = [tuple(sorted(e)) for e in edges]
edges = sorted(set(edges))

write_test("39.in", n, edges)                      # "good" order
edges_shuffled = edges[:]
random.Random(777).shuffle(edges_shuffled)         # deterministic shuffle
write_test("40.in", n, edges_shuffled)             # "shuffled" order

print('Done. Tests written to `tests/`.')
