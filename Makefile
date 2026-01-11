CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -pedantic
TARGET   := vc
SRC      := main.cpp

# default timeout for exact algs in ms (can override: make bench TIMEOUT=5000)
TIMEOUT  ?= 2000

# ---- optional GLPK (LP real) ----
GLPK_CXXFLAGS := -DUSE_GLPK
GLPK_LDLIBS   := -lglpk

.PHONY: all clean bench bench_glpk bb fpt match lp help

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Build with GLPK enabled: make vc_glpk
vc_glpk: $(SRC)
	$(CXX) $(CXXFLAGS) $(GLPK_CXXFLAGS) $(SRC) -o $(TARGET) $(GLPK_LDLIBS)

# Benchmark pe toate testele (fara GLPK): make bench (scrie results.csv)
bench: $(TARGET)
	./$(TARGET) --bench tests --out results.csv --timeout_ms $(TIMEOUT)

# Benchmark pe toate testele (cu GLPK): make bench_glpk
bench_glpk: vc_glpk
	./$(TARGET) --bench tests --out results.csv --timeout_ms $(TIMEOUT)

# Run single test file: make bb T=tests/22.in
T ?= tests/01.in

bb: $(TARGET)
	./$(TARGET) --alg=bb --timeout_ms $(TIMEOUT) < $(T)

fpt: $(TARGET)
	./$(TARGET) --alg=fpt --timeout_ms $(TIMEOUT) < $(T)

match: $(TARGET)
	./$(TARGET) --alg=match --timeout_ms $(TIMEOUT) < $(T)

lp: $(TARGET)
	./$(TARGET) --alg=lp --timeout_ms $(TIMEOUT) < $(T)

clean:
	rm -f $(TARGET) results.csv

help:
	@echo "Targets:"
	@echo "  make / make all           - build vc (no GLPK; lp fallback)"
	@echo "  make vc_glpk              - build vc with GLPK (real LP)"
	@echo "  make bench                - bench all tests (no GLPK)"
	@echo "  make bench_glpk           - bench all tests (with GLPK)"
	@echo "  make bb|fpt|match|lp T=tests/22.in TIMEOUT=2000"
