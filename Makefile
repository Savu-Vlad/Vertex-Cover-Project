CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -pedantic -DUSE_GLPK
TARGET   := vc
SRC      := main.cpp
LDLIBS   := -lglpk

.PHONY: all clean bench help

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

bench: $(TARGET)
	./$(TARGET) --bench tests --out results.csv --timeout_ms 2000 --reps 5

help:
	@echo "Build:"
	@echo "  make                 -> builds ./$(TARGET) with GLPK"
	@echo ""
	@echo "Run (stdin graph):"
	@echo "  ./$(TARGET) --alg=bb|fpt|match|lp < input.in"
	@echo ""
	@echo "Benchmark:"
	@echo "  make bench"
	@echo "  ./$(TARGET) --bench tests --out results.csv --timeout_ms 2000 --reps 5 [--no_cleanup]"
	@echo ""
	@echo "Args:"
	@echo "  --alg=bb|fpt|match|lp"
	@echo "  --bench <folder>     (runs all *.in)"
	@echo "  --out <file.csv>"
	@echo "  --timeout_ms <int>"
	@echo "  --reps <int>         (median by total time)"
	@echo "  --no_cleanup         (skip cleanup_cover)"
	@echo "  --help               (prints CLI help)"
	@echo ""
	@echo "Fedora deps:"
	@echo "  sudo dnf install -y glpk-devel"
