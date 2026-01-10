CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -pedantic
TARGET   := vc
SRC      := main.cpp

.PHONY: all clean run bb fpt match pd bench

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Benchmark pe toate testele: make bench (scrie results.csv)
bench: $(TARGET)
	./$(TARGET) --bench tests --out results.csv

clean:
	rm -f $(TARGET) results.csv
