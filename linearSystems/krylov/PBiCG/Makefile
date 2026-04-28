# Compiler configuration
CXX ?= g++
CXXFLAGS ?= -std=c++11 -Wall -Wextra -O2
BUILD_DIR ?= build
TARGET ?= $(BUILD_DIR)/PBiCG
SOURCE ?= PBiCGexample.cpp PBiCG.cpp

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(SOURCE) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE)

run: $(TARGET)
	./$(TARGET)

debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)

clean:
	rm -rf $(BUILD_DIR)

help:
	@echo "Available commands:"
	@echo "  make"
	@echo "  make run"
	@echo "  make debug"
	@echo "  make clean"
	@echo "  make SOURCE='example.cpp PBiCGsolver.cpp' TARGET=build/PBiCGsolver"

.PHONY: all run debug clean help
