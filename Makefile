# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17
RELEASE_FLAGS := -O3

# Directories
OBJ_DIR := obj
BIN_DIR := bin
DATA_DIR := output
SRC_DIR := src
INLINE_DIR := $(SRC_DIR)/inline

# Source files
SOURCES := $(wildcard $(SRC_DIR)/*.cc)
OBJECTS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SOURCES))
INLINE_HEADERS := $(wildcard $(INLINE_DIR)/*.hh)
EXECUTABLE := $(BIN_DIR)/main

# Colors for pretty output
RED := \033[0;31m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC := \033[0m

.PHONY: all release clean clean_all help

all: release

release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(EXECUTABLE)
	@echo "$(GREEN)Build completed in release mode$(NC)"

$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "$(YELLOW)Linking $(EXECUTABLE)$(NC)"
	@$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INLINE_HEADERS) | $(OBJ_DIR)
	@echo "$(YELLOW)Compiling $<$(NC)"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR) $(DATA_DIR):
	@mkdir -p $@
	@echo "$(GREEN)Created directory: $@$(NC)"

clean:
	@rm -f $(OBJECTS) $(EXECUTABLE)
	@echo "$(RED)Removed object files and executable$(NC)"

clean_all: clean
	@rm -rf $(DATA_DIR)/*.dat $(DATA_DIR)/gif/*
	@echo "$(RED)Removed all generated data and logs$(NC)"

help:
	@echo "Available targets:"
	@echo "  all/release - Build in release mode (default)"
	@echo "  clean       - Remove object files and executable"
	@echo "  clean_all   - Remove all generated files (including data and logs)"
	@echo "  help        - Show this help message"
