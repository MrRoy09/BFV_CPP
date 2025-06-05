CXX = clang++
CXXFLAGS = -O2

# Directories
INCLUDE_DIR = ./include
TEST_DIR = ./test

# Find all source files
INCLUDE_SOURCES = $(wildcard $(INCLUDE_DIR)/*.cpp)
BASE_SOURCES = $(wildcard *.cpp)

# Find all test files and extract test names
TEST_FILES = $(wildcard $(TEST_DIR)/*.cpp)
TEST_NAMES = $(basename $(notdir $(TEST_FILES)))
TEST_TARGETS = $(TEST_NAMES)

# Default target
.PHONY: all clean

all: $(TEST_TARGETS)

# Pattern rule to build each test
$(TEST_TARGETS): %: $(TEST_DIR)/%.cpp $(INCLUDE_SOURCES) $(BASE_SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDE_SOURCES) $(BASE_SOURCES) $(TEST_DIR)/$*.cpp -o $@

# Alternative single test build (backwards compatibility)
test: $(TEST)
ifdef TEST
$(TEST): $(TEST_DIR)/$(TEST).cpp $(INCLUDE_SOURCES) $(BASE_SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDE_SOURCES) $(BASE_SOURCES) $(TEST_DIR)/$(TEST).cpp -o $(TEST)
endif

clean:
	rm -f $(TEST_TARGETS)

# Debug target to show what will be built
debug:
	@echo "Include sources: $(INCLUDE_SOURCES)"
	@echo "Base sources: $(BASE_SOURCES)"  
	@echo "Test files: $(TEST_FILES)"
	@echo "Test targets: $(TEST_TARGETS)"

# Help target
help:
	@echo "Usage:"
	@echo "  make all          - Build all tests"
	@echo "  make <test_name>  - Build specific test"
	@echo "  make TEST=<name>  - Build specific test (old way)"
	@echo "  make clean        - Remove all test executables"
	@echo "  make debug        - Show build information"
	@echo ""
	@echo "Available tests: $(TEST_NAMES)"