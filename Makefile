# Makefile for Grid-based FMM for Superconductors
# Supports gfortran and ifort

# Compiler selection
FC = gfortran
# FC = ifort

# Compiler flags
ifeq ($(FC),gfortran)
    FFLAGS = -O3 -march=native -ffast-math -funroll-loops -Wall
    FFLAGS_DEBUG = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace
    OPENMP = -fopenmp
else ifeq ($(FC),ifort)
    FFLAGS = -O3 -xHost -fast -unroll
    FFLAGS_DEBUG = -g -O0 -check all -traceback -warn all
    OPENMP = -qopenmp
endif

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Source files (order matters for dependencies)
SOURCES = $(SRC_DIR)/types_mod.f90 \
          $(SRC_DIR)/potential_interface_mod.f90 \
          $(SRC_DIR)/fmm_grid_optimized_mod.f90 \
          $(SRC_DIR)/md_mod.f90 \
          $(SRC_DIR)/main.f90

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SOURCES))

# Module files
MODULES = $(BUILD_DIR)/*.mod

# Executable
EXEC = $(BIN_DIR)/fmm_2d

# Default target
all: directories $(EXEC)

# Create directories
directories:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Link executable
$(EXEC): $(OBJECTS)
	$(FC) $(FFLAGS) $(OPENMP) -o $@ $^
	@echo "Build complete: $(EXEC)"

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(OPENMP) -J$(BUILD_DIR) -c $< -o $@

# Dependencies (module dependencies)
$(BUILD_DIR)/potential_interface_mod.o: $(BUILD_DIR)/types_mod.o
$(BUILD_DIR)/fmm_grid_optimized_mod.o: $(BUILD_DIR)/types_mod.o
$(BUILD_DIR)/md_mod.o: $(BUILD_DIR)/types_mod.o
$(BUILD_DIR)/main.o: $(BUILD_DIR)/types_mod.o \
                     $(BUILD_DIR)/potential_interface_mod.o \
                     $(BUILD_DIR)/fmm_grid_optimized_mod.o \
                     $(BUILD_DIR)/md_mod.o

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean all

# Run the program
run: all
	./$(EXEC)

# Test different expansion orders
test_orders: all
	@echo "Testing different expansion orders..."
	@for p in 2 4 6 8 10; do \
		echo ""; \
		echo "========================================"; \
		echo "Testing with p = $$p"; \
		echo "========================================"; \
		sed -i "s/p_order = [0-9]*/p_order = $$p/" $(SRC_DIR)/main.f90; \
		$(MAKE) clean all; \
		./$(EXEC); \
	done

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) *.mod *.o frame_*.csv

# Clean everything including output
distclean: clean
	rm -f frame_*.csv

# Help
help:
	@echo "Available targets:"
	@echo "  all          - Build the executable (default)"
	@echo "  debug        - Build with debug flags"
	@echo "  run          - Build and run the program"
	@echo "  test_orders  - Test with different expansion orders"
	@echo "  clean        - Remove build artifacts"
	@echo "  distclean    - Remove all generated files"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "To use Intel compiler: make FC=ifort"

.PHONY: all directories debug run test_orders clean distclean help
