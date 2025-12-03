# Makefile for Fortran Multipole Simulation
# ===========================================

# Compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall -Wextra -fbacktrace -g
LDFLAGS =

# Python
PYTHON = python3

# Targets
TARGET = test_simulation
MODULES = legendre_data.o multipole_module.o
OBJECTS = test_simulation.o

# Source files
SRC_DATA = legendre_data.f90
SRC_MODULE = multipole_module.f90
SRC_MAIN = test_simulation.f90

# Default target
all: generate_data $(TARGET)

# Generate Legendre data module
generate_data: $(SRC_DATA)

$(SRC_DATA): generate_coefficients.py
	@echo "Generating Legendre polynomial coefficients..."
	$(PYTHON) generate_coefficients.py 20 360
	@echo "Done!"

# Build executable
$(TARGET): $(MODULES) $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete!"

# Compile modules
legendre_data.o: $(SRC_DATA)
	@echo "Compiling legendre_data module..."
	$(FC) $(FFLAGS) -c $<

multipole_module.o: $(SRC_MODULE) legendre_data.o
	@echo "Compiling multipole_module..."
	$(FC) $(FFLAGS) -c $<

# Compile main program
test_simulation.o: $(SRC_MAIN) multipole_module.o
	@echo "Compiling test_simulation..."
	$(FC) $(FFLAGS) -c $<

# Run simulation
run: $(TARGET)
	@echo "Running simulation..."
	./$(TARGET)
	@echo ""
	@echo "Simulation complete! Visualize with:"
	@echo "  python3 visualize.py"

# Visualize results
visualize: run
	@echo "Creating visualization..."
	$(PYTHON) visualize.py

# Plot energy
plot_energy:
	$(PYTHON) visualize.py --energy

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -f *.o *.mod $(TARGET)
	rm -f frame_*.csv energy.log
	rm -f energy_evolution.png animation.mp4
	@echo "Clean complete!"

# Clean data files only
clean_data:
	@echo "Cleaning data files..."
	rm -f frame_*.csv energy.log
	rm -f energy_evolution.png animation.mp4

# Clean everything including generated coefficients
distclean: clean
	@echo "Removing generated files..."
	rm -f $(SRC_DATA)
	@echo "Distclean complete!"

# Help
help:
	@echo "Fortran Multipole Simulation - Makefile"
	@echo "========================================"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build everything"
	@echo "  make run          - Build and run simulation"
	@echo "  make visualize    - Run and visualize results"
	@echo "  make plot_energy  - Plot energy evolution"
	@echo "  make clean        - Remove build and output files"
	@echo "  make clean_data   - Remove only simulation output"
	@echo "  make distclean    - Remove everything including generated code"
	@echo "  make help         - Show this help"
	@echo ""
	@echo "Generated files:"
	@echo "  legendre_data.f90 - Precomputed coefficients"
	@echo "  frame_*.csv       - Simulation frames"
	@echo "  energy.log        - Energy evolution log"
	@echo ""

.PHONY: all generate_data run visualize plot_energy clean clean_data distclean help
