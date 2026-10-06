# Clean first
rm -f *.o *.mod tddvr_so

# Compile potential module (required)
gfortran -cpp -fopenmp -O3 -c potential_module.f90


# Compile potential module (required)
gfortran -cpp -fopenmp -O3 -c energy_calculation_module.f90


# Compile potential module (required)
gfortran -cpp -fopenmp -O3 -c adiabatic_module.f90

# Compile main program with disabled modules
gfortran -cpp -fopenmp -O3 -c SO-TDDVR.f90

# Link
gfortran -fopenmp -o tddvr_so potential_module.o SO-TDDVR.o energy_calculation_module.o adiabatic_module.o

# Run
./tddvr_so>b&
