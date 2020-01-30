FC = gfortran 
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O3
#FCFLAGS = -g 
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package

PROGRAMS = SPH_2018 
MODS = 7_variable_types.o 6_SPH_time_vars_2018.o 5_SPH_global_vars_2018.o \
       4_SPH_material_vars_2018.o 3_SPH_material_2018.o \
       2_SPH_main_2018.o 1_SPH_2018.o
sph: $(MODS)
#
#%: %.o
	$(FC)  $(FCFLAGS) -o sph $(MODS) 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.csv.*	

veryclean: clean
	rm -f *~ $(PROGRAMS)
	

