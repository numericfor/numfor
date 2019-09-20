# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
# FC=gfortran
# FFLAGS=-O3
include /home/fiol/trabajo/programas/numfor/scripts/compiler_gnu.mk

OBJ=$(patsubst %.f90,%.o,$(shell ls *.f90))
all: lib$(LIB).a
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f90,$(@F)) -o $@
lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
clean:
	rm *.o






