# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>
# Pearu
# Modified by J. Fiol	

%.o : %.f90
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f90,$(@F)) -o $@

SRC:=$(shell ls *.f90)

OBJ=$(patsubst %.f90,%.o,$(shell ls *.f90))
all: lib$(LIB).a

# Here are defined the compiler-dependent variables. Change later for relative path	
include /home/fiol/trabajo/programas/numfor/scripts/compiler_gnu.mk

lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
	cp lib$(LIB).a /home/fiol/trabajo/programas/numfor/lib
	# cp *.o ..		

clean:
	rm *.o








