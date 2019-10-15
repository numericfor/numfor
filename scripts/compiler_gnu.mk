# For the gnu-fortran compiler

# Try to guess the architecture	
FC:= gfortran
arch:=-march=native

FFLAGS:= -cpp -ffree-form $(arch) -std=f2008 \
	 -Wall -Werror -pedantic -Wno-trampolines -Wno-unused-function -Wno-maybe-uninitialized

# Option for a shared library
ifdef SHARED
  FFLAGS+=-fPIC
endif

ifdef	PROFILE			# Si queremos hacer profiling (generalmente no)
  FFLAGS_PROF:=-pg -Q  -fprofile-arcs -ftest-coverage
  runprofile:=gprof $(EXECUT) > $(out).dat && cat $(out).dat | gprof2dot | dot -Tpng -o $(out).png
else
  FFLAGS_PROF:= -fomit-frame-pointer -pipe
  runprofile:=
endif

ifdef DEBUG
  FFLAGS_PAR:= 
  FFLAGS+= -g -ffpe-trap=invalid -fbounds-check
else				# Super-optimizado (no sé cuánto afecta)
  FFLAGS_PAR:= 
  FFLAGS+= -O3 -funroll-all-loops
endif

FFLAGS_EXTRA= -J $(MODD) $(FFLAGS_PAR) $(FFLAGS_PROF) $(FFLAGS_USER)

