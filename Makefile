# Makefile, by J. FIOL 

SHELL = /bin/bash

################################################################################
######## TEMPORARY VARIABLES
DEBUG=2
################################################################################

#	Directorios	
# top_dir:=$(realpath $(dir $(PWD)))# Raíz de nuestro proyecto 
top_dir:=$(realpath $(PWD))# Raíz de nuestro proyecto

# Directorios para código principal (SRCD), tests (TSTD)
# objetos *.obj, modulos (*.mod), ejecutables (en BIND)
SRCD:=$(top_dir)/src
TSTD:=$(SRCD)/tests
PRGD:=$(SRCD)/progs

OBJD:=$(top_dir)/lib
MODD:=$(top_dir)/finclude
BIND:=$(top_dir)/bin
# Directorios para la documentación (código y teoría)
DOCDIR:=$(top_dir)/doc

PRJ = "numfor"
VERSION = "(Version $(shell git rev-parse --short --verify HEAD))"

# ######################################################################
LIB_NAMES := utils/strings/strings.f90
XTR_NAMES = 


# Por ahora nada acá		
XTR_SOURCES := $(addprefix $(SRCD)/,$(XTR_NAMES))
LIB_SOURCES := $(addprefix $(SRCD)/,$(LIB_NAMES))
# ######################################################################	
#	
SOURCES:=  $(LIB_SOURCES)

# Los objetos	
LIB_OBJECTS = $(addprefix $(OBJD)/,$(notdir $(LIB_NAMES:f90=o)))
XTR_OBJECTS = $(addprefix $(OBJD)/,$(XTR_NAMES:f90=o))
ALL_OBJECTS =  $(XTR_OBJECTS) $(LIB_OBJECTS)


######################################################################
############## Compiler-dependent commands and options ###############
# Try to guess the architecture	
arch:=-march=native
# arch:=-march=core2
ifeq ($(compiler),intel)
# For the intel compiler
  F95:=ifort
  FFLAGS:= $(arch) -check all -warn all -std03 -fpp -warn errors -warn stderrors -fpe3
  ifdef DEBUG
    FFLAGS= -g -ftrapuv -traceback -ip -pad
  else
    FFLAGS= -O3  -ip -pad
  endif
  FFLAGS_EXTRA= -module $(MODD)  $(FFLAGS_PAR) $(FFLAGS_PROF) $(FFLAGS_USER)
else
# For the gnu-fortran compiler
  F95:= gfortran
  ifdef	PROFILE			# Si queremos hacer profiling (generalmente no)
    FFLAGS_PROF=-pg -Q  -fprofile-arcs -ftest-coverage
    runprofile=gprof $(EXECUT) > $(out).dat && cat $(out).dat | gprof2dot | dot -Tpng -o $(out).png
  else
    FFLAGS_PROF= -fomit-frame-pointer -pipe
    runprofile=
  endif
  FFLAGS:= -ffree-form -ffree-line-length-0 $(arch) -std=f2008 -fPIC -Wall -Werror -pedantic -Wno-trampolines -Wno-unused-function -Wno-maybe-uninitialized

  ifdef DEBUG
    FFLAGS_PAR:= 
    FFLAGS+= -g -ffpe-trap=invalid -fbounds-check
  else				# Super-optimizado (no sé cuánto afecta)
    FFLAGS_PAR:= 
    FFLAGS+= -O3 -funroll-all-loops
  endif
  FFLAGS_EXTRA= -J $(MODD) $(FFLAGS_PAR) $(FFLAGS_PROF) $(FFLAGS_USER)
endif
######################################################################
# Donde buscar *.o y *.mod
# INCLUDES= -I $(MODD) -I $(OBJD)
INCLUDES= -I $(MODD) -I $(SRCD)
# Donde encontrar los *.o no debería ser necesario acá		
LDFLAGS=-L $(OBJD)


FC = $(F95) $(INCLUDES) $(FFLAGS) $(FFLAGS_EXTRA)
# FLINK = $(FC) $(LDFLAGS_EXTRA) $(LDFLAGS)
FLINK = $(F95) $(LDFLAGS)
# 

# ########################################################################
# ##########################  IMPLICIT RULES  ###########################

$(OBJD)/%.o : $(SRCD)/utils/strings/%.f90
	$(FC) -c -o $@ $<

$(OBJD)/%.o : $(SRCD)/%.F90
	$(FC) -c -o $@ $<

$(OBJD)/%.o : $(TSTD)/%.f90
	$(FC) -c -o $@ $<

$(OBJD)/%.o : $(PRGD)/%.F90
	$(FC) -c -o $@ $<

$(TSTD)%.o: $(TSTD)%.f90
	$(FC) -c -o $@ $<

%.f90: %.F90
	$(FC) -E $< > $@  


# ########################################################################
# ##########################  EXPLICIT RULES  ###########################

all: library

library: $(OBJD)/strings.o $(MODD)/strings.mod

#  Esta regla es para todos los tests (excepto test_FForma). Por ejemplo:
#     make -k  exe=test_dubois11 test_dubois11
$(test): $(LIB_OBJECTS) $(TSTD)/test_$(test).o
	$(FLINK) $^ -o $(BIND)/test_$(test)

$(prog): $(LIB_OBJECTS) $(OBJD)/$(prog).o
	$(FLINK) $^ -o $(BIND)/$@


# ########################################################################
# ######## DOCUMENTACIÓN

doc:  $(DOCDIR)/html

$(DOCDIR)/html: Doxyfile $(SOURCES) README.md
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |' Doxyfile	|\
	sed -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 ${VERSION}|' | \
	sed -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${top_dir}|' > $(PRJ).dox
	doxygen $(PRJ).dox

doc-api: Doxyfile $(SOURCES) 
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |' Doxyfile	|\
	sed -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 ${VERSION}|' | \
	sed -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${top_dir}|'	|\
	sed -e 's|\(OUTPUT_DIRECTORY[ ]*=\)\(.*\)|\1 ../doc-api|' | sed -e 's|\(EXTRACT_PRIVATE[ ]*=\)\(.*\)|\1 NO|'	> $(PRJ)_api.dox
	doxygen $(PRJ)_api.dox

view-doc: $(DOCDIR)/html
	gio open $(DOCDIR)/html/index.html || firefox $(DOCDIR)/html/index.html

tags: TAGS

TAGS: Makefile $(SOURCES)
	etags $(SOURCES)

.PHONY: All clean clean-all view-doc doc-api doc

clean:
	rm -f $(SRCD)/*~

clean-doc:
	rm -f doc/html

clean-all: clean
	rm -f $(ALL_OBJECTS) $(OBJD)/*.mod




myobjects: $(SOURCES)
	sfmakedepend -f Makefile $(INCLUDES) $^ $(MNFILE) 
	sed -e 's|$(top_dir)|$$\(top_dir\)|g' Makefile | sed -e 's|[ ]iso_c_binding.mod| |g' > newMakefile
	mv newMakefile Makefile
	mv Makefile.old Makefile.bak


# DO NOT DELETE THIS LINE - used by make depend
