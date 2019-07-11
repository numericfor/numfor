# Makefile, by J. FIOL 

SHELL = /bin/bash


################################################################################
######                Information on the project
#               
PRJ = "numfor"
VERSION = "(Version $(shell git rev-parse --short --verify HEAD))"

########	Variables defining Directory Structure	########################
#	Path where make is invoked	
# top_dir:=$(realpath $(CURDIR))# Root of the project
#	Absolute path of this Makefile
top_dir:=$(realpath $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))
# #	Relative path of this Makefile
# top_dir:=$(dir $(firstword $(MAKEFILE_LIST)))


# Main Path of code (SRCD), tests (TSTD)
# object files *.obj, modules (*.mod), executable files (in BIND), and documentation
SRCD:=$(top_dir)/src
TSTD:=$(SRCD)/tests
PRGD:=$(SRCD)/progs

OBJD:=$(top_dir)/lib
MODD:=$(top_dir)/finclude
BIND:=$(top_dir)/bin
DOCDIR:=$(top_dir)/doc

################################################################################
#	PROGRAMS
depends=$(top_dir)/scripts/sfmakedepend

######################################################################
# Path to *.mod and source files
# INCLUDES:= -I $(MODD) -I $(SRCD)
INCLUDES:= -I $(MODD)

# Path to *.o files
LDFLAGS:=-L $(OBJD)

# # ######################################################################
MODULES := utils

# # ######################################################################
FMODULES =  $(addprefix $(SRCD)/,$(MODULES))
# look for include files in
# each of the modules
INCLUDES += $(patsubst %,-I %, $(FMODULES))

# extra libraries if required
LIBS :=
# each module will add to this
SRC :=

# include the description for each module
# include $(patsubst %,%/module.mk,$(FMODULES))
include $(addsuffix /module.mk,$(FMODULES))

# determine the object files
# OBJ := $(patsubst %.f90,%.o,$(filter %.f90,$(SRC)))
OBJ := $(addprefix $(OBJD)/,$(notdir $(SRC:.f90=.o)))

# VPATH:= $(patsubs %,:%,)
# dependencies
include $(SRC:.f90=.d)
# calculate code dependencies
%.d: %.f90
	$(depends) $(INCLUDES)  $< > $@

tt:
	@echo $(INCLUDES)
	@echo 'SRC: ' $(SRC)
	@echo 'OBJ: ' $(OBJ)


# # # ######################################################################
# # Source files	
# LIB_NAMES := utils/basic.f90 utils/strings.f90 utils/oopstring.f90
# XTR_NAMES := 

# XTR_SOURCES := $(addprefix $(SRCD)/,$(XTR_NAMES))
# LIB_SOURCES := $(addprefix $(SRCD)/,$(LIB_NAMES))
# # ######################################################################	
# #	
# SOURCES:=  $(LIB_SOURCES) $(XTR_SOURCES)

# # Objects files	
# LIB_OBJECTS := $(addprefix $(OBJD)/,$(notdir $(LIB_NAMES:f90=o)))
# XTR_OBJECTS := $(addprefix $(OBJD)/,$(XTR_NAMES:f90=o))
# ALL_OBJECTS :=  $(LIB_OBJECTS) $(XTR_OBJECTS)



################################################################################

######################################################################
######                   TEMPORARY VARIABLES
DEBUG=2
######################################################################

######################################################################
############## Compiler-dependent commands and options ###############
# Try to guess the architecture	
arch:=-march=native

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


# Commands for compiling and linking programs	
FC = $(F95) $(INCLUDES) $(FFLAGS) $(FFLAGS_EXTRA)
FLINK = $(F95) $(LDFLAGS_EXTRA) $(LDFLAGS)
# 

# ########################################################################
# ##########################  IMPLICIT RULES  ###########################

%.o : %.f90
	$(FC) -c -o $@ $<
# $(OBJD)/%.o : $(SRCD)/utils/%.f90
# 	$(FC) -c -o $@ $<

# $(OBJD)/%.o : $(SRCD)/%.F90
# 	$(FC) -c -o $@ $<

# $(OBJD)/%.o : $(TSTD)/%.f90
# 	$(FC) -c -o $@ $<

# $(OBJD)/%.o : $(PRGD)/%.F90
# 	$(FC) -c -o $@ $<

# $(TSTD)%.o: $(TSTD)%.f90
# 	$(FC) -c -o $@ $<

%.f90: %.F90
	$(FC) -E $< > $@  


# ########################################################################
# ##########################  EXPLICIT RULES  ###########################

all: library

library: $(OBJD)/strings.o $(MODD)/strings.mod

# Once the library is working, we write tests, and compile them with this rule	
$(test): $(LIB_OBJECTS) $(OBJD)/$(test).o $(TSTD)/$(test).f90
	$(FLINK) $^ -o $(BIND)/$(test)

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
oopstring.o: strings.mod
fstring.mod: $(top_dir)/src/utils/strings/oopstring.o
strings.mod: $(top_dir)/src/utils/strings/strings.o
