# Makefile, by J. FIOL 

SHELL = /bin/bash


################################################################################
######                Information on the project
# 
PRJ:= numfor
# Version and date of current revision
VERSION:=$(shell git show -s --format="%h")
DATE:=$(shell git show -s --format="%ci"| cut -d " " -f 1)

########	Definition of directory Structure	########################
#	Absolute path of this Makefile
top_dir:=$(realpath $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))

# Main Path of code (SRCD), tests (TSTD)
# object files *.obj (OBJD), modules *.mod (MODD)
# executable files (BIND), and documentation (DOCDIR)
SRCD:=$(top_dir)/src
TSTD:=$(SRCD)/tests

SCRPTD:=$(top_dir)/scripts
OBJD:=$(top_dir)/lib
MODD:=$(top_dir)/finclude
BIND:=$(top_dir)/bin
# Extra documentation sources and examples
DOCSD:=$(top_dir)/doc/sources
EXMPD:=$(top_dir)/doc/examples
# Directory where documentation will be built
DOCDIR:=$(top_dir)/docs


######################################################################
##############		     Some options               ##############

# Choose the compiler. Default is "gnu". Accepted also "intel" (untested)
compiler = gnu
# DEBUG=YES

################################################################################
#	PROGRAMS
depends:=$(SCRPTD)/sfmakedepend
AR:=ar
RM:=rm -f
CP:=cp -f
LN:=ln -s -f
MKDIR:=mkdir -p
######################################################################

######################################################################
# Path to *.mod and source files
INCLUDES:= -I $(MODD)
# INCLUDES += -I $(SRCD)/interpolate/fitpack
# Path to *.o files
LDFLAGS:=-L $(OBJD)

# # ######################################################################
MODULES := utils arrays random randist interpolate integrate

# # ######################################################################
FMODULES =  $(addprefix $(SRCD)/,$(MODULES))
# look for include files in each of the modules
INCLUDES += $(patsubst %,-I %, $(FMODULES))


# extra libraries if required
LIBS := -lnumfor

# This is the main point-of-entry to the library.
# each module will add to this
SRC := 

# include the description for each module
include $(addsuffix /module.mk,$(FMODULES))

SRC += $(SRCD)/$(PRJ).f90

# The object files will be in the same places that the code
OBJ:= $(SRC:.f90=.o)
deps:= $(SRC:.f90=.d)

# ####################################################################
# ############## Compiler-dependent commands and options #############
ifeq ($(compiler),intel)
  include $(SCRPTD)/compiler_intel.mk
else
  include $(SCRPTD)/compiler_gnu.mk
endif
# ####################################################################
# ############# Commands for compiling and linking code     ##########
FLINK = $(FC) $(LDFLAGS_EXTRA) $(LDFLAGS)


# ####################################################################
# ##########################  IMPLICIT RULES  ########################
%.o : %.f90
	$(FC)  $(FFLAGS) $(FFLAGS_EXTRA) $(INCLUDES) -c -o $@ $<

# calculate code dependencies
%.d: %.f90
	$(depends) $(INCLUDES)  $< > $@

# ####################################################################
# ##########################  EXPLICIT RULES  ########################
numlib:= $(OBJD)/libnumfor.a

library: $(numlib)

$(numlib): $(OBJ) | $(OBJD)
	$(AR) rcs $@ $(OBJ) 

shared_library: $(OBJD)/libnumfor.so

$(OBJD)/libnumfor.so: $(OBJ)  | $(OBJD)
	$(FC) $(INCLUDES)  $(FFLAGS) $(FFLAGS_EXTRA) -shared -o $@ $^

$(OBJD) $(BIND):
	$(MKDIR) $@

# dependencies
include $(deps)


# ########################################################################
# Once the library is working, we write and compile tests with this rule	
# Rule used to create the examples. For instance, to make test_strings:
# make tst=oopstring test
# Default test will be test_strings
tst=strings
test: $(BIND)/test_$(tst)
$(BIND)/test_$(tst): $(TSTD)/test_$(tst).f90 $(OBJD)/libnumfor.a | $(BIND)
	$(FC)  $(FFLAGS) $(FFLAGS_EXTRA) $(LDFLAGS)  $< -o $@ $(LIBS)

# Once the library is working, we write and compile examples 
# Rule used to create the examples. For instance, to make ex_fstring:
# make ex=fstring1 example
# Default example will be ex_ftring1
ex=ftring1
example: $(BIND)/ex_$(ex) 
$(BIND)/ex_$(ex): $(EXMPD)/ex_$(ex).f90 $(OBJD)/libnumfor.a | $(BIND)
	$(FC)  $(FFLAGS) $(FFLAGS_EXTRA) $(LDFLAGS)  $< -o $@ $(LIBS)

# ########################################################################
# ######## INSTALLATION

# Defaul value. May be overriden by command-line
prefix:=$(HOME)/.local

include $(SCRPTD)/install.mk

# ########################################################################
# ######## DOCUMENTATION
DOXYF:=$(SCRPTD)/Doxyfile

doc:  $(DOCDIR)/index.html

$(DOCDIR)/index.html: Makefile $(DOXYF) $(SRC) $(top_dir)/README.md  $(DOCSD)/*.md $(EXMPD)/*.f90
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |' $(DOXYF) -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 "${VERSION} ($(DATE))"|'\
            -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${DOCDIR}|' > $(top_dir)/$(PRJ).dox
	cd $(top_dir) && doxygen $(PRJ).dox && cd -


.PHONY: library tags clean clean-backup clean-all view-doc doc-api doc

doc-api: $(DOXYF) $(SRC)
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |'  $(DOXYF) |\
	sed -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 ${VERSION}|' | \
	sed -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${DOCDIR}|'	|\
	sed -e 's|\(OUTPUT_DIRECTORY[ ]*=\)\(.*\)|\1 ../doc-api|' |\
	sed -e 's|\(EXTRACT_PRIVATE[ ]*=\)\(.*\)|\1 NO|' > $(top_dir)/$(PRJ)_api.dox
	cd $(top_dir) && doxygen $(PRJ)_api.dox && cd -


view-doc: doc
	gio open $(DOCDIR)/index.html || firefox $(DOCDIR)/index.html

# ########################################################################
# ######## TAGS (for editing)

tags: TAGS

TAGS: Makefile $(SRC)
	etags $(SRC)

# ########################################################################
# ######## Cleaning

clean: clean-backup clean-obj

clean-backup:
	$(RM) $(top_dir)/*~
	$(RM) $(top_dir)/*/*~
	$(RM) $(top_dir)/*/*/*~

clean-obj:
	$(RM) $(MODD)/*.mod
	$(RM) $(OBJD)/*.*
	$(RM) $(OBJ)
	$(RM) $(deps)

clean-doc:
	$(RM) -r $(DOCDIR)/html
	$(RM) $(top_dir)/$(PRJ).dox $(top_dir)/$(PRJ)_api.dox

clean-all: clean-backup clean-obj clean-doc
