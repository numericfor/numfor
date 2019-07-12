# Makefile, by J. FIOL 

SHELL = /bin/bash


################################################################################
######                Information on the project
#               
PRJ = "numfor"
VERSION = "(Version $(shell git rev-parse --short --verify HEAD))"

########	Definition of directory Structure	########################
#	Absolute path of this Makefile
top_dir:=$(realpath $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))

# Main Path of code (SRCD), tests (TSTD)
# object files *.obj (OBJD), modules *.mod (MODD)
# executable files (BIND), and documentation (DOCDIR)
SRCD:=$(top_dir)/src
TSTD:=$(SRCD)/tests
# PRGD:=$(SRCD)/progs

OBJD:=$(top_dir)/lib
MODD:=$(top_dir)/finclude
BIND:=$(top_dir)/bin
DOCDIR:=$(top_dir)/doc


######################################################################
##############		     Some options               ##############

# Choose the compiler. Default is gnu
compiler = gnu

################################################################################
#	PROGRAMS
depends:=$(top_dir)/scripts/sfmakedepend
AR:=ar
RM:=rm
######################################################################

######################################################################
# Path to *.mod and source files
INCLUDES:= -I $(MODD)

# Path to *.o files
LDFLAGS:=-L $(OBJD)

# # ######################################################################
MODULES := utils

# # ######################################################################
FMODULES =  $(addprefix $(SRCD)/,$(MODULES))
# look for include files in each of the modules
INCLUDES += $(patsubst %,-I %, $(FMODULES))

# extra libraries if required
LIBS := -lnumfor
# each module will add to this
SRC :=

# include the description for each module
include $(addsuffix /module.mk,$(FMODULES))

# The object files will be in the same places that the code
OBJ:= $(SRC:.f90=.o)
deps:= $(SRC:.f90=.d)

# dependencies
include $(deps)


# ####################################################################
# ############## Compiler-dependent commands and options #############
ifeq ($(compiler),intel)
  include compiler_intel.mk
else
  include compiler_gnu.mk
endif
# ####################################################################
# ############# Commands for compiling and linking programs ##########
FC = $(F95) $(INCLUDES) $(FFLAGS) $(FFLAGS_EXTRA)
FLINK = $(F95) $(LDFLAGS_EXTRA) $(LDFLAGS)


# ####################################################################
# ##########################  IMPLICIT RULES  ########################
%.o : %.f90
	$(FC) -c -o $@ $<

# calculate code dependencies
%.d: %.f90
	$(depends) $(INCLUDES)  $< > $@


# ####################################################################
# ##########################  EXPLICIT RULES  ########################

library: $(OBJD)/libnumfor.a

$(OBJD)/libnumfor.a: $(OBJ) 
	$(AR) rcs $@ $<

# Once the library is working, we write tests, and compile them with this rule	
$(test): $(TSTD)/$(test).f90
	$(FC) $(LDFLAGS)  $^ -o $(BIND)/$(test) $(LIBS)

$(prog): $(LIB_OBJECTS) $(OBJD)/$(prog).o
	$(FLINK) $^ -o $(BIND)/$@


# ########################################################################
# ######## DOCUMENTACIÓN

doc:  $(DOCDIR)/html

$(DOCDIR)/html: Doxyfile $(SRC) README.md
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |' Doxyfile	|\
	sed -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 ${VERSION}|' | \
	sed -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${DOCDIR}|' > $(PRJ).dox
	doxygen $(PRJ).dox

.PHONY: library tags TAGS clean clean-all view-doc doc-api doc

doc-api: Doxyfile $(SRC) 
	sed -e 's|\(INPUT[ ]*=\)\(.*\)|\1 ${SUBDIRS} |' Doxyfile	|\
	sed -e 's|\(PROJECT_NUMBER[ ]*=\)\(.*\)|\1 ${VERSION}|' | \
	sed -e 's|\(STRIP_FROM_PATH[ ]*=\)\(.*\)|\1 ${DOCDIR}|'	|\
	sed -e 's|\(OUTPUT_DIRECTORY[ ]*=\)\(.*\)|\1 ../doc-api|' | sed -e 's|\(EXTRACT_PRIVATE[ ]*=\)\(.*\)|\1 NO|'	> $(PRJ)_api.dox
	doxygen $(PRJ)_api.dox

view-doc: $(DOCDIR)/html
	gio open $(DOCDIR)/html/index.html || firefox $(DOCDIR)/html/index.html

tags: TAGS

TAGS: Makefile $(SRC)
	etags $(SRC)

clean:
	$(RM) -f $(top_dir)/*~
	$(RM) -f $(SRCD)/*~
	$(RM) -f $(SRCD)/*/*~
	$(RM) -f $(TSTD)/*~
	$(RM) -f $(top_dir)/docs/sources/*~


clean-doc:
	$(RM) -f $(DOCDIR)/html

clean-all: clean clean-doc
	$(RM) $(MODD)/*.mod
	$(RM) $(OBJ)
	$(RM) $(deps)
