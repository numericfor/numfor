# This is just to avoid writing the path for each file
THISDIR := arrays
FILES := array_utils.f90 grids.f90 sorting.f90 histogram.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
