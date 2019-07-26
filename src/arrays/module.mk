# This is just to avoid writing the path for each file
THISDIR := arrays
FILES := grids.f90 histogram.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
