# This is just to avoid writing the path for each file
THISDIR := utils
FILES := basic.f90 grids.f90 strings.f90 oopstring.f90 histogram.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
