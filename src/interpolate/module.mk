# This is just to avoid writing the path for each file
THISDIR := interpolate
FILES := splines.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
