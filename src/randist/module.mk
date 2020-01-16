# This is just to avoid writing the path for each file
THISDIR := randist
FILES := uniform.f90 gauss.f90 exponential.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
