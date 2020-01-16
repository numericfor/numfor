# This is just to avoid writing the path for each file
THISDIR := random
FILES := mt19937_64.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)
