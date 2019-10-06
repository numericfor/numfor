# This is just to avoid writing the path for each file
THISDIR :=  interpolate
# FILES := polynomial.f90 csplines.f90  fitp.f90 wfitpack.f90
FILES := polynomial.f90 csplines.f90 wfitpack.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)

