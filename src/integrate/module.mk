# This is just to avoid writing the path for each file
THISDIR := integrate
FILES := func_integ.f90 qsimpson.f90 wquadpack.f90 tanhsinh.f90 qadaptive.f90

SRC += $(addprefix $(SRCD)/$(THISDIR)/,$(FILES)) $(addprefix $(SRCD)/$(THISDIR)/,$(THISDIR).f90)

