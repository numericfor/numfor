# Values for the intel compiler
FC:=ifort
FFLAGS:= $(arch) -check all -warn all -std03 -fpp -warn errors -warn stderrors -fpe3
ifdef DEBUG
  FFLAGS= -g -ftrapuv -traceback -ip -pad
else
  FFLAGS= -O3  -ip -pad
endif
FFLAGS_EXTRA= -module $(MODD)  $(FFLAGS_PAR) $(FFLAGS_PROF) $(FFLAGS_USER)
