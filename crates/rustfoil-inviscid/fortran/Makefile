# Makefile for XFOIL inviscid debug test harnesses
#
# This Makefile builds test programs that link against Xfoil-instrumented
# to extract reference data for validating the Rust implementation.

# Paths
XFOIL_INST = ../../../Xfoil-instrumented
XFOIL_SRC = $(XFOIL_INST)/src
XFOIL_BIN = $(XFOIL_INST)/bin
TESTDATA = ../../../testdata

# Compiler settings
FC = gfortran
FFLAGS = -O2 -fdefault-real-8 -fdefault-double-8 -I$(XFOIL_SRC)

# XFOIL object files needed for inviscid solve
XFOIL_OBJS = $(XFOIL_BIN)/xpanel.o \
             $(XFOIL_BIN)/xoper.o \
             $(XFOIL_BIN)/xsolve.o \
             $(XFOIL_BIN)/xfoil_debug.o \
             $(XFOIL_BIN)/naca.o \
             $(XFOIL_BIN)/xutils.o \
             $(XFOIL_BIN)/xgeom.o \
             $(XFOIL_BIN)/spline.o \
             $(XFOIL_BIN)/userio.o

# Targets
all: test_inviscid

# First ensure Xfoil-instrumented is built
$(XFOIL_BIN)/%.o:
	@echo "Building Xfoil-instrumented..."
	$(MAKE) -C $(XFOIL_BIN)

test_inviscid: test_inviscid.f
	$(FC) $(FFLAGS) -o $@ $< $(XFOIL_OBJS) -lm

# Run targets to generate reference data
run_all: test_inviscid
	./test_inviscid > $(TESTDATA)/inviscid_ref.json

run_geometry: test_inviscid
	./test_inviscid geometry > $(TESTDATA)/geometry_ref.json

run_psilin: test_inviscid
	./test_inviscid psilin > $(TESTDATA)/psilin_ref.json

run_ggcalc: test_inviscid
	./test_inviscid ggcalc > $(TESTDATA)/ggcalc_ref.json

run_clcalc: test_inviscid
	./test_inviscid clcalc > $(TESTDATA)/clcalc_ref.json

clean:
	rm -f test_inviscid *.o

.PHONY: all clean run_all run_geometry run_psilin run_ggcalc run_clcalc
