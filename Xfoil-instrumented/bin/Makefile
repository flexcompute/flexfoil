#*********************************************************
# Makefile for XFOIL V6.93 - INSTRUMENTED VERSION
# Based on original XFOIL Makefile
# Modified for debug output instrumentation
#*********************************************************


SHELL = sh
BINDIR = .

PROGS = xfoil_instrumented

SRC = ../src
OSRC = ../osrc

# Add xfoil_debug.o for instrumentation
XFOILOBJ = xfoil.o xpanel.o xoper.o xtcam.o xgdes.o xqdes.o xmdes.o \
xsolve.o xbl.o xblsys.o xpol.o xplots.o pntops.o xgeom.o xutils.o modify.o \
blplot.o polplt.o aread.o naca.o spline.o plutil.o iopol.o gui.o sort.o \
dplot.o profil.o xfoil_debug.o

XUTILOBJ  = userio.o

FTNLIB = 

##----------------------------------------------------
OSOBJ = frplot0.o

# Use this for individual TS-wave frequency plotting
OSOBJ = frplot.o ntcalc.o osmap.o getosfile.o

##----------------------------------------------------
PLTOBJ = ../plotlib/libPlt.a 

# The extra location arg here is for Linux which places X libs in /usr/X11R6
PLTLIB = -L/usr/X11R6/lib -lX11

###================================================
###  Default compilers and flags
FC = f77
FFLAGS  = -O 
FFLOPT  = -O
INSTALLCMD = install -s

CC = cc
CFLAGS = -O -DUNDERSCORE

##--------------------------
### macOS with gfortran (Homebrew)
FC = gfortran
DBL = -fdefault-real-8
FIX = -std=legacy -fallow-argument-mismatch
FFLAGS = -O2 $(DBL) $(FIX)
FFLOPT = -O2 $(DBL) $(FIX)
PLTOBJ = ../plotlib/libPlt.a
PLTLIB = -L/opt/homebrew/opt/libx11/lib -lX11

CC = cc
CFLAGS = -O -DUNDERSCORE -I/opt/homebrew/opt/libx11/include -I/opt/homebrew/opt/xorgproto/include


all:	 $(PROGS)

install: 
	$(INSTALLCMD) $(PROGS) $(BINDIR)

clean:
	-/bin/rm $(PROGS)
	-/bin/rm $(XFOILOBJ) $(XUTILOBJ) $(OSOBJ)


xfoil_instrumented: $(XFOILOBJ) $(XUTILOBJ) $(OSOBJ)
	$(FC) -o xfoil_instrumented $(XFOILOBJ) $(XUTILOBJ) $(OSOBJ) $(PLTOBJ) $(PLTLIB) $(FTNLIB)


# Debug module (new) - needs include path for XBL.INC
xfoil_debug.o: $(SRC)/xfoil_debug.f $(SRC)/XBL.INC $(SRC)/BLPAR.INC
	$(FC) -c $(FFLAGS) -I$(SRC) $(SRC)/xfoil_debug.f

xfoil.o: $(SRC)/xfoil.f  $(SRC)/XFOIL.INC $(SRC)/XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xfoil.f
xpanel.o: $(SRC)/xpanel.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLOPT)  $(SRC)/xpanel.f
xoper.o: $(SRC)/xoper.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xoper.f
xsolve.o: $(SRC)/xsolve.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLOPT)  $(SRC)/xsolve.f
dplot.o: $(SRC)/dplot.f  $(SRC)/XFOIL.INC
	$(FC) -c $(FFLOPT)  $(SRC)/dplot.f
xtcam.o: $(SRC)/xtcam.f  $(SRC)/XFOIL.INC $(SRC)/XDES.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xtcam.f
xgdes.o: $(SRC)/xgdes.f  $(SRC)/XFOIL.INC $(SRC)/XDES.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xgdes.f
xqdes.o: $(SRC)/xqdes.f  $(SRC)/XFOIL.INC $(SRC)/XDES.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xqdes.f
xmdes.o: $(SRC)/xmdes.f  $(SRC)/XFOIL.INC $(SRC)/XDES.INC $(SRC)/CIRCLE.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xmdes.f
xbl.o:   $(SRC)/xbl.f    $(SRC)/XFOIL.INC $(SRC)/XBL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xbl.f
xblsys.o: $(SRC)/xblsys.f                  $(SRC)/XBL.INC $(SRC)/BLPAR.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xblsys.f
xplots.o: $(SRC)/xplots.f $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xplots.f
pntops.o: $(SRC)/pntops.f $(SRC)/XFOIL.INC $(SRC)/XDES.INC
	$(FC) -c $(FFLAGS)  $(SRC)/pntops.f
blplot.o: $(SRC)/blplot.f $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/blplot.f
xpol.o: $(SRC)/xpol.f     $(SRC)/XFOIL.INC
	$(FC) -c $(FFLAGS)  $(SRC)/xpol.f
xgeom.o: $(SRC)/xgeom.f 
	$(FC) -c $(FFLAGS)  $(SRC)/xgeom.f
xutils.o: $(SRC)/xutils.f 
	$(FC) -c $(FFLAGS)  $(SRC)/xutils.f
modify.o: $(SRC)/modify.f 
	$(FC) -c $(FFLAGS)  $(SRC)/modify.f
aread.o: $(SRC)/aread.f 
	$(FC) -c $(FFLAGS)  $(SRC)/aread.f
naca.o: $(SRC)/naca.f 
	$(FC) -c $(FFLAGS)  $(SRC)/naca.f
plutil.o: $(SRC)/plutil.f 
	$(FC) -c $(FFLAGS)  $(SRC)/plutil.f
userio.o: $(SRC)/userio.f 
	$(FC) -c $(FFLAGS)  $(SRC)/userio.f
gui.o: $(SRC)/gui.f 
	$(FC) -c $(FFLAGS)  $(SRC)/gui.f
spline.o: $(SRC)/spline.f 
	$(FC) -c $(FFLAGS)  $(SRC)/spline.f
sort.o: $(SRC)/sort.f
	$(FC) -c $(FFLAGS)  $(SRC)/sort.f
profil.o: $(SRC)/profil.f
	$(FC) -c $(FFLAGS)  $(SRC)/profil.f

polplt.o: $(SRC)/polplt.f $(SRC)/PINDEX.INC
	$(FC) -c $(FFLAGS)  $(SRC)/polplt.f
iopol.o: $(SRC)/iopol.f $(SRC)/PINDEX.INC
	$(FC) -c $(FFLAGS)  $(SRC)/iopol.f

frplot0.o: $(SRC)/frplot0.f
	$(FC) -c $(FFLAGS)  $(SRC)/frplot0.f
frplot.o: $(SRC)/frplot.f
	$(FC) -c $(FFLAGS)  $(SRC)/frplot.f
ntcalc.o: $(SRC)/ntcalc.f
	$(FC) -c $(FFLAGS)  $(SRC)/ntcalc.f

osmap.o: $(OSRC)/osmap.f
	$(FC) -c $(FFLAGS)  $(OSRC)/osmap.f

getosfile.o: $(OSRC)/getosfile.c
	$(CC) -c $(CFLAGS)  $(OSRC)/getosfile.c
