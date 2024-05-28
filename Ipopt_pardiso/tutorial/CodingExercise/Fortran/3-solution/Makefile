# Copyright (C) 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the four CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = TutorialFortran

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = TutorialFortran.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# Fortran Compiler options
F77 = gfortran

# Fotran Compiler options
FFLAGS = -O2 

# additional Fortran Compiler options for linking
F77LINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib

prefix=/usr/local
exec_prefix=${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)
#INCL = -I${prefix}/include/coin -I/usr/local/include/coin-or/hsl -I/usr/local/include/coin-or/mumps    -DIPOPTLIB_BUILD $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt` -lstdc++ -lm
#LIBS = -L${exec_prefix}/lib -lipopt -L/usr/local/lib -lcoinhsl -lcoinmumps -llapack -lblas    -lm  -ldl -lstdc++ -lm

all: $(EXE)

.SUFFIXES: .f .o

$(EXE): $(OBJS)
	$(F77) $(F77LINKFLAGS) $(FFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -f $(EXE) $(OBJS) IPOPT.OUT

.f.o:
	$(F77) $(FFLAGS) $(INCL) -c -o $@ $<
