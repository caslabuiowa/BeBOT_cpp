# Copyright (C) 2005, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Directory to the sources for the (example) problem definition files
SRCDIR = .


# Compiler and options
CXX = g++
CXXFLAGS = -O2 -DNDEBUG 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib

# get some directories, so we can expand ${exec_prefix}/lib for -rpath arg for libtool
prefix=/usr/local
exec_prefix=${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt`
#INCL = -I${prefix}/include/coin-or -I/usr/local/include/coin-or/hsl -I/usr/local/include/coin-or/mumps    -DIPOPTLIB_BUILD

# Linker flags
LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt`
#LIBS = -L${exec_prefix}/lib -lipopt -L/usr/local/lib -lcoinhsl -lcoinmumps -llapack -lblas    -lm  -ldl

# Name of the executable
EXE = solve_problem

# List of all object files
MAINOBJ =  solve_problem.o \
        MittelmannDistCntrlNeumA.o \
        MittelmannDistCntrlNeumB.o \
        MittelmannDistCntrlDiri.o \
        MittelmannBndryCntrlDiri.o \
        MittelmannBndryCntrlDiri3D.o \
        MittelmannBndryCntrlDiri3D_27.o \
        MittelmannBndryCntrlDiri3Dsin.o \
        MittelmannBndryCntrlNeum.o \
        LuksanVlcek1.o \
        LuksanVlcek2.o \
        LuksanVlcek3.o \
        LuksanVlcek4.o \
        LuksanVlcek5.o \
        LuksanVlcek6.o \
        LuksanVlcek7.o \
        RegisteredTNLP.o

all: $(EXE)

# Dependencies
LuksanVlcek1.o : LuksanVlcek1.hpp
LuksanVlcek2.o : LuksanVlcek2.hpp
LuksanVlcek3.o : LuksanVlcek3.hpp
LuksanVlcek4.o : LuksanVlcek4.hpp
LuksanVlcek5.o : LuksanVlcek5.hpp
LuksanVlcek6.o : LuksanVlcek6.hpp
LuksanVlcek7.o : LuksanVlcek7.hpp
MittelmannDistCntrlNeumA.o: MittelmannDistCntrlNeumA.hpp
MittelmannDistCntrlNeumB.o: MittelmannDistCntrlNeumB.hpp
MittelmannDistCntrlDiri.o: MittelmannDistCntrlDiri.hpp
MittelmannBndryCntrlDiri.o: MittelmannBndryCntrlDiri.hpp
MittelmannBndryCntrlDiri3D.o: MittelmannBndryCntrlDiri3D.hpp
MittelmannBndryCntrlDiri3D_27.o: MittelmannBndryCntrlDiri3D_27.hpp
MittelmannBndryCntrlDiri3Dsin.o: MittelmannBndryCntrlDiri3Dsin.hpp
MittelmannBndryCntrlNeum.o: MittelmannBndryCntrlNeum.hpp

solve_problem.o: \
	LuksanVlcek1.hpp \
	LuksanVlcek2.hpp \
	LuksanVlcek3.hpp \
	LuksanVlcek4.hpp \
	LuksanVlcek5.hpp \
	LuksanVlcek6.hpp \
	LuksanVlcek7.hpp \
	MittelmannDistCntrlNeumA.hpp \
	MittelmannDistCntrlNeumB.hpp \
	MittelmannDistCntrlDiri.hpp \
	MittelmannBndryCntrlDiri.hpp \
	MittelmannBndryCntrlDiri3D.hpp \
	MittelmannBndryCntrlDiri3D_27.hpp \
	MittelmannBndryCntrlDiri3Dsin.hpp \
	MittelmannBndryCntrlNeum.hpp \
	MittelmannParaCntrl.hpp

.SUFFIXES: .cpp .o

$(EXE): $(MAINOBJ) $(LIB)
	$(CXX) $(CXXFLAGS) $(CXXLINKFLAGS) -o $@ $(MAINOBJ) $(LIBS)

clean:
	rm -rf $(EXE) $(MAINOBJ)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<
