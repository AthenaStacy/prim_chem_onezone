GCC = g++
CPPFLAGS = -O2

all: PrimDriver 

PrimDriver: PrimChemDriver.o PrimChemBurner.o PrimChemIntegrate.o Jacobian.o 
	$(GCC) $(CPPFLAGS) -o PrimDriver PrimChemDriver.o PrimChemBurner.o PrimChemIntegrate.o Jacobian.o 

PrimChemBurner.o : PrimChemBurner.cpp PrimChemBurner.H PrimChemGlobals.H
	$(GCC) $(CPPFLAGS) -c PrimChemBurner.cpp

PrimChemDriver.o : PrimChemDriver.cpp PrimChemDriver.H PrimChemGlobals.H
	$(GCC) $(CPPFLAGS) -c PrimChemDriver.cpp

PrimChemIntegrate.o : PrimChemIntegrate.cpp PrimChemIntegrate.H PrimChemGlobals.H 
	$(GCC) $(CPPFLAGS) -c PrimChemIntegrate.cpp

Jacobian.o : Jacobian.cpp PrimChemGlobals.H
	$(GCC) $(CPPFLAGS) -c Jacobian.cpp

clean:
	rm *.o PrimDriver
