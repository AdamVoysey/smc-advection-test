## Makefile created to compile the spectral propagation test
## program.  Options are special for this program  as it is 
## too large for normal compilation.   JGLi18Dec2017
## Syntax:  make -f Makefile 
##

FC = ftn
FCFLAGS = 
CPPFLAGS= 

WRfrctnHyb: Constants.o  W3DispMD.o  W3PSMCMD.o  WRefractnHyb.o 
	@ echo "Linking $@ $^"
	$(FC) $(FCFLAGS) -o $@ $^

Constants.o: Constants.f90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c Constants.f90

W3DispMD.o: Constants.o W3DispMD.f90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c W3DispMD.f90

W3PSMCMD.o: Constants.o W3PSMCMD.f90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c W3PSMCMD.f90

WRefractnHyb.o: Constants.o W3DispMD.o W3PSMCMD.o WRefractnHyb.f90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c WRefractnHyb.f90

%.o: %.F90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

%.o: %.f90
	@ echo "Compiling $<"
	$(FC) $(FCFLAGS) $(CPPFLAGS) -c $<

clean:
	@ echo "Cleaning up ..."
	rm -f *.o *.mod 

