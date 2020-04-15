##################################
# makefile: otpoc stands for Optical and Transport Properties of Organic Microcavities
##################################

FC 	= mpif90  #gfortran
FFLAGS= -O3 -g -fbounds-check -fbacktrace
LOADLIBES = libarpack_mac.a fftlib.a -llapack

# to link statically we might need something like: 
#FFLAGS += -static
#LOADLIBES += -lblas -lgfortran -lm


SRC = main.f basis.f modmain.f maps.f hamiltonian.f diag.f correlation.f hamiltonian1.f correlationd.f dm.f

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

otpoc: $(OBJ)
	$(FC) -o otpoc $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp otpoc ~/bin/.


# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.
 
modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 

maps.o: maps.f modmain.o
	 $(FC) -c $(FFLAGS) maps.f 

basis.o: basis.f modmain.o maps.o
	 $(FC) -c $(FFLAGS) basis.f 	

hamiltonian.o: hamiltonian.f modmain.o maps.o basis.o 
	 $(FC) -c $(FFLAGS) hamiltonian.f 	

hamiltonian1.o: hamiltonian1.f hamiltonian.f modmain.o maps.o basis.o 
	 $(FC) -c $(FFLAGS) hamiltonian1.f 	


diag.o: diag.f modmain.o
	 $(FC) -c $(FFLAGS) diag.f 	

correlation.o: correlation.f modmain.o 
	 $(FC) -c $(FFLAGS) correlation.f 	

correlationd.o: correlationd.f modmain.o correlation.o
	 $(FC) -c $(FFLAGS) correlationd.f 	

dm.o: dm.f modmain.o
	 $(FC) -c $(FFLAGS) dm.f 


main.o: main.f modmain.o maps.o basis.o hamiltonian.o diag.o correlation.o  hamiltonian1.o correlationd.o dm.o
	 $(FC) -c $(FFLAGS) main.f	

clean:
	rm -rf *.o *.mod otpoc
  
