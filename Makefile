CC=gcc
MPIFC=mpif90
FC=gfortran-4.8
FFLAGS= -cpp -O2 -ffast-math -ftree-vectorize  
CFLAGS= -O3 -std=c99
#FFLAGS= -cpp -O3 -DDEBUG 
#FFLAGS= -cpp -O3 -DUNFORMATTED 
#FFLAGS= -cpp -O3 -DUNFORMATTED -DBUFFERED
OPENMP=-fopenmp 
#FLIBS= -L$(EXTRAE_HOME)/lib -lomptrace 


C_OBJECTS=  srchndpushc.o 

GSLIB_HOME=gslib
GSLIB=$(GSLIB_HOME)/gslib.a


all: seq seq-lev-for par-lev-for 

seq: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -
	$(FC) $(FFLAGS) -c sgsim.for
	$(FC) $(FFLAGS) -I. sgsim.o -o sgsimFortranSeq.exe $(GSLIB)

seq-lev-for: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -
	$(FC) $(FFLAGS) -c sgsimLevels.for
	gcc-4.8 -c -O3 srchndpushc.c
	$(FC) $(FFLAGS) -I. sgsimLevels.o $(C_OBJECTS) -o sgsimFortranLevelsSeq.exe $(GSLIB)
#	$(FC) $(FFLAGS) -I. sgsimLevels.o -o sgsimFortranLevelsSeq.exe $(GSLIB)


par-lev-for: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(FC) $(FFLAGS) $(OPENMP) -c sgsimLevels.for
	gcc-4.8 -c -O3 srchndpushc.c
	$(FC) $(FFLAGS) $(OPENMP) -I. sgsimLevels.o $(C_OBJECTS) -o sgsimFortranLevelsPar.exe $(GSLIB)
#	$(FC) $(FFLAGS) $(OPENMP) -I. sgsimLevels.o -o sgsimFortranLevelsPar.exe $(GSLIB)

#openmp_mpi_par: 
#	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(MPIFC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -c sgsim.fpp
#	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -I. sgsim.o -o sgsim_OpenMP_MPI $(GSLIB)
#	cp sgsim_OpenMP_MPI ../bin 

#ins:
#	cd $(HOME)/gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -c extrae_module.f90
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -c sgsim.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE sgsim.o -o sgsim.exe $(GSLIB) $(FLIBS) -lgomp  

clean:
	rm *.exe *.o *.mod gslib/*.o gslib/*.a gslib/*.mod
