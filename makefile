MKLROOT=/export/home/jjren/apps/intel/mkl
MKLLIB=$(MKLROOT)/lib/intel64
mklinc=$(MKLROOT)/include/intel64/lp64 
mklinc1=$(MKLROOT)/include

FCCFLAG= -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread -lm

FC=mpiifort
FCCOMPILEOPTS= -O3

.SUFFIXES: .f90 .f

%.o : %.f90
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
%.o : %.f
	$(FC) -c $(FCCOMPILEOPTS) -I/$(mklinc) -I/$(mklinc1) $<
object = kinds_mod.o para_mod.o array_mod.o communicate.o exit_mod.o   \
	   readinput.o Ising.o main.o 
Ising : $(object)
	$(FC) -o $@ $^ -L$(MKLLIB) $(FCCFLAG)

clean:
	rm -f *.o *.mod Ising

	

