

FC= gfortran -std=f2008 
FC=gfortran  -O2 # -Warray-bounds -O0 -g -fbounds-check
BLAS= -L/usr/lib -lblas -llapack
NETCDF_FLAG=-I/opt/local/include/ -L/opt/local/lib -lnetcdff

LIB_DIR = /Users/evenhuis/Dropbox/fortran_programs/lib/
LIB_O   = -I$(LIB_DIR) -L$(LIB_DIR) -lmylib

% : %.f90 
	$(FC) -o $* $*.f90  $(BLAS) #$(NETCDF_FLAG)
	cp $* ../bin

carbon_test:            carbon_test.f90         carbonate_chemistry_mod.o PBR_mod.o ODE_mod.o ExpoKitModule.o 
	$(FC) $(BLAS) -o carbon_test carbon_test.f90 carbonate_chemistry_mod.o PBR_mod.o ODE_mod.o ExpoKitModule.o
	cp carbon_test ../bin

# Library modules
carbonate_chemistry_mod.o: carbonate_chemistry_mod.f90
	$(FC) -c                carbonate_chemistry_mod.f90
	cp $*.o   ../lib
	cp $*.mod ../lib
	../lib/archive_library.sh

ODE_mod.o: ODE_mod.f90  ExpoKitModule.o
	$(FC) -c ODE_mod.f90 ExpoKitModule.o

ExpoKitModule.o : ExpoKitModule.f
	$(FC) $(BLAS) -c ExpoKitModule.f

PBR_mod.o :          PBR_mod.f90 carbonate_chemistry_mod.o ODE_mod.o ExpoKitModule.o
	$(FC) $(BLAS) -c  PBR_mod.f90 carbonate_chemistry_mod.o ODE_mod.o ExpoKitModule.o

linest_mod.o: linest_mod.f90 
	$(FC) -c linest_mod.f90  

m_rnkpar.o : m_rnkpar.f90
	$(FC) -c  m_rnkpar.f90

m_mrgrnk.o : m_mrgrnk.f90
	$(FC) -c  m_mrgrnk.f90

