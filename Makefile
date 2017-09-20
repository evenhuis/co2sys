# makefile


#FC= gfortran-mp-4.5 -O2   -Warray-bounds -O0 -g -fbounds-check
FC= gfortran -std=f2008 -O2 # -Warray-bounds -O0 -g -fbounds-check


test_co2 : test_co2.f90 carbonate_chemistry_mod.f90
	$(FC) $(BLAS)  -c carbonate_chemistry_mod.f90
	$(FC) $(BLAS)  -o test_co2 test_co2.f90 carbonate_chemistry_mod.o

cc_py_mod.so : cc_py_mod.f90 
	touch cc_py_mod.so
	rm    cc_py_mod.so
	$(FC) $(BLAS)  -c carbonate_chemistry_mod.f90
	f2py -c -m cc_py_mod cc_py_mod.f90 carbonate_chemistry_mod.o


