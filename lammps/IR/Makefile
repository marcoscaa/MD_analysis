all:
	gfortran -I../src/ -fopenmp -O2 Dipole_time_corr.f90 -o dtc.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace Dipole_time_corr.f90 -o dtc.x ../src/libanalysis.o
