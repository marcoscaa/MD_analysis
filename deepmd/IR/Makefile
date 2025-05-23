all:
	ifort -I../src -O2 -qopenmp Dipole_time_corr.f90 -o dac.x ../src/libanalysis.o
debug:
	ifort -I../src -qopenmp -O0 -g -traceback -check all Dipole_time_corr.f90 -o dac.x ../src/libanalysis.o
