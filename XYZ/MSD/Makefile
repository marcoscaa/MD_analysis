all:
	gfortran -I../src/ -O2 msd.f90 -o msd.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all msd.f90 -o msd.x ../src/libanalysis.o
