all:
	gfortran -I../src/ -O2 coord.f90 -o coord.x  ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 coord.f90 -o coord.x  ../src/libanalysis.o
