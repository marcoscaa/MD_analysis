all:
	gfortran -I../src/ -O2 ordpar_q6.f90 -o q6.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -O0 -g -fcheck=all -fbacktrace ordpar_q6.f90 -o q6.x ../src/libanalysis.o
