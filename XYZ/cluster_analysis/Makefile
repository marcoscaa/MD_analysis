all:
	gfortran -I../src/ -O2 cluster.f90 -o cluster.x  ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 cluster.f90 -o cluster.x  ../src/libanalysis.o
