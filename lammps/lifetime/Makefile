all:
	ifort -O2 lifetime.f90 -o lt.x
gnu:
	gfortran -O2 -fopenmp lifetime.f90 -o lt.x
debug:
	ifort -O0 -g -check all -traceback lifetime.f90 -o lt.x
