all:
	ifort -I../src/ -O2 op.f90 -o op.x ../src/libmdanalysis.o
gnu:
	gfortran -I../src/ -O2 -fopenmp op.f90 -o op.x ../src/libmdanalysis.o
debug: 
	ifort -I../src/ -O0 -g -traceback -check all -fpe0 op.f90 -o op.x ../src/libmdanalysis.o
