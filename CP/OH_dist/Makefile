all:
	ifort -I../src/ -O0 -qopenmp  OH_dist.f90 -o oh_dist.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback OH_dist.f90 -o oh_dist.x ../src/libmdanalysis.o
