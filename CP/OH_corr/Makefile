all:
	ifort -I../src/ -O0 -qopenmp  OH_corr.f90 -o oh_corr.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all OH_corr.f90 -o oh_corr.x ../src/libmdanalysis.o
