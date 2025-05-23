all:
	ifort -I../src/ -O2 restime.f90 -o restime.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all restime.f90 -o restime.x ../src/libmdanalysis.o
