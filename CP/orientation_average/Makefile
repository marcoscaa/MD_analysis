all:
	ifort -I../src/ -O2 orientation.f90 -o orient.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all orientation.f90 -o orient.x ../src/libmdanalysis.o
