all:
	ifort -I../src/ -O0 -qopenmp SFG.f90 -o sfg.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all SFG.f90 -o sfg.x ../src/libanalysis.o
