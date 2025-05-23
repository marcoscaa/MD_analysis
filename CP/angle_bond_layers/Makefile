all:
	ifort -I../src/ -O2 angbond.f90 -o angbond.x ../src/libmdanalysis.o
	#ifort -O2 angbond_gro.f90 -o angbond_gro.x
debug:
	ifort -I../src/ -O0 -g -traceback -check all angbond.f90 -o angbond.x  ../src/libmdanalysis.o
	#ifort -O0 -g -traceback -check all angbond_gro.f90 -o angbond_gro.x
