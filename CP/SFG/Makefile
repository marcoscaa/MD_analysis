all:
	ifort -I../src/ -O0 -qopenmp sfg2.f90 -o sfg.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback sfg2.f90 -o sfg.x ../src/libmdanalysis.o
old:
	ifort -I../src/ -O0 -g -traceback -check all sfg.f90 -o sfg.x ../src/libmdanalysis.o
