all:
	gfortran -O2 -I../src polyhedral.f90 -o poly.x ../src/libanalysis.o 
debug: 
	gfortran -O0 -I../src  -g polyhedral.f90 -o poly.x ../src/libanalysis.o
