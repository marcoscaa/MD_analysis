all:
	ifort -I../src -O2  wan2dip.f90 -o wan2dip.x ../src/libanalysis.o
gnu:
	gfortran -I../src -O2  wan2dip_eigen.f90 -o wan2dip_eigen.x ../src/libanalysis.o
debug:
	ifort -I../src -O0 -g -traceback -check all wan2dip.f90 -o wan2dip.x ../src/libanalysis.o
