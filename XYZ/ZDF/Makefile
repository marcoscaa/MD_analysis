all:
	gfortran -I../src/ -O2 ZDF.f90 -o zdf.x ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 ZDF.f90 -o zdf.x ../src/libanalysis.o
