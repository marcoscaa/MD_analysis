all:
	gfortran -I../src -O2  wan_hbond.f90 -o wanhb.x ../src/libanalysis.o
debug:
	gfortran -I../src -O0 -fbacktrace wan_hbond.f90 -o wanhb.x ../src/libanalysis.o
