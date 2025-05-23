all:
	ifort -I../src/ -O2 ZDF.f90 -o zdf.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 ZDF.f90 -o zdf.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback ZDF.f90 -o zdf.x ../src/libanalysis.o
gnu-debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 ZDF.f90 -o zdf.x ../src/libanalysis.o
