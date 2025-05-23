all:
	ifort -I../src/ -O2 lmp_to_raw.f90 -o lmp_to_raw.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 lmp_to_raw.f90 -o lmp_to_raw.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback lmp_to_raw.f90 -o lmp_to_raw.x ../src/libanalysis.o
