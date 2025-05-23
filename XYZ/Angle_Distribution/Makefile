all:
	ifort -I../src/ -O2 angle_dist.f90 -o angdist.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 angle_dist.f90 -o angdist.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback angle_dist.f90 -o angdist.x ../src/libanalysis.o
gnu-debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 angle_dist.f90 -o angdist.x ../src/libanalysis.o
