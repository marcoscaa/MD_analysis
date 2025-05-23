all:
	ifort -I../src/ -O2 orientation_distribution.f90 -o orient_dist.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 orientation_distribution.f90 -o orient_dist.x ../src/libanalysis.o
	gfortran -I../src/ -O2 average_orientation_distribution.f90 -o ave_orient_dist.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback orientation_distribution.f90 -o orient_dist.x ../src/libanalysis.o
gnu-debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 orientation_distribution.f90 -o orient_dist.x ../src/libanalysis.o
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 average_orientation_distribution.f90 -o ave_orient_dist.x ../src/libanalysis.o
