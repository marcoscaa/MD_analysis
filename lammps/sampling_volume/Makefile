all:
	gfortran -I../src/ -O2 -fopenmp Sampling_volume.f90 -o sv.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -O0 -g -fopenmp Sampling_volume.f90 -o sv.x ../src/libanalysis.o
