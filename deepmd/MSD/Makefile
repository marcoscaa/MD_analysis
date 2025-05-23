all:
	ifort -O2 -I../src/ -qopenmp msd.f90 -o msd.x ../src/libanalysis.o
gnu:
	gfortran -O2 -I../src/ -fopenmp msd.f90 -o msd.x ../src/libanalysis.o
