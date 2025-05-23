all:
	gfortran -c -cpp libanalysis.f90 
	gfortran -I./ -O2 msd.f90 -o msd.x ./libanalysis.o
