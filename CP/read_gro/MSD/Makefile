all:
	ifort -qopenmp -xHost -g -O3 msd.f90 -o msd.x 
	ifort -qopenmp -xHost -g -O3 msd_filter.f90 -o filter_msd.x 
	#gfortran -fopenmp -g -O3 read_gro.f90 
