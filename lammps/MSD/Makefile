all:
	gfortran -I../src/ -O2 msd_radial_distance.f90 -o msdr.x ../src/libanalysis.o
	gfortran -I../src/ -O2 msd.f90 -o msd.x ../src/libanalysis.o
	gfortran -I../src/ -O2 msd_xy.f90 -o msd_xy.x ../src/libanalysis.o
	gfortran -I../src/ -O2 msd_disctrtized_molecular_plane.f90 -o msd_dmp.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all msd_radial_distance.f90 -o msdr.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all msd_disctrtized_molecular_plane.f90 -o msd_dmp.x ../src/libanalysis.o
