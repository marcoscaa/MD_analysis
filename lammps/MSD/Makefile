all:
	ifort -I../src/ -O2 -fopenmp msd.f90 -o msd.x ../src/libanalysis.o
	ifort -I../src/ -O2 -fopenmp msd_xy_nolayers.f90 -o msdxy.x ../src/libanalysis.o
	ifort -I../src/ -O2 -fopenmp msd_z.f90 -o msdz.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 msd_radial_distance.f90 -o msdr.x ../src/libanalysis.o
	gfortran -I../src/ -O2 msd.f90 -o msd.x ../src/libanalysis.o
	gfortran -I../src/ -O2 msd_disctrtized_molecular_plane.f90 -o msd_dmp.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -fbacktrace -fdump-core msd.f90 -o msd.x ../src/libanalysis.o
gnu-debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all msd_radial_distance.f90 -o msdr.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all msd_disctrtized_molecular_plane.f90 -o msd_dmp.x ../src/libanalysis.o
