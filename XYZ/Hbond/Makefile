all:
	gfortran -I../src/ -O2 hb_distribution_oligoamines.f90 -o hb.x ../src/libanalysis.o
	gfortran -I../src/ -O2 hb_total_distribution_oligoamines.f90 -o hbt.x ../src/libanalysis.o
	gfortran -I../src/ -O2 hb_average_surface_distance.f90 -o hbsurf.x ../src/libanalysis.o
	gfortran -I../src/ -O2 hb_average_surface_distance_bivo4.f90 -o hbsurfbvo.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all hb_distribution_oligoamines.f90 -o hb.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all hb_total_distribution_oligoamines.f90 -o hbt.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all hb_average_surface_distance.f90 -o hbsurf.x ../src/libanalysis.o
