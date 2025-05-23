all:
	gfortran -I../src/ -O2 rotcorr_radial_distance.f90 -o rotcorr_r.x ../src/libanalysis.o
	gfortran -I../src/ -O2 rotcorr.f90 -o rotcorr.x ../src/libanalysis.o
	gfortran -I../src/ -O2 rotcorr_disctrtized_molecular_plane.f90 -o rotcorr_dmp.x ../src/libanalysis.o
	gfortran -I../src/ -O2 rotcorr_layers.f90 -o rotcorr_l.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all rotcorr.f90 -o rotcorr.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all rotcorr_radial_distance.f90 -o rotcorr_r.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all rotcorr_disctrtized_molecular_plane.f90 -o rotcorr_dmp.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace -fcheck=all rotcorr_layers.f90 -o rotcorr_l.x ../src/libanalysis.o
