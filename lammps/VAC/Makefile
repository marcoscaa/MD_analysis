gnu:
	gfortran -I../src/ -fopenmp -O2 VAC.f90 -o vac.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O2 VAC_layers.f90 -o vac_l.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -fopenmp -O0 -g -fbacktrace VAC_layers.f90 -o vac_l.x ../src/libanalysis.o
