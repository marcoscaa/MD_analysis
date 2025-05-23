all:
	gfortran -I../src/ -fopenmp -O2 SFG.f90 -o sfg.x ../src/libanalysis.o
	gfortran -I../src/ -fopenmp -O2 SFG_layers.f90 -o sfg_l.x ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -fopenmp -I../src/ -O0 SFG.f90 -o sfg.x ../src/libanalysis.o
	gfortran -g -fbacktrace -fcheck=all -fopenmp -I../src/ -O0 SFG_layers.f90 -o sfg_l.x ../src/libanalysis.o
