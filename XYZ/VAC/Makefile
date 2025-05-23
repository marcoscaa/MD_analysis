all:
	 gfortran -I../src/ -O2 VAC_from_pos.f90 -o vac.x ../src/libanalysis.o
	 gfortran -I../src/ -O2 VAC_layers.f90 -o vacl.x ../src/libanalysis.o
debug:
	 gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 VAC_from_pos.f90 -o vac.x ../src/libanalysis.o
	 gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 VAC_layers.f90 -o vacl.x ../src/libanalysis.o
