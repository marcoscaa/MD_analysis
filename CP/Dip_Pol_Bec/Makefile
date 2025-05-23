all:
	ifort -I../src/ -O2 Wannier.f90 -o wannier.x ../src/libmdanalysis.o
	#ifort -I../src/ -O2 Wannier_centraldiff.f90 -o wannier_cd.x ../src/libmdanalysis.o
debug:
	#ifort -I../src/ -O0 -g -traceback -check all Wannier.f90 -o wannier.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -g  Wannier.f90 -o wannier.x ../src/libmdanalysis.o
