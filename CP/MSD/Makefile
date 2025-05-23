all:
	ifort -I../src/ -qopenmp -O0 msd_cp.f90 -o msd_cp.x ../src/libmdanalysis.o
	ifort -I../src/ -qopenmp -O0 msd_cp_xy.f90 -o msd_cp_xy.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all  msd_cp.f90 -o msd_cp.x ../src/libmdanalysis.o
