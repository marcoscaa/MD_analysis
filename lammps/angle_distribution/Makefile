all:
	ifort -I../src/ -O2 angle_dist.f90 -o angle_dist.x ../src/libanalysis.o
	ifort -I../src/ -O2 angle_dist_cnt.f90 -o angle_dist_cnt.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 angle_dist.f90 -o angle_dist.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback angle_dist.f90 -o angle_dist.x ../src/libanalysis.o
