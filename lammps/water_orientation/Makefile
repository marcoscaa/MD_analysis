all:
	gfortran -I../src/ -O2 orientation_cnt.f90 -o ohangle.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 orientation_cnt_2D.f90 -o ohangle2d.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 average_orientation_distribution.f90 -o averorient.x ../src/libanalysis.o
debug:
	gfortran -I../src/ -O0 -g -traceback -check all orientation_cnt.f90 -o ohangle.x ../src/libanalysis.o 
