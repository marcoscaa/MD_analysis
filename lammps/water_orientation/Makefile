all:
	ifort -I../src/ -O2 orientation_cnt.f90 -o ohangle.x ../src/libanalysis.o 
	ifort -I../src/ -O2 orientation_cnt_2D.f90 -o ohangle2d.x ../src/libanalysis.o 
debug:
	ifort -I../src/ -O0 -g -traceback -check all orientation_cnt.f90 -o ohangle.x ../src/libanalysis.o 
