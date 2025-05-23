all:
	ifort -I../src/ -O2 op.f90 -o op.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -traceback -check all -g op.f90 -o op.x ../src/libanalysis.o
