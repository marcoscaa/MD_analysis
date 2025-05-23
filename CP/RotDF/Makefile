all:
	ifort -I../src/ -O2 Rot_DF.f90 -o rotdf.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback -check all Rot_DF.f90 -o rotdf.x ../src/libmdanalysis.o
