all:
	#ifort -O2 ZDF.f90 -o zdf.x
	ifort -I../src/ -O2 ZDF_cp.f90 -o zdf_cp.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -traceback -check all -g ZDF_cp.f90 -o zdf_cp.x ../src/libmdanalysis.o
