all:
	ifort -I../src -O2 ZDF.f90 -o zdf.x ../src/libanalysis.o
	ifort -I../src -O2 ZDF_npy.f90 -o zdf_npy.x ../src/libanalysis.o
gnu:
	gfortran -I../src -O2 ZDF.f90 -o zdf.x ../src/libanalysis.o
	#gfortran -I../src -O2 ZDF_npy.f90 -o zdf_npy.x ../src/libanalysis.o
debug:
	ifort -I../src -O0 -g -traceback -check all ZDF.f90 -o zdf.x ../src/libanalysis.o
	ifort -I../src -O0 -g -traceback -check all ZDF_npy.f90 -o zdf_npy.x ../src/libanalysis.o
