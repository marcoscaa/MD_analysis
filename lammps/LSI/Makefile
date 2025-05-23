all:
	ifort -O2 LSI_water_ion.f90 -o lsi_wi.x 
	ifort -O2 LSI.f90 -o lsi.x 
gnu:
	gfortran -O2 LSI_water_ion.f90 -o lsi_wi.x 
	gfortran -O2 LSI.f90 -o lsi.x 
debug:
	ifort -O0 -g -traceback -check all LSI_water_ion.f90 -o lsi_wi.x
	ifort -O0 -g -traceback -check all LSI.f90 -o lsi.x
debug-gnu:
	gfortran -O0 -g -fbacktrace -fcheck=all LSI_water_ion.f90 -o lsi.x 
