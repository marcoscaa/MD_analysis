all:
	ifort -qopenmp -O0 vac_cp_bulk.f90 -o vac.x 
bulk:
	ifort -qopenmp -O2 vac_cp_bulk.f90 -o vac.x 
debug:
	#ifort -qopenmp -xHost -g -O2 vac_cp.f90 -o vac.x 
	#ifort -I../src/ -qopenmp -O0 -g -traceback -check all vac_cp_bulk.f90 -o vac.x ../src/libmdanalysis.o
	ifort -qopenmp -O0 -g -traceback -check all vac_cp_bulk.f90 -o vac.x
