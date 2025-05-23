all:
	ifort -O2 Dipole_DF_gro.f90 -o dipole.x 
debug:
	ifort -O0 -g -traceback -check all Dipole_DF_gro.f90 -o dipole.x
