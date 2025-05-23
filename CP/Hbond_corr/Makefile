all:
	#ifort -O2 Hbond.f90 -o hbond.x
	#ifort -O2 Hbond_surf_cp.f90 -o hbond_surf_cp.x
	ifort -I../src/ -O0 -qopenmp Hbond_cp.f90 -o hbond_cp.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -qopenmp Hbond_cp_nodiff.f90 -o hbond_cp_nodiff.x ../src/libmdanalysis.o
debug:
	#ifort -O0 -g -traceback -check all Hbond.f90 -o hbond.x
	#ifort -O0 -g -traceback -check all Hbond_surf_cp.f90 -o hbond_surf_cp.x
	ifort -I../src/ -O0 -g -traceback -check all Hbond_cp.f90 -o hbond_cp.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -g -traceback -check all Hbond_cp_nodiff.f90 -o hbond_cp_nodiff.x ../src/libmdanalysis.o
