all:
	#gfortran -I../../lammps/src/ -O2 op.f90 -o op.x ../../lammps/src/libanalysis.o
	gfortran -I../src/ -O2 compute_water_ion_position.f90 -o pos_wi.x ../src/libanalysis.o
	#gfortran -I../src/ -O2 Hbond.f90 -o hb_wi.x ../src/libanalysis.o
	#gfortran -I../../lammps/src/ -O2 msd.f90 -o msd.x ../../lammps/src/libanalysis.o
	gfortran -I../src/ -O2 P_OOO.f90 -o pooo.x ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -fopenmp -I../src/ -O0 Hbond.f90 -o hb_wi.x ../src/libanalysis.o
