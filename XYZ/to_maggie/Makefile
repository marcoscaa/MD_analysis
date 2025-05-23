all:
	gfortran -c -cpp libanalysis.f90
	gfortran -I./ -O2 compute_water_ion_position.f90 -o pos_wi.x ./libanalysis.o
	gfortran -I./ -O2 Hbond.f90 -o hb_wi.x ./libanalysis.o
