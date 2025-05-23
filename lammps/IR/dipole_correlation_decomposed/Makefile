all:
	gfortran -c -cpp libanalysis.f90 
	gfortran -I./ -O2 Dipole_time_corr.f90 -o dtc.x ./libanalysis.o
