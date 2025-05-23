gnu:
	#gfortran -I../src/ -O2 dipole_distribution.f90 -o dipdist.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 dipole_distribution_radial.f90 -o dipdistR.x ../src/libanalysis.o 
	#gfortran -I../src/ -O2 mean_dipole_distribution_radial.f90 -o dipdistR_mean.x ../src/libanalysis.o 
	#gfortran -I../src/ -O2 select_polar_configurations.f90 -o select_dip.x ../src/libanalysis.o 
gnu-debug:
	gfortran -I../src/ -O0 -g -fbacktrace -fcheck=all dipole_distribution.f90 -o dipdist.x ../src/libanalysis.o 
	gfortran -I../src/ -O0 -g -fbacktrace -fcheck=all dipole_distribution_radial.f90 -o dipdistR.x ../src/libanalysis.o 
