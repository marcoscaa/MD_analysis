all:
	ifort -O2 -I../src/ -qopenmp Dipole_Polarizability_time_correlation_cutoff.f90 -o dptc.x ../src/libanalysis.o
	ifort -O2 -I../src/ -qopenmp Global_Dipole_Polarizability_time_correlation_cutoff.f90 -o gdptc.x ../src/libanalysis.o
	ifort -O2 -I../src/ -qopenmp Compute_total_dipole.f90 -o tdip.x ../src/libanalysis.o
gnu:
	gfortran -O2 -I../src/ -fopenmp Dipole_Polarizability_time_correlation_cutoff.f90 -o dptc.x ../src/libanalysis.o
	gfortran -O2 -I../src/ -fopenmp Compute_total_dipole.f90 -o tdip.x ../src/libanalysis.o
