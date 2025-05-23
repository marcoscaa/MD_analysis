all:
	ifort -I../src -O2 Polarizability_time_correlation.f90 -o ptc.x ../src/libanalysis.o
	ifort -I../src -O2 Total_polarizability_time_correlation.f90 -o tptc.x ../src/libanalysis.o
	ifort -I../src -O2 Polarizability_time_correlation_cutoff.f90 -o ptcc.x ../src/libanalysis.o
	#ifort -I../src -O2 Polarizability_time_correlation_cutoff_select.f90 -o ptccs.x ../src/libanalysis.o
	ifort -I../src -O2 Polarizability_time_correlation_no_time_average.f90 -o ptccs.x ../src/libanalysis.o
debug:
	ifort -I../src -O0 -g -traceback -check all Polarizability_time_correlation.f90 -o ptc.x ../src/libanalysis.o
