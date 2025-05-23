all:
	#ifort -qopenmp -xHost -g -O3 oh_corr.f90 -o oh_corr.x 
	ifort -I../src/ -qopenmp -g -O2 rotcorr_cp.f90 -o rotcorr.x ../src/libmdanalysis.o
	#ifort -qopenmp -O0 -g -traceback -check all rotcorr_cp.f90 -o dipole_corr.x
debug:
	ifort -I../src/ -O0 -g -traceback -check all -fpe0 rotcorr_cp.f90 -o rotcorr.x ../src/libmdanalysis.o
