all:
	gfortran -c -cpp libanalysis.f90
	gfortran -I./ -O2 ZDF.f90 -o zdf.x ./libanalysis.o
	gfortran -I./ -O2 VAC_layers.f90 -o vac_l.x ./libanalysis.o
debug:
	gfortran -c -cpp -O0 -g -fcheck=all -fbacktrace libanalysis.f90
	gfortran -I./ -O0 -g -fbacktrace ZDF.f90 -o zdf.x ./libanalysis.o
	gfortran -I./ -fopenmp -O0 -g -fbacktrace VAC_layers.f90 -o vac_l.x ./libanalysis.o
