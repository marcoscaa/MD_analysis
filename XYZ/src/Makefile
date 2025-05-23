gnu:
	gfortran -c -cpp libanalysis.f90
debug:
	gfortran -g -fbacktrace -fcheck=all -c -cpp libanalysis.f90
cuda: 
	pgif90 -cpp -D__CUDA -Mcuda=cc70 libanalysis.f90  
intel:
	ifort -c -cpp libanalysis.f90
