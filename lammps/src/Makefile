gnu:
	gfortran -c -cpp libanalysis.f90
gnu-debug:
	gfortran -c -cpp -O0 -g -fcheck=all -fbacktrace libanalysis.f90
cuda: 
	pgif90 -cpp -D__CUDA -Mcuda=cc70 libanalysis.f90  
intel:
	ifort -c -cpp libanalysis.f90
