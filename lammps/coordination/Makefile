all:
	gfortran -O2 -fopenmp coord_number.f90 -o coord.x 
	gfortran -O2 print_local_environment.f90 -o printl.x
debug: 
	gfortran -O0 -g coord_number.f90 -o coord.x
cuda:
	pgf90 -Mcuda -Mcuda=cc70,cuda10.0 coord_number_cuda.f90 -o coord_cuda.x
