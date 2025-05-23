all:
	gfortran -O2 SofK.f90 -o sofk.x 
cuda:
	pgf90 -O2 -Mcuda=cc60 SofK_cuda_Faber_Ziman.f90 -o sofk_c.x 
debug:
	gfortran -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan SofK.f90 -o sofk.x
