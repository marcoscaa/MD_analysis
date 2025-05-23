all:
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan define_molecule.f90 -o dm.x ../src/libanalysis.o
