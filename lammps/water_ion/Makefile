all:
	gfortran -I../src/ -O2 Hbond.f90 -o hb.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 P_OOO.f90 -o pooo.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 print_environment_around_water_ion.f90 -o print.x ../src/libanalysis.o 
debug:
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan Hbond.f90 -o hb.x ../src/libanalysis.o
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan P_OOO.f90 -o pooo.x ../src/libanalysis.o
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan print_environment_around_water_ion.f90 -o print.x ../src/libanalysis.o
