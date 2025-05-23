all:
	ifort -I../src/ -O2 hb_angle.f90 -o hbang.x ../src/libanalysis.o 
	ifort -I../src/ -O2 hb_angle_decomposed.f90 -o hbangdec.x ../src/libanalysis.o 
	ifort -I../src/ -O2 hb_distribution.f90 -o hbdist.x ../src/libanalysis.o 
gnu:
	gfortran -I../src/ -O2 hb_distribution.f90 -o hbdist.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 hb_distribution_oligoamines.f90 -o hbdist_oa.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 hb_lifetime_oligoamines.f90 -o hblife_oa.x ../src/libanalysis.o 
	gfortran -I../src/ -O2 hb_lifetime.f90 -o hblife.x ../src/libanalysis.o 
debug:
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan hb_angle.f90 -o hbang.x ../src/libanalysis.o
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan hb_lifetime.f90 -o hblife.x ../src/libanalysis.o
	gfortran -I../src/ -O0 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan hb_lifetime_oligoamines.f90 -o hblife_oa.x ../src/libanalysis.o
debug-intel:
	ifort -I../src/ -O0 -g -traceback -check all hb_distribution.f90 -o hbdist.x ../src/libanalysis.o 
