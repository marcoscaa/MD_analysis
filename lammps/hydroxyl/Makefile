all:
	ifort -O2 coord_numb.f90 -o cn.x
	ifort -O2 count_hydroxyl_groups.f90 -o ch.x
gnu:
	gfortran -O2 coord_numb.f90 -o cn.x
	gfortran -O2 count_hydroxyl_groups.f90 -o ch.x
debug:
	ifort -O0 -g -traceback -check all coord_numb.f90 -o cn.x
