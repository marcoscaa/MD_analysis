gnu:
	gfortran -c -cpp libanalysis.f90
	gfortran -I./ -O2 proton_transfer.f90 -o pt.x ./libanalysis.o
