all:
	ifort -O2 -c module_pt.f90
	ifort -O2 proton_transfer.f90 -o pt.x module_pt.o
	ifort -O2 proton_transport.f90 -o ptt.x module_pt.o
	ifort -O2 coord_numb.f90 -o cn.x 
debug:
	ifort -O0 -g -traceback -check all -nofor-main -c module_pt.f90  
	ifort -O0 -g -traceback -check all proton_transfer.f90 -o pt.x module_pt.o 
	ifort -O0 -g -traceback -check all proton_transport.f90 -o ptt.x module_pt.o 
	ifort -O0 -g -traceback -check all coord_numb.f90 -o cn.x 
debug-gnu:
	gfortran -I../src/ -O0 -g -fbacktrace -fcheck=all proton_transfer.f90 -o pt.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 proton_transfer.f90 -o pt.x ../src/libanalysis.o
