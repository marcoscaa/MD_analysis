all:
	ifort -O2 -c module_pt.f90
	ifort -O2 proton_transfer.f90 -o pt.x module_pt.o
debug:
	ifort -O0 -g -traceback -check all -nofor-main -c module_pt.f90  
	ifort -O0 -g -traceback -check all proton_transfer.f90 -o pt.x module_pt.o 
test:
	ifort -O2 -c module_pt.f90
	ifort -O2 proton_transfer.f90 -o pt.x2 module_pt.o
