all:
	#ifort -O2 ordpar_cp.f90 -o ordpar.x 
	#ifort -I../src/ -O2 ordpar_cp_dist.f90 -o ordpar_dist.x ../src/libmdanalysis.o
	ifort -I../src/ -O2 -qopenmp ordpar_q6.f90 -o op_ql.x ../src/libmdanalysis.o
	ifort -I../src/ -O2 ordpar_tetra.f90 -o op_tetra.x ../src/libmdanalysis.o
debug:
	#ifort -O0 -g -traceback -check all ordpar_cp.f90 -o ordpar.x
	#ifort -I../src/ -O0 -g -traceback -check all ordpar_cp_dist.f90 -o ordpar_dist.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -g -traceback -check all -fpe0 ordpar_q6.f90 -o op_ql.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -g -traceback -check all -fpe0 ordpar_tetra.f90 -o op_tetra.x ../src/libmdanalysis.o
