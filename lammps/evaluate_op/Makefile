all:
	gfortran -I../src/ -O2 -fopenmp op.f90 -o op.x ../src/libanalysis.o
	gfortran -I../src/ -O2 -fopenmp msd.f90 -o msd.x ../src/libanalysis.o
	gfortran -I../src/ -O2 -fopenmp print_water_ions.f90 -o print_wi.x ../src/libanalysis.o
	gfortran -I../src/ -O2 -fopenmp hydrogen_bond_water_ions.f90 -o hb_wi.x ../src/libanalysis.o
	gfortran -I../src/ -O2 -fopenmp S_and_proton_wire_CVs.f90 -o s_and_pw.x ../src/libanalysis.o
	gfortran -I../src/ -O2 -fopenmp environment_around_water_reacting.f90 -o env_react.x ../src/libanalysis.o
debug: 
	gfortran -I../src/ -O0 -g op.f90 -o op.x ../src/libanalysis.o
	gfortran -I../src/ -O0 -g -fbacktrace S_and_proton_wire_CVs.f90 -o s_and_pw.x ../src/libanalysis.o
