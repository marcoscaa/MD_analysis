all:
	#ifort -I../src/ -O2 Hbond_layers_cp.f90 -o hb.x ../src/libmdanalysis.o
	ifort -I../src/ -O0 -g -traceback Hbond_layers_cp_2ndlayer.f90 -o hb_2l.x ../src/libmdanalysis.o
	#ifort -O2 Hbond_layers.f90 -o hb_cp.x
debug:
	ifort -I../src/ -O0 -g -traceback -check all Hbond_layers_cp.f90 -o hb.x ../src/libmdanalysis.o
