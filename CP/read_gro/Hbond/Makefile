all:
	ifort -O2 Hbond_layers.f90 -o hb.x 
	ifort -O2 Hbond_layers_orient.f90 -o hbor.x 
debug:
	ifort -O0 -g -traceback -check all Hbond_layers.f90 -o hb.x
	ifort -O0 -g -traceback -check all Hbond_layers_orient.f90 -o hbor.x 
