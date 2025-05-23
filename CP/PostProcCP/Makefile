all:
	ifort -O2 postproc.f90 -o postproc.x
debug:
	ifort -O0 -g -traceback -check all postproc.f90 -o postproc.x
