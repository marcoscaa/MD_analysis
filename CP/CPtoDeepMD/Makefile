all:
	ifort -O2 cptodeepmd.f90 -o cptodeepmd.x 
	ifort -O2 cptodeepmd_virial.f90 -o cptodeepmd_v.x 
wfc:
	ifort -I../src -O0 -g cptodeepmd_wfc.f90 -o cptodeepmd_wfc.x ../src/libmdanalysis.o
wfcgnu:
	gfortran -I../src -O0 -g cptodeepmd_wfc.f90 -o cptodeepmd_wfc.x ../src/libmdanalysis.o
