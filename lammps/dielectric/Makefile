all:
	ifort -I../src/ -O2 dielectric_water_CNT.f90 -o dielectric.x ../src/libanalysis.o
	ifort -I../src/ -O2 dielectric_water.f90 -o dielectricw.x ../src/libanalysis.o
gnu:
	gfortran -I../src/ -O2 dielectric_water_CNT.f90 -o dielectric.x ../src/libanalysis.o
	gfortran -I../src/ -O2 dielectric_water_CNT_radial_dependence.f90 -o dielectricr.x ../src/libanalysis.o
	gfortran -I../src/ -O2 dielectric_water.f90 -o dielectricw.x ../src/libanalysis.o
debug:
	ifort -I../src/ -O0 -g -traceback dielectric_water_CNT.f90 -o dielectric.x ../src/libanalysis.o
