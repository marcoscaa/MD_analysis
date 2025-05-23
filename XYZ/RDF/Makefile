all:
	gfortran -I../src/ -O2 RDF.f90 -o rdf.x  ../src/libanalysis.o
	gfortran -I../src/ -O2 RDF_xy.f90 -o rdfxy.x  ../src/libanalysis.o
	#gfortran -I../src/ -O2 RDF_water_ion.f90 -o rdf_wi.x  ../src/libanalysis.o
debug:
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 RDF.f90 -o rdf.x  ../src/libanalysis.o
	gfortran -g -fbacktrace -fcheck=all -I../src/ -O0 RDF_xy.f90 -o rdfxy.x  ../src/libanalysis.o
