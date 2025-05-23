all:
	ifort -I../src/ -O2 RDF.f90 -o rdf.x ../src/libmdanalysis.o
	ifort -O2 RDF_water_ion.f90 -o rdf_wi.x
debug:
	ifort -I../src/ -O0 -traceback -check all -g RDF.f90 -o rdf.x ../src/libmdanalysis.o
	ifort -O0 -O0 -traceback -check all RDF_water_ion.f90 -o rdf_wi.x
