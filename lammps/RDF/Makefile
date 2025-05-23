all:
	ifort -I../src/ -O2 RDF_decomposed.f90 -o rdf_dec.x ../src/libanalysis.o
	ifort -O2 RDF.f90 -o rdf.x 
	ifort -O2 RDF_number_density.f90 -o rdf_nd.x 
	ifort -O2 RDF_xy.f90 -o rdfxy.x 
gnu:
	gfortran -O2 RDF.f90 -o rdf.x 
	gfortran -O2 RDF_water_ion.f90 -o rdf_wi.x 
	gfortran -O2 RDF_number_density.f90 -o rdf_nd.x 
	gfortran -O2 RDF_oligoamines.f90 -o rdf_oa.x
debug:
	ifort -O0 -g -traceback -check all RDF.f90 -o rdf.x
	ifort -O0 -g -traceback -check all RDF_tio2.f90 -o rdf_tio2.x
