all:
	#ifort -O2 ZDF.f90 -o zdf.x
	ifort -I../src/ -O2 RDF.f90 -o rdf.x ../src/libmdanalysis.o
debug:
	ifort -I../src/ -O0 -traceback -check all -g RDF.f90 -o rdf.x ../src/libmdanalysis.o
