all:
	gfortran -c -cpp libanalysis.f90
	gfortran -I./ -fopenmp -O2 VAC_from_pos.f90 -o vac.x ./libanalysis.o
	gfortran -O2 RDF.f90 -o rdf.x
	gfortran -I./ -O2 msd.f90 -o msd.x ./libanalysis.o
