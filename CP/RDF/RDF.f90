MODULE histogram
  IMPLICIT NONE 
  INTEGER, PARAMETER                   :: nbins = 200
  INTEGER, ALLOCATABLE                 :: hist(:,:)
  INTEGER                              :: iat1, iat2
  DOUBLE PRECISION                     :: maxd
  DOUBLE PRECISION, ALLOCATABLE        :: number_density(:)
END MODULE histogram 

PROGRAM RDF
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS
    CALL READ_CEL
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(5)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, number_density, nbins, iat1, iat2, maxd 
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         ind_atom, pos, nlayers 
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, cel_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, cel_file)
  CALL getarg(3, index_file)
  READ *, iat1, iat2

  CALL READ_INDEX ( index_file )
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(hist(nlayers,nbins)); hist = 0
  ALLOCATE(number_density(nlayers)); number_density = 0

  maxd = minval( (/ (box(i,i), i=1,3) /) ) / 2.d0

  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit  =  5,file  =  cel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram,  ONLY : hist, number_density, iat1, iat2, nbins, maxd
  USE parameters, ONLY : pos, natoms, ind_atom, box, nlayers
  IMPLICIT NONE
  INTEGER                    :: i, j, ind, layer
  INTEGER                    :: layer_index, npairs(nlayers)
  DOUBLE PRECISION           :: d, Dist, z, pbc_z, det_3x3

  npairs=0

  DO i = 1,natoms

    IF ( ( ind_atom(i) == iat1 ) .or. ( iat1 == 0 ) )  THEN

      layer = layer_index(pos(i,3))

      DO j = 1,natoms

        IF (( i .ne. j ) .and. ( ( ind_atom(j) == iat2 ) .or. ( iat2 == 0 ))) THEN

            d = Dist(i,j)
            ind = int( dble(nbins) * d / maxd ) + 1 
            IF ( ind <= nbins ) hist(layer,ind) = hist(layer,ind) + 1
            npairs(layer) = npairs(layer) + 1

        END IF 

      END DO

    END IF

  END DO

  number_density = number_density + float(npairs)/det_3x3(box)

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist, maxd, nbins, number_density
  USE parameters, ONLY : nlayers, layers, box
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: dV(nlayers), avg_rdf(nlayers)
  DOUBLE PRECISION           :: r, dr, min_box
  DOUBLE PRECISION, PARAMETER:: Pi=3.14159265d0

  dr = maxd / dble(nbins) 
  min_box = min(box(1,1), box(2,2))/2.d0

  IF (allocated(layers)) THEN
  DO i = 1, nlayers
    dV(i) = Pi * dr * dr * ( layers(i+1) - layers(i) )
  END DO
  ELSE
    dV = (4d0/3d0)*Pi * dr * dr * dr
  END IF

  OPEN(unit = 2,file = "RDF.dat")

  DO i = 1,nbins

    r = (DBLE(i)-0.5d0)*dr
    avg_rdf = dble(hist(:,i)) / dV / dble( i**3 - (i-1)**3)  / number_density 

    IF ( r < min_box ) THEN
      WRITE(2,fmt = "(E12.5, *(3X, E14.7))"), r, avg_rdf 
    END IF

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
