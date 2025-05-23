MODULE histogram
  IMPLICIT NONE 
  INTEGER                              :: nbins
  INTEGER, ALLOCATABLE                 :: hist(:),neq
  INTEGER                              :: iat1, iat2
  DOUBLE PRECISION                     :: maxd
END MODULE histogram 

PROGRAM RDF
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_RAW_POS_BOX
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, neq, nbins, iat1, iat2, maxd 
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos 
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file 

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  OPEN(unit=1, file=index_file)
  READ(1,*) natoms, nframes, nbins, maxd
  READ(1,*) iat1, iat2
  CLOSE(1)

  !Read atype from type.raw
  ALLOCATE(atype(natoms))
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  ALLOCATE(pos(3,natoms))
  ALLOCATE(hist(nbins)); hist = 0

  neq = 0

  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit  =  2,file  =  'box.raw')

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram,  ONLY : hist, neq, iat1, iat2, nbins, maxd
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: i, j, ind
  DOUBLE PRECISION           :: d, Dist

  DO i = 1,natoms

    IF ( ( atype(i) == iat1 ) .or. ( iat1 == -1 ) )  THEN

      neq = neq + 1

      DO j = 1,natoms

        IF (( i .ne. j ) .and. ( ( atype(j) == iat2 ) .or. ( iat2 == -1 ))) THEN

            d = Dist(i,j)
            ind = int( dble(nbins) * d / maxd ) + 1 
            IF ( ind <= nbins ) hist(ind) = hist(ind) + 1

        END IF 

      END DO

    END IF

  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist, maxd, nbins, neq
  USE parameters, ONLY : box
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: dV, avg_rdf
  DOUBLE PRECISION           :: r, dr, min_box
  DOUBLE PRECISION, PARAMETER:: Pi=3.14159265d0

  dr = maxd / dble(nbins) 
  min_box = minval( (/ (box(i,i), i=1,3) /) ) / 2.d0
  dV = (4d0/3d0)*Pi * dr * dr * dr

  OPEN(unit = 2,file = "RDF.dat")

  DO i = 1,nbins

    r = (DBLE(i)-0.5d0)*dr
    avg_rdf = dble(hist(i)) / dV / dble( i**3 - (i-1)**3)  / DBLE(neq)

    IF ( r < min_box ) THEN
      WRITE(2,fmt = "(E12.5, *(3X, E14.7))"), r, avg_rdf 
    END IF

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
