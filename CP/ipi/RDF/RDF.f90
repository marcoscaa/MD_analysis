MODULE histogram
  IMPLICIT NONE 
  INTEGER                              :: nbins
  INTEGER, ALLOCATABLE                 :: hist(:)
  CHARACTER(5)                         :: iat1, iat2
  DOUBLE PRECISION                     :: maxd, number_density
END MODULE histogram 

PROGRAM RDF
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS_CEL
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, nbins, iat1, iat2, maxd, number_density 
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos 
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file 

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  OPEN(unit=1, file=index_file)
  READ(1,*) natoms, nframes, nequil, nbins, maxd
  READ(1,*) iat1, iat2
  CLOSE(1)

  ALLOCATE(pos(3,natoms))
  ALLOCATE(hist(nbins)); hist = 0
  ALLOCATE(atype(natoms))
  number_density=0.d0

  OPEN(unit  =  1,file  =  pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram,  ONLY : hist, iat1, iat2, nbins, maxd, number_density
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: i, j, ind, n
  DOUBLE PRECISION           :: d, Dist

  n = 0

  DO i = 1,natoms

    IF ( ( TRIM(atype(i)) == TRIM(iat1) ) .or. ( TRIM(iat1) == 'all' ) )  THEN

      DO j = 1,natoms

        IF (( i .ne. j ) .and. ( ( TRIM(atype(j)) == TRIM(iat2) ) .or. ( TRIM(iat2) == 'all' ))) THEN

            d = Dist(i,j)
            ind = int( dble(nbins) * d / maxd ) + 1 
            IF ( ind <= nbins ) hist(ind) = hist(ind) + 1
            n = n + 1

        END IF 

      END DO

    END IF

  END DO

  number_density = number_density + dble(n)/box(1,1)/box(2,2)/box(3,3)

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist, maxd, nbins, number_density
  USE parameters, ONLY : box
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: dV, avg_rdf
  DOUBLE PRECISION           :: r, dr, min_box
  DOUBLE PRECISION, PARAMETER:: Pi=3.14159265d0

  dr = maxd / dble(nbins) 
  dV = (4d0/3d0)*Pi * dr * dr * dr

  OPEN(unit = 2,file = "RDF.dat")

  DO i = 1,nbins

    r = (DBLE(i)-0.5d0)*dr
    avg_rdf = dble(hist(i)) / dV / dble( i**3 - (i-1)**3)  / number_density

    IF ( r < maxd ) THEN
      WRITE(2,fmt = "(E12.5, *(3X, E14.7))"), r, avg_rdf 
    END IF

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
