
PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_EXTXYZ
    if (MOD(frame,stride)==0) CALL COMPUTE_RDF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, number_density, &
                         nhist, gofr, nframes, nequil, &
                         box, ind_rdf, maxr, atype, stride
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist, maxr
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !All means all species

  number_density = 0.d0

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(gofr(nhist)); gofr=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_RDF
  USE parameters, ONLY : natoms, pos, ind_rdf, box, &
                         nhist, gofr, atype, maxr, number_density
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: n_coord 
  INTEGER                     :: ind, gofr_tmp(nhist), n_tmp
  INTEGER                     :: n
  REAL*8                      :: Dist, d, vol, M33DET

  n = 0 !total number of pairs
  DO iat = 1, natoms

    IF ( ( trim(atype(iat)) == trim(ind_rdf(1)) ) .or. ( trim(ind_rdf(1)) == "All" ) ) THEN

      gofr_tmp=0
      n_tmp=0
!$omp parallel do private(d,ind) reduction(+:gofr_tmp,n_tmp)
      DO iat2 = 1, natoms

        IF ( ( ( trim(atype(iat2)) == trim(ind_rdf(2)) ) .or. ( trim(ind_rdf(2)) == "All" ) ) &
             .and. ( iat.ne.iat2 ) ) THEN

          d = Dist(iat,iat2)
          ind = int(nhist*d/maxr) + 1
          IF (ind<=nhist) gofr_tmp(ind) = gofr_tmp(ind) + 1
          n_tmp = n_tmp + 1

        END IF

      END DO
!$omp end parallel do
      gofr = gofr + gofr_tmp
      n = n + n_tmp 

    END IF

  END DO

  vol = M33DET(box)
  number_density = number_density + float(n)/vol

END SUBROUTINE COMPUTE_RDF

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, gofr, maxr, number_density, &
                         atype, natoms, ind_rdf, nframes
  IMPLICIT NONE
  REAL*8                     :: bin, aver, vol
  REAL*8, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist, N 
  INTEGER                    :: get_natoms_type
  
  OPEN(unit = 2,file = "rdf.dat")

  write(2,*) '#', 'Z', 'RDF'

  bin =  maxr/float(nhist)
  number_density = number_density !/ float(nframes)

  DO ihist = 1, nhist

    vol = FourPioverThree * ( (float(ihist)*bin)**3 - ((float(ihist)-1)*bin)**3 )
    aver = float(gofr(ihist)) / number_density / vol
    WRITE(2,*) bin*(float(ihist)-0.5), aver 

  END DO

  CLOSE(2)

END SUBROUTINE

