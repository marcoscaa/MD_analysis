!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE gofr_m
  IMPLICIT NONE 
  INTEGER                            :: ind_rdf(2)
  REAL*8                             :: maxr, number_density(2)
  REAL*8, ALLOCATABLE                :: gofr(:,:)
END MODULE gofr_m

PROGRAM RDF
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) CALL COMPUTE_RDF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, nhist, nframes, nequil, &
                         box, atype, stride
  USE gofr_m
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist
  READ(1,*) box !3 dimensions
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species

  maxr=minval(box)/2. !Max distance considered in g(r)
  number_density = 0.d0

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(gofr(2,nhist)); gofr=0.d0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_RDF
  USE parameters, ONLY : natoms, pos, box, &
                         nhist, atype
  USE gofr_m
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: ind, n_coord
  REAL*8                      :: n(2), n_tmp(2)
  REAL*8                      :: d(4), gofr_tmp(2,nhist),costheta

  n = 0.d0 !total number of pairs
  DO iat = 1, natoms

    IF ( ( atype(iat) == ind_rdf(1) ) .or. ( ind_rdf(1) == 0 ) ) THEN

      gofr_tmp=0.d0
      n_tmp=0
!$omp parallel do private(d,ind) reduction(+:gofr_tmp,n_tmp)
      DO iat2 = 1, natoms

        IF ( ( ( atype(iat2) == ind_rdf(2) ) .or. ( ind_rdf(2) == 0 ) ) &
             .and. ( iat.ne.iat2 ) ) THEN

          CALL DISTANCE_VECTOR_IND(iat,iat2,d)
          ind = int(nhist*d(4)/maxr) + 1
          costheta=d(3)/d(4)
          costheta=costheta*costheta
          IF (ind<=nhist) THEN
            gofr_tmp(1,ind) = gofr_tmp(1,ind) + (1.d0-costheta)
            gofr_tmp(2,ind) = gofr_tmp(2,ind) + costheta
          END IF
          n_tmp(1) = n_tmp(1) + 1.d0 !(1.d0-costheta)
          n_tmp(2) = n_tmp(2) + 1.d0 !costheta

        END IF

      END DO
!$omp end parallel do
      gofr = gofr + gofr_tmp
      n = n + n_tmp 

    END IF

  END DO

  number_density = number_density + n!/product(box)

END SUBROUTINE COMPUTE_RDF

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, atype, natoms, nframes
  USE gofr_m
  IMPLICIT NONE
  REAL*8                     :: bin, aver(2), vol
  REAL*8, PARAMETER          :: FourPioverThree = 4.18879
  INTEGER                    :: ihist, N 
  
  OPEN(unit = 2,file = "rdf_decomposed.dat")

  write(2,*) '#', 'R(angstrom)', 'RDF_XY', 'RDF_Z' 

  bin =  maxr/float(nhist)
  number_density = number_density !/ float(nframes)

  DO ihist = 1, nhist

    vol = FourPioverThree * ( (float(ihist)*bin)**3 - ((float(ihist)-1)*bin)**3 )
    aver = gofr(:,ihist) / number_density / vol
    WRITE(2,fmt='(3(F12.8,3X))') bin*(float(ihist)-0.5), aver(1), aver(2) 

  END DO

  CLOSE(2)

END SUBROUTINE
