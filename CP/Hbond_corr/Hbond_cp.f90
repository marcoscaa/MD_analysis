!Computes the H-bond time correlation function 

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE        :: hcorr(:,:)
  INTEGER, ALLOCATABLE                 :: hbond_matrix(:,:,:)
END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS
    CALL ASSIGN_OH_BOND
    CALL MAKE_HBOND_MATRIX (frame) 
    CALL COARSE_GRAIN_POS (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL CORR_FUNCT 
  CALL PRINT_RESULTS 

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : hcorr, hbond_matrix
  USE parameters, ONLY : natoms, nframes, box, nlayers, &
                         OH_bond, nwater, coarse, pos
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file ) 

  ALLOCATE(pos(natoms,3))
  ALLOCATE(hcorr(nframes,nlayers));hcorr=0.d0
  ALLOCATE(OH_bond(natoms,2));OH_bond=0
  ALLOCATE(hbond_matrix(nframes,natoms,natoms));hbond_matrix=0
  ALLOCATE(coarse(nframes,nwater));coarse=1

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HBOND_MATRIX (frame)
  !Determine what pair of atoms are Hbonded or not
  USE histogram,  ONLY : hbond_matrix
  USE parameters, ONLY : natoms, nwater
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: i1, i2, ind1, ind2
  LOGICAL                    :: is_oxygen, is_water_oxygen
  LOGICAL                    :: is_hbonded

  ind1 = 0

  DO i1 = 1, natoms

    IF ( is_water_oxygen(i1) ) THEN

      ind1 = ind1 + 1
      ind2 = 0

      DO i2 = 1,natoms

        IF ( is_oxygen(i2) ) THEN

          ind2 = ind2 + 1

          IF ( is_hbonded(i1,i2) ) THEN
            hbond_matrix(frame, ind1, ind2) = 1
          ENDIF

        END IF

      END DO

    END IF

  END DO

END SUBROUTINE MAKE_HBOND_MATRIX

SUBROUTINE CORR_FUNCT 
  !Make the h-bond correlation function from the coarse grained function
  USE histogram,  ONLY : hcorr, hbond_matrix
  USE parameters, ONLY : coarse,nlayers, nwater, nframes
  IMPLICIT NONE
  DOUBLE PRECISION           :: surv_prob(nlayers)
  INTEGER                    :: t,t0,iat,layer
  INTEGER*8                  :: counter_in(nlayers), counter_all(nlayers)
  INTEGER*8                  :: hcorr_tmp(nlayers)
  LOGICAL                    :: in_layer

  !Avoiding poor statistics from time intervals larger than half 
  !the simulation length
  DO t=0,nframes-1

    hcorr(t+1,:) = 0
    surv_prob = 0.d0

    DO t0=1,nframes-t
    
      counter_in(:) = 0
      counter_all(:) = 0
      hcorr_tmp(:) = 0

!$omp parallel do private(layer) reduction(+:counter_in, counter_all, hcorr_tmp)
      DO iat=1,nwater
 
        layer = coarse(t0,iat)
        counter_all(layer) = counter_all(layer) + 1
 
        !Absorbing boundaries
        IF ( in_layer(iat,t0,t0+t) ) THEN
          hcorr_tmp(layer) = hcorr_tmp(layer) &
                           + sum(hbond_matrix(t0,iat,:) & 
                           * hbond_matrix(t0+t,iat,:))
          counter_in(layer) = counter_in(layer) + 1
        END IF

      END DO
!$omp end parallel do
  
      DO layer = 1,nlayers
        IF ( counter_all(layer) /= 0 ) THEN
          surv_prob(layer) = surv_prob(layer) &
                           + dble(counter_in(layer)) / dble(counter_all(layer))
          hcorr(t+1,layer) = hcorr(t+1,layer) & 
                       + dble(hcorr_tmp(layer)) / dble( counter_all(layer) )
        END IF
      END DO

    END DO !t0

    surv_prob = surv_prob / dble(t0-1) 

    DO layer = 1, nlayers
      IF ( surv_prob(layer) > 0.d0 ) THEN
        hcorr(t+1,layer) = hcorr(t+1,layer) / dble( t0-1 ) / surv_prob(layer)
      END IF
    END DO

  END DO !t

END SUBROUTINE CORR_FUNCT

SUBROUTINE PRINT_RESULTS 
  USE parameters, ONLY : dt, nframes
  USE histogram,  ONLY : hcorr
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: stime 

  OPEN(unit = 2,file = "hcorr_layer.dat")

  stime = dt*2.418884326509D-2 !fs

  DO i = 1,nframes

    WRITE(2,fmt = "(E12.5, *(3X, E15.7))"), dble(i-1)*stime, hcorr(i,:)/hcorr(1,:) 

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
