
MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: velt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS_VEL 
    CALL READ_CEL
    !CALL SET_CENTER_OF_MASS_TO_ZERO
    !CALL SET_CENTER_OF_MASS_VELOCITY_TO_ZERO
    CALL IDENTIFY_OH_GROUPS
    CALL ASSIGN_VEL (frame) 
    CALL COARSE_GRAIN_POS (frame)
    CALL COARSE_GRAIN_HYDROGEN (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL VEL_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : velt, autocorr 
  USE parameters, ONLY : pos, vel, natoms, nwater, &
                         nframes, coarse, nlayers,&
                         OH_index
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, file_vel, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, file_vel)
  CALL getarg(3, index_file)

  CALL READ_INDEX( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(vel(natoms,3)) ! velocities of water atoms
  ALLOCATE(velt(nframes,3*nwater,3)) ! velocities of water atoms
  ALLOCATE(coarse(nframes,3*nwater)); coarse = 1
  ALLOCATE(autocorr(nframes,nlayers)); autocorr = 0
  ALLOCATE(OH_index(natoms)); OH_index = 0

  OPEN(unit = 1,file = file_pos)
  OPEN(unit = 2,file = file_vel)
  
END SUBROUTINE INITIALIZE

SUBROUTINE COARSE_GRAIN_HYDROGEN (frame)
  !Coarse grained position of Hydrogen is assigned to the 
  !same coarse grained position of their oxygens
  USE parameters, ONLY : nwater, coarse, OH_index, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat2, iat, ih, iw 
  INTEGER                    :: ind_layer
  LOGICAL                    :: is_hydrogen, is_oxygen

  ih = nwater 

  DO iat = 1, natoms

    IF ( is_hydrogen(iat) ) THEN

      ih = ih+1
      iw = 0
      
      DO iat2 = 1, natoms

        IF ( is_oxygen(iat2) ) THEN

          iw = iw + 1

          IF ( OH_index(iat).eq.iat2 ) THEN
            coarse(frame,ih) = ind_layer(iat2) 
            GO TO 11
          END IF

        END IF

      END DO

      PRINT *, 'This H is not bound ', iat
      STOP

11  END IF

  END DO

END SUBROUTINE COARSE_GRAIN_HYDROGEN

SUBROUTINE ASSIGN_VEL( frame )
  !Keep velt vector only with OW and HW atoms
  USE parameters, ONLY : vel, natoms
  USE histogram,  ONLY : velt
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  LOGICAL                    :: is_hydrogen
  LOGICAL                    :: is_water_oxygen

  iw = 1

  DO iat=1,natoms
    IF ( is_water_oxygen(iat) ) THEN
      velt(frame,iw,:) = vel(iat,:)
      iw = iw+1
    END IF
  END DO

  DO iat=1,natoms
    IF ( is_hydrogen(iat) ) THEN
      velt(frame,iw,:) = vel(iat,:)
      iw = iw+1
    END IF
  END DO

END SUBROUTINE ASSIGN_VEL

SUBROUTINE VEL_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : velt, autocorr
  USE parameters, ONLY : nframes, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_vel(nlayers)
  DOUBLE PRECISION           :: dotprod, surv_prob(nlayers)
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: layer, i_bin
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1

    surv_prob = 0.d0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,20

      counter_in = 0
      counter_all = 0
      sum_vel = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, counter_all, sum_vel)
      DO i_at = 1,3*nwater

        layer = coarse(frame2,i_at)
        counter_all(layer) = counter_all(layer) + 1

        IF ( in_layer(i_at,frame2,frame2+frame1) ) THEN

          counter_in(layer) = counter_in(layer) + 1
          dotprod = 0.d0

          DO i_coord = 1,3
            dotprod = dotprod + &
                    & velt(frame2, i_at, i_coord) * &
                    & velt(frame2+frame1, i_at, i_coord)
          END DO

          sum_vel(layer) = sum_vel(layer) + dotprod

        ENDIF
      ENDDO
!$omp end parallel do

      DO i_bin = 1,nlayers
        IF ( counter_all(i_bin) /= 0 ) THEN
          surv_prob(i_bin) = surv_prob(i_bin) &
                           + dble(counter_in(i_bin)) / dble(counter_all(i_bin))
          !Averaged correlation function
          autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) + &
                                sum_vel(i_bin) / dble( counter_all(i_bin) ) 
        END IF
      END DO

    ENDDO

    surv_prob = surv_prob / dble(frame2-1) 

    DO i_bin = 1, nlayers

      IF ( surv_prob(i_bin) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(frame1+1,i_bin) = autocorr(frame1+1,i_bin) / dble( frame2-1 ) / surv_prob(i_bin) 
      END IF

    END DO

  ENDDO

END SUBROUTINE VEL_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,   ONLY : autocorr
  USE parameters, ONLY : box, dt, nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2

  OPEN(unit = 3,file = 'vac.dat')

  DO i = 0,nframes-1

    write(3,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*dt*au_to_fs, autocorr(i+1,:) / autocorr(1,:) 

  ENDDO

  CLOSE(3)

END SUBROUTINE PRINT_RESULTS

