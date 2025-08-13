
MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: velt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, pos
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  CALL REMOVE_EQUIL_VEL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_VELOCITIES
    CALL COARSE_GRAIN_POS (frame)
    CALL IDENTIFY_OH_GROUPS
    CALL COARSE_GRAIN_HYDROGEN (frame)
    CALL ASSIGN_VEL (frame) 
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL VEL_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : velt, autocorr 
  USE parameters, ONLY : pos, vel, coarse, natoms, nwater, nlayers, &
                         OH_index, nframes, layers, corr_time, dt, &
                         nequil, atype, zoffset, atypeO, atypeH
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, file_vel, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, file_vel)
  CALL getarg(3, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil, nlayers, corr_time, dt
  READ(1,*) zoffset, atypeO, atypeH
  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(atype(natoms)) 
  ALLOCATE(vel(3,natoms)) 
  ALLOCATE(velt(3,3*nwater,nframes)) ! velocities of all atoms
  ALLOCATE(coarse(nframes,3*nwater)); coarse = 1
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  ALLOCATE(OH_index(natoms)); OH_index = 0

  OPEN(unit = 1,file = file_pos)
  OPEN(unit = 2,file = file_vel)
  
END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_VEL(frame)
  USE histogram, ONLY : velt
  USE parameters, ONLY : vel, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  LOGICAL                    :: is_hydrogen
  LOGICAL                    :: is_water_oxygen
  
  iw = 1

  DO iat=1,natoms
    IF ( is_water_oxygen(iat) ) THEN
      velt(:,iw,frame) = vel(:,iat)
      iw = iw+1
    END IF
  END DO

  DO iat=1,natoms
    IF ( is_hydrogen(iat) ) THEN
      velt(:,iw,frame) = vel(:,iat)
      iw = iw+1
    END IF
  END DO

END SUBROUTINE ASSIGN_VEL

SUBROUTINE VEL_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : velt, autocorr
  USE parameters, ONLY : nframes, nwater, coarse, nlayers, corr_time
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_vel(nlayers)
  DOUBLE PRECISION           :: dotprod, surv_prob(nlayers)
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: layer, i_bin, total_initial
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers) 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,corr_time

    surv_prob = 0.d0
    total_initial = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,100

      counter_in = 0
      counter_all = 0
      sum_vel = 0.D0

!$omp parallel do private(layer,dotprod) reduction(+:counter_in, counter_all, sum_vel)
      DO i_at = 1,3*nwater

        !layer = coarse(frame2,iwater(i_at,frame2))
        layer = coarse(frame2,i_at)
        counter_all(layer) = counter_all(layer) + 1
        !print *, i_at, frame2,frame2+frame1

        IF ( in_layer(i_at,frame2,frame2+frame1,1) ) THEN

          counter_in(layer) = counter_in(layer) + 1
          dotprod = 0.d0

          DO i_coord = 1,3
            dotprod = dotprod + &
                    & velt(i_coord,i_at,frame2) * &
                    & velt(i_coord,i_at,frame2+frame1)
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
          autocorr(i_bin,frame1+1) = autocorr(i_bin,frame1+1) + &
                                sum_vel(i_bin) !/ dble( counter_all(i_bin) ) 
        END IF
      END DO

      total_initial = total_initial+1

    ENDDO

    surv_prob = surv_prob / dble(total_initial) 

    DO i_bin = 1, nlayers

      IF ( surv_prob(i_bin) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(i_bin,frame1+1) = autocorr(i_bin,frame1+1) / dble( total_initial ) / surv_prob(i_bin) 
      END IF

    END DO

  ENDDO

END SUBROUTINE VEL_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,   ONLY : autocorr
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN(unit = 3,file = 'vac.dat')

  DO i = 0,nframes-1

    write(3,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*dt, autocorr(:,i+1) / autocorr(:,1) 

  ENDDO

  CLOSE(3)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE COARSE_GRAIN_HYDROGEN (frame)
  !Coarse grained position of Hydrogen is assigned to the 
  !same coarse grained position of their oxygens
  USE parameters, ONLY : nwater, coarse, OH_index, natoms, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat2, iat, ih, iw 
  INTEGER                    :: layer_index
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
            coarse(frame,ih) = layer_index(pos(3,iat2)) 
            GO TO 11
          END IF

        END IF

      END DO

      PRINT *, 'This H is not bound ', iat
      STOP

11  END IF

  END DO

END SUBROUTINE COARSE_GRAIN_HYDROGEN

