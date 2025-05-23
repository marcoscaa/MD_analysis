
MODULE histogram 
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE            :: velt(:,:,:)
  REAL*8, ALLOCATABLE            :: autocorr(:,:)
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, pos
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ_POS_VEL
    CALL COARSE_GRAIN_POS (frame)
    CALL ASSIGN_VEL (frame) 
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL VEL_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : velt, autocorr 
  USE parameters, ONLY : pos, vel, coarse, natoms, nlayers, &
                         nframes, layers, corr_time, dt, &
                         nequil, atype, zoffset
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nlayers, corr_time, dt
  READ(1,*) zoffset
  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(atype(natoms)) 
  ALLOCATE(vel(3,natoms)) 
  ALLOCATE(velt(3,natoms,nframes)); velt=0.d0 ! velocities of all atoms
  ALLOCATE(coarse(nframes,natoms)); coarse = 1
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE ASSIGN_VEL(frame)
  USE histogram, ONLY : velt
  USE parameters, ONLY : vel, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat
  
  DO iat=1,natoms
    velt(:,iat,frame) = vel(:,iat)
  END DO

END SUBROUTINE ASSIGN_VEL

SUBROUTINE VEL_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : velt, autocorr
  USE parameters, ONLY : nframes, natoms, coarse, nlayers, corr_time
  IMPLICIT NONE
  REAL*8                     :: sum_vel(nlayers)
  REAL*8                     :: dotprod, surv_prob(nlayers)
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
      DO i_at = 1,natoms

        !layer = coarse(frame2,iwater(i_at,frame2))
        layer = coarse(frame2,i_at)
        IF ( layer .ne. 0 ) THEN
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
  USE parameters, ONLY : dt, corr_time
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN(unit = 3,file = 'vac.dat')

  DO i = 0,corr_time-1

    write(3,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*dt, autocorr(:,i+1) / autocorr(:,1) 

  ENDDO

  CLOSE(3)

END SUBROUTINE PRINT_RESULTS

