!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: msd(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: post(:,:,:)

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM
    CALL ASSIGN_POS (frame)
    CALL COARSE_GRAIN_POS (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL MSD_CORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : post, msd
  USE parameters, ONLY : natoms, nframes, nwater,&
                         pos, coarse, nlayers, nequil,&
                         nhist, nskip, atype, layers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nskip, nhist, nlayers
  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(post(nframes,nwater,3)); post=0 ! coordinates
  ALLOCATE(coarse(nframes,nwater)); coarse = 1
  ALLOCATE(msd(nframes,nlayers)); msd = 0
  ALLOCATE(atype(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)") atype(i)  
  END DO
  CLOSE(1)

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, post
  USE parameters, ONLY : nframes, box, nwater, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_disp(nlayers), dist2
  DOUBLE PRECISION           :: surv_prob(nlayers) 
  INTEGER                    :: frame1, frame2, i_at, i_bin, ix
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), layer 
  LOGICAL                    :: in_layer

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1

    surv_prob = 0.d0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1

      counter_in = 0
      counter_all = 0
      sum_disp = 0.D0

!$omp parallel do private(layer,dist2) reduction(+:counter_in,sum_disp,counter_all)
      DO i_at = 1,nwater

        layer = coarse(frame2,i_at)
        counter_all(layer) = counter_all(layer) + 1

        IF ( in_layer(i_at,frame2,frame2+frame1) ) THEN

          counter_in(layer) = counter_in(layer) + 1

          !Consider only the displacement in XY, no PBC
          DO ix = 1,2
            dist2 = post(frame2, i_at, ix) - post(frame2+frame1, i_at, ix)
            !dist2 = dist2 - nint( dist2 / box(ix,ix) ) * box(ix,ix)
            dist2 = dist2*dist2

            sum_disp(layer) = sum_disp(layer) + dist2 
          END DO

        ENDIF
      ENDDO
!$omp end parallel do

      DO i_bin = 1,nlayers
        IF ( counter_all(i_bin) /= 0 ) THEN
          surv_prob(i_bin) = surv_prob(i_bin) &
                           + dble(counter_in(i_bin)) / dble(counter_all(i_bin))
          !Averaged correlation function
          msd(frame1+1,i_bin) = msd(frame1+1,i_bin) + &
                                sum_disp(i_bin) / dble( counter_all(i_bin) ) 
        END IF
      END DO

    ENDDO

    surv_prob = surv_prob / dble(frame2-1) 

    DO i_bin = 1, nlayers
      IF ( surv_prob(i_bin) > 0.d0 ) THEN
        msd(frame1+1,i_bin) = msd(frame1+1,i_bin) / dble( frame2-1 ) / surv_prob(i_bin)
      END IF
    END DO

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE ASSIGN_POS ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : post
  USE parameters, ONLY : pos, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  LOGICAL                    :: is_water_oxygen

  iw = 1
 
  DO iat=1,natoms

    IF ( is_water_oxygen(iat) ) THEN
      post(frame,iw,:) = pos(iat,:)
      iw = iw +1
    END IF

  END DO

END SUBROUTINE ASSIGN_POS

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : msd
  USE parameters, ONLY : dt, box, nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: i, k
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION           :: tconst 

  tconst = au_to_fs * dt 
  OPEN(unit = 2,file = "msd_xy.dat")

  DO i = 0,nframes-1

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*tconst, &
           ( msd(i+1,k), k=1,nlayers ) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
