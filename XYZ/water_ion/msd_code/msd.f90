!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE                      :: msd(:)
  REAL*8, ALLOCATABLE                      :: post(:,:)

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  CALL READ_FILE

  CALL MSD_CORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : post, msd
  USE parameters, ONLY : nframes, dt, box
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) nframes, dt, box
  CLOSE(1)

  ALLOCATE(post(3,nframes)); post=0 ! coordinates
  ALLOCATE(msd(nframes)); msd = 0

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, post
  USE parameters, ONLY : nframes, box
  IMPLICIT NONE
  REAL*8                     :: sum_disp, dist2
  INTEGER                    :: frame1, frame2, i_at, i_bin, ix
  INTEGER                    :: n0

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes-1!int(nframes/4)

    n0=0
    sum_disp = 0.D0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1

      DO ix = 1,3
        dist2 = post(ix,frame2) - post(ix, frame2+frame1)
        dist2 = dist2 - nint( dist2 / box(ix) ) * box(ix)
        dist2 = dist2*dist2

        sum_disp = sum_disp + dist2 
      END DO

      n0=n0+1

    ENDDO

    msd(frame1+1) = sum_disp / float( n0 ) 

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE READ_FILE
  USE histogram, ONLY : post
  USE parameters, ONLY : nframes
  INTEGER   :: i, junk

  DO i=1,nframes
    READ(1,*) junk, post(1,i), post(2,i), post(3,i)
  END DO

END SUBROUTINE READ_FILE

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : msd
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 2,file = "msd.dat")

  DO i = 0,nframes-1!int(nframes/4)

      write(2,fmt = '(E19.12,3X,*(E19.12,3X))'), float(i)*dt, msd(i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
