!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: msd(:)
  DOUBLE PRECISION, ALLOCATABLE            :: post(:,:,:)
  INTEGER                                  :: stride, typmsd

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL ASSIGN_POS (frame)
  END DO

  CALL MSD_CORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : post, msd, stride, typmsd
  USE parameters, ONLY : natoms, nframes, nattype,&
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nattype, nframes, nequil, dt, stride, typmsd
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(post(3,nattype,nframes)); post=0 ! coordinates
  ALLOCATE(msd(nframes)); msd = 0
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, post, stride
  USE parameters, ONLY : nframes, box, nattype
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_disp, dist2
  INTEGER                    :: frame1, frame2, i_at, i_bin, ix
  INTEGER                    :: n0

  !frame1 is the duration of the interval 
  DO frame1 = 0,int(nframes/2)

    n0=0
    sum_disp = 0.D0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

!$omp parallel do private(dist2) reduction(+:sum_disp)
      DO i_at = 1,nattype

          DO ix = 1,2
            dist2 = post(ix,i_at,frame2) - post(ix, i_at, frame2+frame1)
            !dist2 = dist2 - nint( dist2 / box(ix) ) * box(ix)
            dist2 = dist2*dist2

            sum_disp = sum_disp + dist2 
          END DO

      ENDDO
!$omp end parallel do

    n0=n0+1

    ENDDO

    msd(frame1+1) = sum_disp / dble( n0 ) / dble(nattype)

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE ASSIGN_POS ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : post, typmsd
  USE parameters, ONLY : pos, natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  LOGICAL                    :: is_water_oxygen

  iw = 1
 
  DO iat=1,natoms

    IF ( atype(iat)==typmsd ) THEN 
      post(:,iw,frame) = pos(:,iat)
      iw = iw +1
    END IF

  END DO

END SUBROUTINE ASSIGN_POS

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : msd
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: tconst 

  tconst = dt 
  OPEN(unit = 2,file = "msd.dat")

  DO i = 0,int(nframes/2)

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*tconst, msd(i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
