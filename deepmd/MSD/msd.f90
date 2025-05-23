!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE            :: msd(:)
  INTEGER, ALLOCATABLE                     :: msdind(:)
  INTEGER                                  :: stride, typmsd, nmsd

END MODULE histogram

PROGRAM MSD
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL READ_BIN

  CALL SELECT_ATOM_TYPE

  CALL MSD_CORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM MSD

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : msd, stride, typmsd
  USE parameters, ONLY : natoms, nframes, &
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, dt, stride, typmsd
  CLOSE(1)

  ALLOCATE(msd(nframes)); msd = 0
  ALLOCATE(atype(natoms))

  OPEN(unit=1, file='type.raw')
  READ(1,*) atype
  CLOSE(1)
  
END SUBROUTINE INITIALIZE

SUBROUTINE SELECT_ATOM_TYPE
  USE histogram, only  : typmsd, msdind, nmsd
  USE parameters, only : natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iat, tmp(natoms)
  
  nmsd = 0
  tmp=0

  DO iat=1,natoms
    IF (atype(iat)==typmsd) THEN
      nmsd=nmsd+1
      tmp(nmsd)=iat
    END IF
  END DO

  ALLOCATE(msdind(nmsd))
  msdind=tmp(1:nmsd)  

  IF (ANY(msdind==0)) THEN
    WRITE(*,*) 'Something wrong with the assignment of atom types for MSD'
    WRITE(*,*) 'STOP!!!'
    STOP
  END IF

END SUBROUTINE SELECT_ATOM_TYPE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, stride, nmsd, msdind
  USE parameters, ONLY : nframes, boxt, post
  IMPLICIT NONE
  REAL*8                     :: sum_disp, dist2
  INTEGER                    :: frame1, frame2, i_at, i_bin, ix
  INTEGER                    :: n0

  !frame1 is the duration of the interval 
  DO frame1 = 0,int(nframes/4)

    n0=0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,1 !nframes-frame1,stride

      sum_disp = 0.D0

!$omp parallel do private(dist2) reduction(+:sum_disp)
      DO i_at = 1,nmsd

        DO ix = 1,3
          dist2 = post(ix,msdind(i_at),frame2) - post(ix, msdind(i_at), frame2+frame1)
          !Warning: only true for orthorhombic cells!!!
          dist2 = dist2 - nint( dist2 / boxt(ix,ix,frame2) ) * boxt(ix,ix,frame2)
          dist2 = dist2*dist2

          sum_disp = sum_disp + dist2 
        END DO

      ENDDO
!$omp end parallel do

    n0=n0+1

    ENDDO

    msd(frame1+1) = sum_disp / dble( n0 ) / dble(nmsd)

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : msd
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 2,file = "msd.dat")

  DO i = 0,int(nframes/4)

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*dt, msd(i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

