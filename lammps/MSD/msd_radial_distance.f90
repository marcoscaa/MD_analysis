!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: msd(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: post(:,:,:)
  INTEGER                                  :: stride, typmsd, nhist,tcorr
  INTEGER                                  :: central_atom
  INTEGER, ALLOCATABLE                     :: rindex(:,:)

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
  USE histogram,  ONLY : post, msd, stride, typmsd, central_atom, &
                         rindex, tcorr, nhist
  USE parameters, ONLY : natoms, nframes, nattype,&
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nattype, nframes, nequil, dt
  READ(1,*) stride, typmsd, nhist, tcorr
  READ(*,*) central_atom
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(post(3,nattype,nframes)); post=0 ! coordinates
  ALLOCATE(msd(nhist,tcorr)); msd = 0
  ALLOCATE(atype(natoms))
  ALLOCATE(rindex(nattype,nframes)); rindex=0

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, post, stride, tcorr, nhist, rindex
  USE parameters, ONLY : nframes, box, nattype
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_disp(nhist), dist2
  INTEGER                    :: frame1, frame2, i_at, ix
  INTEGER                    :: ind, ninlayer(nhist), ihist

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    sum_disp = 0.D0
    ninlayer=0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

      DO i_at = 1,nattype

          ind=rindex(i_at,frame2)

          if (ind <= nhist) then

            DO ix = 1,3
              dist2 = post(ix,i_at,frame2) - post(ix, i_at, frame2+frame1)
              dist2 = dist2 - nint( dist2 / box(ix) ) * box(ix)
              dist2 = dist2*dist2

              sum_disp(ind) = sum_disp(ind) + dist2 
            END DO

            ninlayer(ind)=ninlayer(ind)+1

          end if

      ENDDO

    ENDDO

    DO ihist=1,nhist
      if (ninlayer(ihist).gt.0) then
        msd(ihist,frame1+1) = sum_disp(ihist) / dble(ninlayer(ihist))
      end if
    END DO

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE ASSIGN_POS ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : post, typmsd, central_atom, rindex, nhist
  USE parameters, ONLY : pos, natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  INTEGER                    :: assign_hist_index
  REAL*8                     :: Dist

  iw = 1
 
  DO iat=1,natoms

    IF ( atype(iat)==typmsd ) THEN 
      post(:,iw,frame) = pos(:,iat)
      rindex(iw,frame) = assign_hist_index(Dist(iat,central_atom),nhist)
      iw = iw +1
    END IF

  END DO

END SUBROUTINE ASSIGN_POS

INTEGER FUNCTION assign_hist_index(d,nhist)
  IMPLICIT NONE
  REAL*8,INTENT(IN) :: d
  INTEGER, INTENT(IN) :: nhist
  REAL*8, PARAMETER :: dmax=6.0d0

  assign_hist_index = int(float(nhist)*d/dmax) + 1

END FUNCTION assign_hist_index

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : msd, tcorr
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 2,file = "msd.dat")

  DO i = 0,tcorr-1

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))') dble(i)*dt, msd(:,i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
