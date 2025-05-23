!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: rotcorr(:)
  DOUBLE PRECISION, ALLOCATABLE            :: vect(:,:,:)
  INTEGER                                  :: stride, tcorr
  INTEGER                                  :: typerot1,typerot2
  INTEGER                                  :: nattype1,nattype2
  INTEGER, ALLOCATABLE                     :: is_bound(:,:,:)

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL ASSIGN_VECTOR (frame)
  END DO

  CALL ROT_CORR 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : vect, rotcorr, stride, tcorr, is_bound, &
                         nattype1, nattype2, typerot1, typerot2
  USE parameters, ONLY : natoms, nframes, nattype,&
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nattype1, nattype2, nframes, nequil, dt, stride
  READ(1,*) typerot1, typerot2, tcorr
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(vect(3,nframes,nattype2)); vect=0 ! normal bond vector
  ALLOCATE(is_bound(nframes,nattype1,nattype2)); is_bound=0
  ALLOCATE(rotcorr(tcorr)); rotcorr = 0
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE ROT_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : rotcorr, vect, stride, tcorr, nattype1, nattype2, is_bound
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  DOUBLE PRECISION           :: dotprod, dist2
  INTEGER                    :: frame1, frame2, i_at, j_at
  INTEGER                    :: c

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    c=0
    dotprod = 0.D0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

!$omp parallel do reduction(+:dotprod,c)
      DO i_at = 1,nattype1

        DO j_at = 1,nattype2

          IF (ALL(is_bound(frame2:frame2+frame1,i_at,j_at)==1)) THEN

            dotprod = dotprod + sum(vect(:,frame2,j_at)*vect(:,frame2+frame1,j_at))
            c = c + 1

          END IF

        END DO

      ENDDO
!$omp end parallel do

    ENDDO

    rotcorr(frame1+1) = dotprod / dble( c ) 

  ENDDO

END SUBROUTINE ROT_CORR

SUBROUTINE ASSIGN_VECTOR ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : vect, typerot1, typerot2, is_bound
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, jat, iw, jw
  REAL*8                     :: vec(3), d
  REAL*8, PARAMETER          :: rcut=1.8d0

  iw = 1
  DO iat=1,natoms

    IF ( atype(iat)==typerot1 ) THEN 

      jw=1
      DO jat=1,natoms

        IF ( atype(jat)==typerot2 ) THEN

          CALL DIST_UNIT_VECTOR(iat,jat,vec,d)

          IF (d<rcut) THEN

            vect(:,frame,jw) = vec 
            is_bound(frame,iw,jw)=1

          END IF
          jw=jw+1

        END IF

      END DO
      iw=iw+1

    END IF

  END DO

END SUBROUTINE ASSIGN_VECTOR

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : rotcorr, tcorr
  USE parameters, ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: tconst 

  tconst = dt 
  OPEN(unit = 2,file = "rotcorr.dat")

  DO i = 0,tcorr-1

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*tconst, rotcorr(i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
