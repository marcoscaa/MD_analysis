!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: rotcorr(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: vect(:,:,:)
  INTEGER                                  :: stride, tcorr
  INTEGER                                  :: typerot1,typerot2
  INTEGER                                  :: nattype1,nattype2
  INTEGER, ALLOCATABLE                     :: is_bound(:,:,:), rindex(:,:)
  INTEGER                                  :: nhist, central_atom

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
                         nattype1, nattype2, typerot1, typerot2, &
                         nhist, central_atom, rindex
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
  READ(1,*) typerot1, typerot2, tcorr, nhist, central_atom
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(vect(3,nframes,nattype2)); vect=0 ! normal bond vector
  ALLOCATE(is_bound(nframes,nattype1,nattype2)); is_bound=0
  ALLOCATE(rindex(nattype1,nframes)); rindex=0
  ALLOCATE(rotcorr(nhist,tcorr)); rotcorr = 0
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE ROT_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : rotcorr, vect, stride, tcorr, nattype1, &
                         nattype2, is_bound, rindex, nhist
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  DOUBLE PRECISION           :: dotprod(nhist)
  INTEGER                    :: frame1, frame2, i_at, j_at, i
  INTEGER                    :: c(nhist), ind

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    c=0
    dotprod = 0.D0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

!$omp parallel do private(ind) reduction(+:dotprod,c)
      DO i_at = 1,nattype1

        ind = rindex(i_at,frame2)

        if (ind > 0) then

          DO j_at = 1,nattype2

            IF (ALL(is_bound(frame2:frame2+frame1,i_at,j_at)==1)) THEN

              dotprod(ind) = dotprod(ind) + &
                             sum(vect(:,frame2,j_at)*vect(:,frame2+frame1,j_at))
              c(ind) = c(ind) + 1

            END IF

          END DO

        end if

      ENDDO
!$omp end parallel do

    ENDDO

    DO i=1,nhist
      IF (c(i)>0) THEN
        rotcorr(i,frame1+1) = dotprod(i) / dble( c(i) ) 
      END IF
    END DO

  ENDDO

END SUBROUTINE ROT_CORR

SUBROUTINE ASSIGN_VECTOR ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : vect, typerot1, typerot2, is_bound, rindex, &
                         central_atom, nhist
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, jat, iw, jw
  INTEGER                    :: assign_hist_index
  REAL*8                     :: vec(3), d, Dist
  REAL*8, PARAMETER          :: rcut=1.4d0

  iw = 1
  DO iat=1,natoms

    IF ( atype(iat)==typerot1 ) THEN 

      jw=1
      rindex(iw,frame) = assign_hist_index(Dist(iat,central_atom),nhist)

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

INTEGER FUNCTION assign_hist_index(d,nhist)
  IMPLICIT NONE
  REAL*8,INTENT(IN) :: d
  INTEGER, INTENT(IN) :: nhist
  REAL*8, PARAMETER :: dmax=6.0d0

  if (d <= dmax) then
    assign_hist_index = int(float(nhist)*d/dmax) + 1
  else
    assign_hist_index = 0 
  end if

END FUNCTION assign_hist_index

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

      write(2,fmt = '(E11.4,3X,*(E12.5,3X))'), dble(i)*tconst, rotcorr(:,i+1)

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
