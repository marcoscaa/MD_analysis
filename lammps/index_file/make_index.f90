MODULE parameters
  IMPLICIT NONE 
  REAL*8, PARAMETER                        :: cutoff_OH=1.2, cutoffOwTi=2.6
  REAL*8, PARAMETER                        :: cutoff_TiOti=2.6 
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(3)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes
  INTEGER                                  :: nTi5c 
  INTEGER, ALLOCATABLE                     :: index_ti5c(:), index_o2c(:)

END MODULE parameters

PROGRAM ZDF
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE

  CALL READ_FRAME
  CALL PRINT_INDEX 

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, pos, &
                         ind_atom
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(ind_atom(natoms))

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, box, ind_atom
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  REAL*8                     :: box_tmp(2)
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, ind_atom(ind), pos(:,ind)
    pos(:,ind) = pos(:,ind) * box !Non-reduced coordinates
  END DO

END SUBROUTINE

SUBROUTINE PRINT_INDEX
  USE parameters, ONLY : natoms, ind_atom, cutoff_TiOTi
  IMPLICIT NONE
  INTEGER                    :: iat, CordNumb, cn
 
  do iat=1,natoms

    if ( (ind_atom(iat).eq.1) .or. (ind_atom(iat).eq.2) ) then

      WRITE(3, fmt = "(I1)") , ind_atom(iat)

    elseif ( ind_atom(iat).eq.3 ) then

      cn = CordNumb(iat,1,cutoff_TiOTi)
      if ( cn .ge. 2 ) then
        WRITE(3, fmt = "(I1)") , 3
      else
        if ( cn .eq. 0 ) then
          WRITE(3, fmt = "(I1)"), 4
        else
          WRITE(3, fmt = "(I1)"), 5
        end if
      end if

    end if

  end do
 
END SUBROUTINE PRINT_INDEX 

REAL*8 FUNCTION Dist(inda, indb)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) - pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

SUBROUTINE DistVect(inda, indb, vect)
  ! Vector from a to b including pbc
  USE parameters, ONLY : pos,box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz, vect
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) - pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  vect = xyz 

END SUBROUTINE DistVect

SUBROUTINE mindist_index_vector(ind_center, type_mindist, ind_mindist, vec_mindist)
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind_center, type_mindist
  INTEGER             :: ind_mindist, iat
  REAL*8              :: vec_mindist(3)
  REAL*8              :: dmin, Dist, d

  dmin=100.

  DO iat=1,natoms
    IF ( (ind_atom(iat)==type_mindist ) .and. (iat.ne.ind_center) ) THEN
      d=Dist(ind_center, iat) 
      IF ( d < dmin ) THEN
        dmin=d 
        ind_mindist=iat
      ENDIF
    ENDIF
  ENDDO

  IF (dmin>10) THEN
    print *, 'Too large mindist. Stop!'
    STOP
  ENDIF

  CALL DistVect(ind_center,ind_mindist, vec_mindist)

END SUBROUTINE mindist_index_vector

INTEGER FUNCTION CordNumb(ind,atype,cutoff)
  ! Number of nearest atoms with type atype 
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind, atype
  INTEGER                    :: iat, cn
  REAL*8, INTENT(IN)         :: cutoff !angstrom
  REAL*8                     :: Dist

  cn=0

  DO iat = 1,natoms
    IF (ind_atom(iat)==atype) THEN
      IF (Dist(ind,iat)<=cutoff) cn=cn+1
    ENDIF
  ENDDO
 
  CordNumb=cn 

END FUNCTION CordNumb

REAL*8 FUNCTION mindist(ind, type_mindist)
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind, type_mindist
  INTEGER             :: iat
  REAL*8              :: dmin, Dist, d

  dmin=100.

  DO iat=1,natoms
    IF ( (ind_atom(iat)==type_mindist ) .and. (iat.ne.ind) ) THEN
      d=Dist(ind, iat) 
      IF ( d < dmin ) THEN
        dmin=d 
      ENDIF
    ENDIF
  ENDDO

  mindist=dmin 

END FUNCTION mindist
