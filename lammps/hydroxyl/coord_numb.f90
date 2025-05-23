!This is the code I used to check proton transfer at the
!TiO2/water interface. 4 output files will be generated.
!With those, it is possible to compute statistics and 
!dynamics of proton transfer on TiO2. NOTE: This is not a 
!general code!!!

MODULE parameters
  IMPLICIT NONE 
  REAL*8, PARAMETER                        :: cutoff_OH=1.2, cutoffOwTi=2.6
  REAL*8, PARAMETER                        :: cutoff_TiOti=2.6 
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(3)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes, tstep, nequil
  INTEGER                                  :: nTi5c 
  INTEGER, ALLOCATABLE                     :: index_ti5c(:), index_o2c(:)

END MODULE parameters

PROGRAM Hydrox
  USE parameters, ONLY : nequil, nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL ANALYSIS(frame)
  END DO

  CLOSE(1);CLOSE(3);CLOSE(4);CLOSE(5);CLOSE(6)

END PROGRAM Hydrox

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         ind_atom
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(ind_atom(natoms))

  CLOSE(1)

  CALL Find_Ti5c_and_O2c(pos_file)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 3, file = "n_hydroxyl.dat")
  OPEN(unit = 4, file = "hydroxyl_per_site.dat")
  OPEN(unit = 5, file = "connected_hydroxyl_per_site.dat")
  OPEN(unit = 6, file = "non_connected_hydroxyl_per_site.dat")
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, box, ind_atom
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
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

SUBROUTINE ANALYSIS(step)
  USE parameters, ONLY : index_ti5c, nTi5c, index_O2c, natoms, &
                         ind_atom, cutoffOwTi
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: iat, ih, iOw, iho, ttype
  INTEGER                    :: n_Hydroxyl, ind_H(2), iOti
  INTEGER                    :: hydroxyls(nTi5c), hydroxyls_connected(nTi5c)
  REAL*8                     :: mindist, Dist, tmp(3)
 
  CALL Count_Hydroxyl(index_O2c, hydroxyls, n_hydroxyl)

  hydroxyls_connected = 0

  do iat=1, nTi5c
    if ( hydroxyls(iat) == 1 ) then
      CALL mindist_index_vector( index_O2c(iat), 2, ih, tmp ) 
      do iOw = 1, natoms
        if ( ( ind_atom(iOw)==4 ) .and. ( Dist(ih, iOw) < 2.2 ) ) then
          CALL get_H(iOw, ind_H)
          do iho = 1,2 
            do iOti = 1,natoms
              if ( (ind_atom(iOti) == 4) .and. ( Dist( ind_H(iho), iOti ) < 2.2 ) ) then  
                if ( mindist( iOti, 1 ) < cutoffOwTi ) then 
                  hydroxyls_connected(iat) = 1
                  go to 130
                end if
              end if
            end do
          end do
        end if
      end do
130 end if
  end do

  WRITE(3, fmt = "(I8,3X,I2)") step, n_Hydroxyl
  WRITE(4, fmt = "(I8,3X,*(I1,2X))") step, hydroxyls
  WRITE(5, fmt = "(I8,3X,*(I1,2X))") step, hydroxyls_connected
  WRITE(6, fmt = "(I8,3X,*(I1,2X))") step, hydroxyls - hydroxyls_connected
 
END SUBROUTINE ANALYSIS

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

SUBROUTINE Find_Ti5c_and_O2c (coord_file)
  USE parameters, ONLY : natoms, index_ti5c, index_O2c, &
                         cutoff_TiOti, ind_atom, nTi5c
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: coord_file 
  INTEGER                    :: iat, iTi, iO
  INTEGER                    :: CordNumb
  INTEGER                    :: bufferTi(natoms), bufferO2c(natoms)
  REAL(8)                    :: d, Dist
  !
  !Read 1 frame to get box and coordinates
  OPEN(unit  =  1,file  =  coord_file)
  CALL READ_FRAME
  CLOSE(1)
  !
  iO = 0
  iTi= 0
  !
  DO iat=1,natoms
    !
    !Find 5 coordinated Ti atoms
    IF ( ind_atom(iat) == 1 ) THEN
      !
      IF (CordNumb(iat,3,cutoff_TiOti)==5) THEN
        !
        iTi=iTi+1
        bufferTi(iTi)=iat
        !
      END IF
      !
    !Find 2 coordinated O atoms
    ELSE IF ( ind_atom(iat) == 3 ) THEN
      !
      IF (CordNumb(iat,1,cutoff_TiOti)==2) THEN
        !
        iO=iO+1
        bufferO2c(iO)=iat
        !
      END IF
      !
    END IF
    !
  END DO
  !
  ALLOCATE(index_ti5c(iTi)); index_Ti5c=bufferTi(1:iTi)
  ALLOCATE(index_O2c(iO));   index_O2c=bufferO2c(1:iO)
  nTi5c=iTi
  !
  IF (iTi.ne.iO) THEN
    PRINT *, "STOP!!!"
    PRINT *, "Inconsistent assignment of undercoordinated Ti and O"
    STOP
  END IF
  !
END SUBROUTINE Find_Ti5c_and_O2c

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

SUBROUTINE Count_Hydroxyl(index_O2c, hydroxyls, cn)
  USE parameters, ONLY : cutoff_OH, nTi5c
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: index_O2c(nTi5c)
  INTEGER              :: hydroxyls(nTi5c), iat, cn
  INTEGER              :: CordNumb

  cn = 0
  hydroxyls=0

  DO iat=1,nTi5c

    IF ( CordNumb(index_O2c(iat),2,cutoff_OH)==1 ) then 
      cn = cn + 1 
      hydroxyls(iat)=1
    END IF

  END DO

END SUBROUTINE Count_Hydroxyl

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

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

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

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, atype => ind_atom 
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( atype(i) == 2 ) THEN

      dist_OH = Dist(ind_O,i)

      !Select only H's within a threshold
      IF ( dist_OH < cutoff_OH ) THEN
 
        !Very basic sorting of all OH distances within the threshold
        IF ( dist_OH < mindist(1) ) THEN 

          ind_H(2) = ind_H(1)
          ind_H(1) = i
          mindist(2) = mindist(1)
          mindist(1) = dist_OH

        ELSEIF ( dist_OH < mindist(2) ) THEN
 
          ind_H(2) = i
          mindist(2) = dist_OH

        END IF
 
      END IF
 
    END IF

  END DO

  !For consitency: higher index comes first
  IF ( ind_H(1) < ind_H(2) ) THEN
    i = ind_H(2)
    ind_H(2) = ind_H(1)
    ind_H(1) = i
  END IF

END SUBROUTINE get_H
