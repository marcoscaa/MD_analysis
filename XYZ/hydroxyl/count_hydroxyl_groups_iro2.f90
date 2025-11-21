!This is the code I used to check proton transfer at the
!TiO2/water interface. 4 output files will be generated.
!With those, it is possible to compute statistics and 
!dynamics of proton transfer on TiO2. NOTE: This is not a 
!general code!!!

PROGRAM Hydrox
  USE parameters, ONLY : nequil, nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = nequil+1,nframes
    CALL READ_EXTXYZ
    CALL ANALYSIS(frame)
  END DO

  CLOSE(1);CLOSE(3);CLOSE(4)

END PROGRAM Hydrox

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype 
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
  ALLOCATE(atype(natoms))

  CLOSE(1)

  CALL Find_O2c(pos_file)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 3, file = "n_hydroxyl.dat")
  OPEN(unit = 4, file = "hydroxyl_per_site.dat")
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS(step)
  USE parameters, ONLY : nTi5c, index_O2c, natoms, &
                         atype, cutoffOwTi
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: iat, ih, iOw, iho, ttype
  INTEGER                    :: n_Hydroxyl, ind_H(2), iOti
  INTEGER                    :: hydroxyls(nTi5c), hydroxyls_connected(nTi5c)
  REAL*8                     :: mindist, Dist, tmp(3)
 
  CALL Count_Hydroxyl(index_O2c, hydroxyls, n_hydroxyl)

  WRITE(3, fmt = "(I8,3X,I2)") step, n_Hydroxyl
  WRITE(4, fmt = "(I8,3X,*(I1,2X))") step, hydroxyls
 
END SUBROUTINE ANALYSIS

SUBROUTINE Find_O2c (coord_file)
  USE parameters, ONLY : natoms, index_O2c, &
                         cutoff_TiO, atype, nTi5c
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: coord_file 
  INTEGER                    :: iat, iTi, iO
  INTEGER                    :: CordNumb
  INTEGER                    :: bufferTi(natoms), bufferO2c(natoms)
  REAL(8)                    :: d, Dist
  !
  !Read 1 frame to get box and coordinates
  OPEN(unit  =  1,file  =  coord_file)
  CALL READ_EXTXYZ
  CLOSE(1)
  !
  iO = 0
  iTi= 0
  !
  DO iat=1,natoms
    !
    !Find 2 coordinated O atoms
    IF ( TRIM(atype(iat)) == "O" ) THEN
      !
      IF (CordNumb(iat,"Ru   ",cutoff_TiO)==2) THEN
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
  ALLOCATE(index_O2c(iO));   index_O2c=bufferO2c(1:iO)
  nTi5c=iO
  !
  !
END SUBROUTINE Find_O2c

SUBROUTINE mindist_index_vector(ind_center, type_mindist, ind_mindist, vec_mindist)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: ind_center
  CHARACTER(5), INTENT(IN) :: type_mindist
  INTEGER                  :: ind_mindist, iat
  REAL*8                   :: vec_mindist(3)
  REAL*8                   :: dmin, Dist, d

  dmin=100.

  DO iat=1,natoms
    IF ( (TRIM(atype(iat))==TRIM(type_mindist) ) .and. (iat.ne.ind_center) ) THEN
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

  CALL DISTANCE_VECTOR_IND(ind_center,ind_mindist, vec_mindist)

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

    IF ( CordNumb(index_O2c(iat),"H    ",cutoff_OH)==1 ) then 
      cn = cn + 1 
      hydroxyls(iat)=1
    END IF

  END DO

END SUBROUTINE Count_Hydroxyl

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, atype
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( TRIM(atype(i)) == "H" ) THEN

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
