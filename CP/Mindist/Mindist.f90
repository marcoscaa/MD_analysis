!Computes the average and instantaneous Hbond

MODULE parameters
  IMPLICIT NONE 
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(9)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes
  INTEGER                                  :: frame, nequil

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil*(natoms+2)
    READ(1, *)
  END DO
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL GET_MINDIST 
  END DO

  CLOSE(1);CLOSE(2)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  OPEN(unit = 1,file = file_name)
  OPEN(unit = 2,file = "Mindist.dat")
 

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: junk2

  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, *), junk2, pos(i,:)
  END DO

END SUBROUTINE

SUBROUTINE GET_MINDIST
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j
  REAL*8                     :: d_OOw, Dist, min_dist

  min_dist = 1000. 

  DO i = 1,natoms

    !Selecting only surface O from TiO2
    IF ( ind_atom(i) == 3 )  THEN

      DO j = 1,natoms

        !Selecting only Ow
        IF ( ind_atom(j) == 3 ) THEN

          !Distance between the two atoms
          d_OOw = Dist(pos(i,:), pos(j,:))

          IF ( d_OOw < min_dist ) min_dist = d_OOw 
 
        ENDIF ! ind_atom(j)

      END DO !j

    END IF !ind_atom

  END DO !i 

  WRITE(2,fmt = "(I10, 3X, F10.5)"), frame, min_dist 

END SUBROUTINE GET_MINDIST
 
REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist
