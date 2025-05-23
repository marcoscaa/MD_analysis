!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: force
  DOUBLE PRECISION, DIMENSION(3,3)                   :: box
  INTEGER                                            :: natoms, nframes, nequil
  LOGICAL, PARAMETER                                 :: apply_pbc=.true.

END MODULE parameters

PROGRAM CONVERT
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL READ_CEL
    IF (apply_pbc) CALL PBC_COORD
    if (mod(frame,1)==0) CALL PRINT_DeepMD
  END DO

  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(4);CLOSE(5);CLOSE(6)

END PROGRAM CONVERT

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, file_force
  CHARACTER(100)             :: index_file, file_cel

  CALL getarg(1, file_pos)
  CALL getarg(2, file_force)
  CALL getarg(3, file_cel)
  CALL getarg(4, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(force(natoms,3))

  CLOSE(1)

  OPEN(unit  =  1,file  =  file_pos)
  OPEN(unit  =  2,file  =  file_force)
  OPEN(unit  =  3,file  =  file_cel)
  OPEN(unit  =  4,file  =  "coord.raw")
  OPEN(unit  =  5,file  =  "force.raw")
  OPEN(unit  =  6,file  =  "box.raw")
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j

  READ(1,*)
  READ(2,*)

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(i,:)
    READ(2, fmt='(3(3X,E22.14))'), force(i,:)
  END DO

  !Bohr to angstrom
  !Force to eV/Angstrom
  pos = pos * 0.529177
  force = force * 51.4219

END SUBROUTINE

SUBROUTINE READ_CEL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  READ(3,*)

  DO i = 1,3
    READ(3, fmt='(3(F14.8))'), box(i,:)
  END DO

  !Bohr to angstrom
  box = box * 0.529177

END SUBROUTINE

SUBROUTINE PBC_COORD
  USE parameters 
  IMPLICIT NONE
  INTEGER                    :: iat, ipol

  DO iat = 1,natoms
    DO ipol = 1,3

      pos(iat,ipol) = pos(iat,ipol) - nint(pos(iat,ipol)/box(ipol,ipol))*box(ipol,ipol) 
      pos(iat,ipol) = pos(iat,ipol) + box(ipol,ipol)/2.d0

    END DO
  END DO

END SUBROUTINE PBC_COORD

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
    READ(2,*)
    READ(3,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_DeepMD
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  WRITE(4, fmt='(*(E14.7,3X))'), ( pos(i,1), pos(i,2), pos(i,3), i=1,natoms )
  WRITE(5, fmt='(*(E14.7,3X))'), ( force(i,1), force(i,2), force(i,3), i=1,natoms )
  WRITE(6, fmt='(9(E14.7,3X))'), ( box(i,1), box(i,2), box(i,3), i=1,3 )

END SUBROUTINE
