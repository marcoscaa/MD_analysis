!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: pos
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: force
  DOUBLE PRECISION, DIMENSION(3,3)                   :: box
  DOUBLE PRECISION, DIMENSION(3,3)                   :: stress
  INTEGER                                            :: natoms, nframes, nequil
  INTEGER                                            :: natoms_scell 
  INTEGER, PARAMETER                                 :: scell=1
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
    IF (scell > 1) CALL MAKE_SUPERCELL
    if (mod(frame,1)==0) CALL PRINT_DeepMD
  END DO

  CLOSE(1);CLOSE(2);CLOSE(3);CLOSE(4)
  CLOSE(5);CLOSE(6);CLOSE(7);CLOSE(8)

END PROGRAM CONVERT

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, file_force
  CHARACTER(100)             :: index_file, file_cel
  CHARACTER(100)             :: file_stress 

  CALL getarg(1, file_pos)
  CALL getarg(2, file_force)
  CALL getarg(3, file_cel)
  CALL getarg(4, file_stress)
  CALL getarg(5, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  natoms_scell = natoms * (scell**3)

  ALLOCATE(pos(natoms_scell,3)); pos=0.d0
  ALLOCATE(force(natoms_scell,3)); force=0.d0

  CLOSE(1)

  OPEN(unit  =  1,file  =  file_pos)
  OPEN(unit  =  2,file  =  file_force)
  OPEN(unit  =  3,file  =  file_cel)
  OPEN(unit  =  4,file  =  file_stress)
  OPEN(unit  =  5,file  =  "coord.raw")
  OPEN(unit  =  6,file  =  "force.raw")
  OPEN(unit  =  7,file  =  "box.raw")
  OPEN(unit  =  8,file  =  "virial.raw")
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j

  READ(1,*)
  READ(2,*)

  DO i = 1,natoms
    READ(1, *), pos(i,:)
    READ(2, *), force(i,:)
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
  DOUBLE PRECISION           :: det

  READ(3,*)
  READ(4,*)

  DO i = 1,3
    READ(3, *), box(i,:)
    READ(4, *), stress(i,:)
  END DO

  print *, stress
  !GPa to au
  stress = stress * 3.398931486995655e-05
  !Stress to virial
  stress = stress * det( box ) 
  !Au to eV
  stress = stress * 27.211399
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
    READ(4,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_DeepMD
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, tot

  tot = natoms_scell

  WRITE(5, fmt='(*(E14.7,3X))'), ( pos(i,1), pos(i,2), pos(i,3), i=1,tot )
  WRITE(6, fmt='(*(E14.7,3X))'), ( force(i,1), force(i,2), force(i,3), i=1,tot )
  WRITE(7, fmt='(9(E14.7,3X))'), ( box(i,1), box(i,2), box(i,3), i=1,3 )
  WRITE(8, fmt='(9(E14.7,3X))'), ( stress(i,1), stress(i,2), stress(i,3), i=1,3 )

END SUBROUTINE

DOUBLE PRECISION FUNCTION det( A )
  !Determinant of a 3X3 matrix
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)  :: A(3,3)

  det = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)- &
        A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+ &
        A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)

END FUNCTION det

SUBROUTINE MAKE_SUPERCELL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: iat,iat2, a, b, c 
  DOUBLE PRECISION           :: shift(3)

  iat = 1

  DO a = 0,scell-1
    DO b = 0,scell-1
      DO c = 0,scell-1

        shift(1) = box(1,1)*a
        shift(2) = box(2,2)*b
        shift(3) = box(3,3)*c

        DO iat2 = 1,natoms
       
          pos(iat,:) = pos(iat2,:) + shift(:) 
          force(iat,:) = force(iat2,:)

          iat = iat + 1

        END DO

      END DO
    END DO
  END DO

  !Supercell
  box = box * dble(scell)
  !Supercell virial
  stress = stress * (scell**3)

END SUBROUTINE MAKE_SUPERCELL
