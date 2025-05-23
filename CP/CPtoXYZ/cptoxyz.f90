!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: pos
  DOUBLE PRECISION, DIMENSION(3,3)                   :: box
  CHARACTER(5), ALLOCATABLE                          :: ind_atom(:)
  LOGICAL                                            :: apply_pbc
  INTEGER                                            :: natoms, nframes, nequil

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL PRINT_XYZ
    !CALL PRINT_GRO
  END DO

  CLOSE(1)
  CLOSE(2)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box(1,:), box(2,:), box(3,:)
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(A5)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  apply_pbc = .true.
  IF ( sum( box ) == 0. ) apply_pbc = .false. 

  OPEN(unit  =  1,file  =  file_name)
  OPEN(unit  =  2,file  =  'out.xyz')
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j

  READ(1,*)
  !READ(1,*)

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(i,:)
  END DO

  !Bohr to angstrom
  DO i = 1,natoms
    DO j = 1,3
      pos(i,j) = pos(i,j) * 0.529177
      IF ( apply_pbc ) THEN
        pos(i,j) = pos(i,j) - nint( ( pos(i,j)/box(j,j) ) ) * box(j,j) 
        pos(i,j) = pos(i,j) + box(j,j)/2.d0
      END IF
    END DO
  END DO

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_XYZ
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  WRITE(2, fmt='(I5)'), natoms
  WRITE(2, *)

  DO i = 1, natoms
    WRITE(2, fmt='(A5,3(3X,F12.8))'), ind_atom(i), pos(i,:)
  END DO

END SUBROUTINE

SUBROUTINE PRINT_GRO
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  WRITE(2, fmt='(a5)'), "blabl"
  WRITE(2, fmt='(i5)'), natoms

  !Angstrom to nm
  pos = pos / 10.d0
  box = box / 10.d0

  DO i = 1, natoms
    WRITE(2, fmt='(i5,2a5,i5,3f8.3)'), 0, "MOL  ", ind_atom(i), i, pos(i,:)
  END DO

  WRITE(2, fmt='(3f10.5)'), box(1,1), box(2,2), box(3,3) 

END SUBROUTINE
