!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, PARAMETER                        :: bta = 0.529177
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: elup, eldwn
  DOUBLE PRECISION, DIMENSION(3,3)                   :: box
  INTEGER                                            :: nelup, neldwn, nframes, nequil

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

  !Number of electrons with spin up and down
  READ(1, *), nelup, neldwn, nframes, nequil
  READ(1, *), box(1,:), box(2,:), box(3,:)
  
  !From angstrom to bohr
  box = box / bta

  !Wannier center positions
  ALLOCATE(elup(nelup,3))
  ALLOCATE(eldwn(neldwn,3))

  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)
  OPEN(unit  =  2,file  =  'wfc.xyz')
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j

  READ(1,*)
  !READ(1,*)

  !Read positions of WFC with spin up
  DO i = 1,nelup
    READ(1, fmt='(3(3X,E22.14))'), elup(i,:)
  END DO

  READ(1,*)

  !Read positions of WFC with spin down
  DO i = 1,neldwn
    READ(1, fmt='(3(3X,E22.14))'), eldwn(i,:)
  END DO

  !Bohr to angstrom
  elup  = elup  * bta 
  eldwn = eldwn * bta

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(nelup+neldwn+2)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_XYZ
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  WRITE(2, fmt='(I5)'), nelup+neldwn
  WRITE(2, *)

  !Write WFC with spins up first
  DO i = 1, nelup
    WRITE(2, fmt='(A5,3(3X,F12.8))'), "U", elup(i,:)
  END DO
  !Write WFC with spins down
  DO i = 1, neldwn
    WRITE(2, fmt='(A5,3(3X,F12.8))'), "D", eldwn(i,:)
  END DO

END SUBROUTINE
