!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:)
  DOUBLE PRECISION                         :: box(9)
  DOUBLE PRECISION, PARAMETER              :: dt = 2.
  INTEGER, ALLOCATABLE                     :: ind_atom(:)
  INTEGER                                  :: natoms, nframes, nequil
  INTEGER                                  :: ind
  INTEGER, PARAMETER                       :: iprint = 100

END MODULE parameters

PROGRAM ZofT
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE 
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_FRAME
    CALL PRINT_Z(frame)
  END DO

  CLOSE(1);CLOSE(2)

END PROGRAM ZofT

SUBROUTINE INITIALIZE 
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  OPEN(unit = 2,file = 'zoft.dat')

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil, ind
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  nframes = nframes - nequil

  OPEN(unit  =  1,file  =  file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  READ(1,*)

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(i,1), pos(i,2), pos(i,3)
  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_Z ( frame )
  USE parameters, ONLY : pos, box, dt, iprint, ind
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  DOUBLE PRECISION           :: z, time
  DOUBLE PRECISION, PARAMETER :: au_to_ps = 2.418884326509D-5 

  z = pos(ind,3)
  z = z - nint( z/box(9) ) * box(9)
  z = z + box(9)/2.

  time = iprint * dt * au_to_ps * (frame-1)

  WRITE(2, fmt="(E12.5,3X,E12.5)") time, z

END SUBROUTINE PRINT_Z
