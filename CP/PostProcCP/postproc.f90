
MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:,:)
  DOUBLE PRECISION, PARAMETER              :: dt = 10.
  INTEGER                                  :: natoms, nframes, nequil
  INTEGER                                  :: ind_equil, sep, nsample
  INTEGER                                  :: iprint

END MODULE parameters

PROGRAM PPCP
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE 
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_FRAME 
  END DO

  CALL PRINT_POSTPROC

  CLOSE(1)

END PROGRAM PPCP

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
  READ(1, *), natoms, nframes, nequil, iprint, sep
  nframes = nframes - nequil

  ALLOCATE(pos(nframes,natoms,3))

  CLOSE(1)

  nsample = 0

  OPEN(unit  =  1,file  =  file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME 
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: ind_frame
  INTEGER                    :: i, ind
  REAL*8                     :: junk

  READ(1,*) ind_frame, junk

  ind_frame = ind_frame - ind_equil

  IF ( MOD( ind_frame, iprint*sep ) == 0 ) THEN

    ind = ind_frame / iprint / sep

    DO i = 1,natoms
      READ(1, fmt='(3(3X,E22.14))'), pos(ind,i,1), pos(ind,i,2), pos(ind,i,3)
    END DO

    nsample = ind 

  ELSE

    DO i = 1,natoms
      READ(1,*)
    END DO
 
  END IF

END SUBROUTINE READ_FRAME

SUBROUTINE REMOVE_EQUIL 
  USE parameters, ONLY : natoms, nequil, ind_equil
  IMPLICIT NONE
  INTEGER                    :: ind_frame
  INTEGER                    :: i,j
  REAL*8                     :: junk

  ind_frame = 0

  DO i = 1,nequil

    READ(1,*) ind_frame, junk

    DO j = 1,natoms
      READ(1,*)
    END DO

  END DO

  ind_equil = ind_frame

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_POSTPROC
  USE parameters, ONLY : pos, nsample, natoms, sep, iprint
  IMPLICIT NONE
  INTEGER                    :: i,j

  OPEN(unit = 2,file = 'postproc')

  DO i = 1,nsample

    WRITE(2, *), i*sep*iprint, '.001'

    DO j = 1,natoms

      WRITE(2, fmt='(3(3X,E22.14))'), pos(i,j,1), pos(i,j,2), pos(i,j,3)
 
    END DO

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_POSTPROC
