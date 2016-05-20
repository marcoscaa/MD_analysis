MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                   :: nbins = 350
  REAL*8, DIMENSION(:,:), ALLOCATABLE  :: pos
  REAL*8, DIMENSION(9)                 :: box
  INTEGER, DIMENSION(nbins)            :: hist
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ind_atom
  INTEGER                              :: natoms, nframes

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: old_i,new_i
  CHARACTER(100)             :: file_name

  PRINT *, "File Name: "
  READ(*,fmt = "(A100)"), file_name
  PRINT *, "Number of frames: "
  READ(*,fmt = "(I10)"), nframes

  !Read number of atoms and box size
  OPEN(unit = 1,file = file_name)

  READ(1,*), natoms
  READ(1,fmt = '(9(F10.4,2X))'), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms,2))

  !Make the index array
  old_i = "null"; ind_atom = 0

  DO i = 1,natoms
    READ(1, *), new_i, pos(i,:)

    IF (old_i .eq. new_i) THEN
      ind_atom(i,1) = ind_atom(i-1,1) 
    ELSE
      ind_atom(i,1) = ind_atom(i-1,1) + 1
    END IF

    !Separe TiO2 from H2O. TiO2 is labeld with 1
    IF ( pos(i,3) < 13.1 ) THEN
      ind_atom(i,2) = 1
    END IF

    old_i = new_i
    
  END DO
 
  !Restart file reading from first line
  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)

  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  REAL                       :: junk
  CHARACTER(5)               :: junk2

  READ(1,*), junk
  READ(1,*), junk

  DO i = 1,natoms
    READ(1,*), junk2, pos(i,:)
  END DO

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, ind
  REAL*8                     :: z

  DO i = 1,natoms

    !Considering only oxygen from water
    IF ((ind_atom(i,1) == 3) .and. (ind_atom(i,2) == 0)) THEN
      z = pos(i,3)  - int(pos(i,3)/box(9))*box(9)  !PBC
      ind = int( z * DBLE(nbins) / box(9) ) + 1
      hist(ind) = hist(ind) + 1
    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: bin,dV

  OPEN(unit = 2,file = "ZDF.dat")

  bin = box(9)/DBLE(nbins)
  dV = box(1)*box(5)*bin

  DO i = 1,nbins
    WRITE(2,fmt = "(F10.5, 3X, F12.7)"), bin/2 + DBLE(i-1)*bin, DBLE(hist(i))/(dV*DBLE(nframes))
  END DO

  CLOSE(2)

END SUBROUTINE



