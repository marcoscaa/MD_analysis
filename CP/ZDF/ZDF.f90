!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 250
  REAL*8, PARAMETER                        :: D_conv=0.03316
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(9)                     :: box
  INTEGER, DIMENSION(:,:), ALLOCATABLE     :: hist
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes, ntype, nequil
  REAL*8                                   :: z0 !quick and dirty way to remove
                                                 !cell drift

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
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
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)
  ntype = 4

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))
  ALLOCATE(hist(ntype,nbins))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  !z0 = pos(1,3) ! System specific !!!
 
  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)
 
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters, ONLY : natoms, pos
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: junk2

  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, *), junk2, pos(i,1), pos(i,2),pos(i,3)
  END DO

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, ind
  REAL*8                     :: z!, minz

  DO i = 1,natoms

      !Applying PBC
      z = pos(i,3)  - nint(pos(i,3)/box(9))*box(9) 
      !Shit it up - all atoms with positive z
      z = z + box(9)/2.

      !Assigning the index to the histogram
      ind = int( z * float(nbins) / box(9) ) + 1
      hist(ind_atom(i),ind) = hist(ind_atom(i),ind) + 1

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+2)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: bin,dV

  OPEN(unit = 2,file = "ZDF.dat")

  bin = box(9)/DBLE(nbins)
  dV = box(1)*box(5)*bin

  DO i = 1,nbins
    WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2 + DBLE(i-1)*bin, &
                         & DBLE(hist(:,i))/(0.03316*dV*DBLE(nframes - nequil))
  END DO

  CLOSE(2)

END SUBROUTINE



