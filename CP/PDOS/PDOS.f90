!Print the density of states from the PDOS output file of CP2K
!It also averages over several trajectories, if more than one
!is available
!USAGE: pdos.x file_name 


MODULE parameters
  IMPLICIT NONE 
  REAL*8, DIMENSION(:), ALLOCATABLE         :: pdos
  REAL*8, PARAMETER                         :: H_to_eV = 27.2114
  REAL*8                                    :: bin_size
  INTEGER                                   :: norb, nbands
  INTEGER, PARAMETER                        :: nframes=1, nbins=400 ! 0.05 eV
  REAL*8, DIMENSION(2)                      :: range_E = (/ -10., 10. /) / H_to_eV

END MODULE parameters

PROGRAM PDOS
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_FRAME
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM PDOS

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: i,eof
  CHARACTER(10)                             :: ind=""
  CHARACTER(100)                            :: file_name
  CHARACTER(1000)                           :: line

  CALL getarg(1,file_name)

  !READ the number of columns in the file
  OPEN(unit  =  1,file  =  file_name)

  do i=1,2
    READ(1,*)
  end do

  READ(1, fmt = '(A500)'), line

  norb=0

  do i=1,499
    IF ( line(i:i) .eq. "." ) norb = norb+1
  end do

  !With this I disconsider the first 3 columns
  norb = norb - 2 

  !Now, get the number of bands in a frame
  nbands=-2
  DO WHILE ( (ind .NE. "1") .and. (eof == 0) )
    READ(1, *, IOSTAT=eof), ind
    nbands = nbands + 1
  END DO

  IF ( eof < 0 ) nbands = nbands + 2
 
  bin_size = (range_E(2) - range_E(1)) / DBLE(nbins)
 
  ALLOCATE(pdos(nbands))
  pdos = 0.D0

  REWIND(1) 

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  !Sum the pdos contributions for each orbital
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: i, bin, junk
  REAL                                      :: occ, E
  REAL*8, DIMENSION(norb)                   :: orbitals 

  !Do nothing with the header lines
  do i =1,2
    READ(1,*)
  end do

  DO i = 1,nbands
    READ(1, *), junk, E, occ, orbitals

    IF ( ( E >= range_E(1) ) .AND. ( E <= range_E(2) ) ) THEN
      bin = int( (E - range_E(1)) / bin_size ) 
      pdos(bin) = pdos(bin) + sum(orbitals)
    END IF

  END DO

END SUBROUTINE READ_FRAME

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: i

  OPEN(unit = 2,file = "pdos.dat")

  DO i = 1,nbins
    WRITE(2,fmt = "(F10.5, 3X, F12.7)"), (bin_size*DBLE(i) + range_E(1)) * H_to_eV , pdos(i) / DBLE(nframes)
  END DO

  CLOSE(2)

END SUBROUTINE
