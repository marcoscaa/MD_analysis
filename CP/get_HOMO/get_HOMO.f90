!Print the HOMO energy of a specific atom type
!from a CP2K pdos output. Number of frames to be
!read is hardcoded in this script
!USAGE: pdos.x file_name > out.dat 


MODULE parameters
  IMPLICIT NONE 
  REAL*8                                    :: homo
  REAL*8, PARAMETER                         :: H_to_eV = 27.2114
  INTEGER                                   :: norb, nbands
  INTEGER, PARAMETER                        :: nframes=3473, nequil=1

END MODULE parameters

PROGRAM PDOS
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: frame

  CALL INITIALIZE

  !Discard equilibration
  !DO frame = 1,nequil*(nbands+2)
  !  READ(1,*)
  !END DO
  
  DO frame = 1,nframes
    CALL READ_FRAME
    !CALL PRINT_RESULTS
  END DO

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
  homo=0.

  !Header of output
  print *, '# HOMO energy (eV) | HOMO position | Fraction occupied'
 
  REWIND(1) 

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  !Sum the pdos contributions for each orbital
  USE parameters
  IMPLICIT NONE
  INTEGER                                   :: i, junk, homo_step
  REAL                                      :: occ, E, sum_orb 
  REAL*8, DIMENSION(norb)                   :: orbitals 

  occ = 1.
  orbitals=0.

  !Do nothing with the header lines
  do i =1,2
    READ(1,*)
  end do

  !Interested only in occupied orbitals
  DO i = 1, nbands 

    READ(1, *), junk, E, occ, orbitals
    !Orbitals are already sorted by energy 
    !print *, sum(orbitals), E
    IF (( sum(orbitals) > 0.5 ) .and. ( occ >= 1 )) THEN
      homo = E
      homo_step = i

    ElSEIF ( occ < 1 ) THEN
      print *, homo*H_to_eV, homo_step - i + 1, sum_orb 
    
      CALL Jump_Lines ( nbands - i )
      EXIT 

    END IF

  sum_orb = sum(orbitals)

  END DO
  
END SUBROUTINE READ_FRAME

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE

  PRINT *, homo * H_to_eV

END SUBROUTINE

SUBROUTINE Jump_Lines (nlines)
  IMPLICIT NONE
  INTEGER                                   :: nlines, i

  DO i = 1,nlines
    READ(1,*)
  END DO

END SUBROUTINE Jump_Lines
  
