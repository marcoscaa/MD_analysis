!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 350, nequil = 2149
  REAL*8, PARAMETER                        :: Pi=3.1415, D_conv=0.03316
  REAL*8                                   :: bin
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(9)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: hist
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes, neq = 0

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

  !User inputs the file names
  !PRINT *, "Index File: "
  !READ(*,fmt = "(A100)"), index_file
  !PRINT *, "Trajectory (.xyz format): "
  !READ(*,fmt = "(A100)"), file_name
  !PRINT *, "Number of frames: "
  !READ(*,fmt = "(I10)"), nframes
 
  index_file = "index_cp2k"
  file_name = "001_wat-pos-1.xyz"
  nframes = 5236

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))
  ALLOCATE(hist(nbins))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)
 
  hist = 0
  bin = 0.5*box(1)/DBLE(nbins)

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: junk2

  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, *), junk2, pos(i,:)
  END DO

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, ind
  REAL*8                     :: d, Dist, z0

  neq = 0

  DO i = 1,natoms

    z0 = pos(i,3) - pos(1,3) !Only for this case

    !Selecting only Ow from the region between TiO2 layers ("bulk" water)
    IF (( ind_atom(i) == 3 ) .and. ( z0 > 9.28 ) .and. ( z0 < 15.28 ))  THEN

      neq = neq + 1 ! number of equivalent atoms to be averaged in RDF

      DO j = 1,natoms

        IF (( i .ne. j ) .and. ( ind_atom(j) == 3 )) THEN

          d = Dist(pos(i,:),pos(j,:))
         
          IF ( d < box(1)/2. ) THEN 
            !Assigning the index to the histogram
            ind = int( d/bin ) + 1 
            hist(ind) = hist(ind) + 1
          END IF

        END IF 

      END DO

    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: dV, avg_rdf

  OPEN(unit = 2,file = "RDF.dat")

  DO i = 1,nbins
    dV = (4./3.)* Pi * DBLE(i**3 - (i-1)**3) * bin**3
    avg_rdf = hist(i) / ( dV * DBLE( (nframes - nequil) * neq) )
    WRITE(2,fmt = "(F10.5, 3X, F12.7)"), bin/2 + DBLE(i-1)*bin, avg_rdf / D_conv
  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+2)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist
