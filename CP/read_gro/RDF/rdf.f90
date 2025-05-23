!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  INTEGER*8, DIMENSION(:), ALLOCATABLE     :: rdf
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE    :: pos
  INTEGER, DIMENSION(:,:), ALLOCATABLE     :: coarse
  REAL*8, DIMENSION(3)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, natO, counter
  INTEGER,PARAMETER                        :: nframes = 300, nbins = 1000
  REAL*8, PARAMETER                        :: Pi=3.1415, D_conv=33.16
  REAL*8                                   :: bin

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_FRAME (frame)
    CALL COARSE_GRAIN (frame)
  END DO

  CALL GET_RDF 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file
  CHARACTER(5)               :: moltype,atomtype
  INTEGER                    :: junk1

  ! init
  natO = 0

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)

  !Assuming constant number of particles during the simulation
  !Getting the number of atoms only once
  OPEN(unit = 1, file = file_name)
  READ(1,*)
  READ(1,*), natoms

  !get the number of O atoms 
  DO i=1,natoms
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype
    if (adjustl(trim(atomtype)).eq.'OH2') natO = natO + 1
  ENDDO
  READ(1,*), box(1), box(2), box(3)
  CLOSE(1)

  !Allocating all the arrays to be used
  ALLOCATE(pos(3,natO,nframes))
  ALLOCATE(coarse(nframes,natO))
  ALLOCATE(rdf(nbins)); rdf = 0.

  bin = 0.5*box(1)/DBLE(nbins)

  OPEN(unit = 1,file = file_name)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (frame)
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, junk1, frame, io, index_
  REAL*8                     :: vtmp(3), z_pbc, z
  CHARACTER(5)               :: moltype,atomtype

  io = 0
  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype, junk1, vtmp(:)
    if ( adjustl(trim(atomtype)).eq.'OH2' ) then
      io = io + 1
      pos(:,io,frame) = vtmp(:)
    
    endif
  END DO

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE COARSE_GRAIN (frame)
  !Create a coarse grain function to determine (with time) if an atom is
  !or not inside a specific window
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, frame
  REAL*8                     :: z, z_pbc
  REAL*8                     :: lower, upper 

  !Initializing the coarse graining function
  coarse(frame,:) = 0

  DO j = 1, natO

      lower = 1. 
      upper = 1.5
      z = z_pbc( pos(3,j,frame), box(3) )

      !If molecule is inside the window, set value to 1
      IF ( ( z >= lower ) .and. ( z < upper ) ) THEN
        coarse(frame,j) = 1
      ENDIF
  ENDDO

END SUBROUTINE COARSE_GRAIN 

SUBROUTINE GET_RDF 
  !Calculate the RDF for each 
  !of the layers in z direction
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: j,k,l, ind
  REAL*8                     :: Dist, d

  counter = 0

  DO j = 1,nframes

    DO k = 1, natO-1

      !Using absorbing boundaries, considering only particles
      !That remains in the window during k to k+j step
      IF (  coarse(j,k) == 1 ) THEN

        counter = counter + 1
        DO l = 1,natO
     
          d = Dist(pos(:,k,j),pos(:,l,j))
          IF ( d < box(1)/2. ) THEN 
            !Assigning the index to the histogram
            ind = int( d/bin ) + 1 
            rdf(ind) = rdf(ind) + 1 
          END IF

        ENDDO

      ENDIF
 
    ENDDO

  ENDDO

END SUBROUTINE GET_RDF

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: dV, avg_rdf
  CHARACTER(100)             :: fmt

  fmt = '(f10.4,3X,f10.4)'
  OPEN(unit = 2,file = "rdf.dat")

  !Write the histogram of N distribution
  DO i = 1,nbins

    dV = (4./3.)* Pi * DBLE(i**3 - (i-1)**3) * bin**3
    avg_rdf = dble(rdf(i)) / (dV * dble(counter)) 
    write(2,fmt), bin/2. + DBLE(i-1)*bin, avg_rdf / D_conv 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Dist2(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist2 = SUM(xyz*xyz)

END FUNCTION Dist2

REAL*8 FUNCTION Dist2xy(a, b)
  ! Distance between two points, in xy only, including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,2

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist2xy = xyz(1)*xyz(1) + xyz(2)*xyz(2)

END FUNCTION Dist2xy

REAL*8 FUNCTION Angle(C, a1, a2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: C, a1, a2
  REAL*8                     :: Dist2, Dist

  Angle = dot_product((C - a1),(C - a2)) / ( norm2(C - a1) * norm2(C - a2) )
  !Angle = Dist2(C,a1) + Dist2(C,a2) - Dist2(a1,a2)
  !Angle = Angle / ( 2 * Dist(C,a1) * Dist(C,a2) )

END FUNCTION Angle

REAL*8 FUNCTION z_pbc(z, box_z)
  !Apply PBC for the z coordinate olny
  IMPLICIT NONE
  REAL*8                     :: z, box_z

  IF ( z > 0 ) THEN
    z_pbc = z - int( z/box_z ) * box_z 
  ELSE
    z_pbc = z - int( -1 + z/box_z ) * box_z
  ENDIF

END FUNCTION z_pbc
