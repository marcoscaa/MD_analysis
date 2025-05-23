!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: msd
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE    :: pos
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE   :: coarse
  INTEGER, DIMENSION(:,:), ALLOCATABLE     :: hist_lt
  REAL*8, DIMENSION(9)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, natO, nframes, nequil
  INTEGER,PARAMETER                        :: nbins = 5
  REAL*8, PARAMETER                        :: dt = 0.001935108 !in ps

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes-nequil
    CALL READ_FRAME (frame)
    CALL COARSE_GRAIN (frame)
  END DO

  CALL GET_MSD 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  ! init
  natO = 0

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Assuming constant number of particles during the simulation
  !Getting the number of atoms only once
  OPEN(unit = 1,file = index_file)
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  box(1:3) = (/ box(1), box(5), box(9) /)

  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
    if (ind_atom(i) == 3) natO = natO + 1
  END DO

  CLOSE(1)

  !Allocating all the arrays to be used
  ALLOCATE(pos(3,natO,nframes-nequil))
  ALLOCATE(coarse(nframes-nequil,nbins,natO))
  ALLOCATE(msd(nframes-nequil,nbins)); msd = 0.

  OPEN(unit = 1,file = file_name)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (frame)
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, frame, io
  CHARACTER(5)               :: junk

  io = 0
  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    if ( ind_atom(i) == 3 ) then
      io = io + 1
      READ(1, *), junk, pos(:,io,frame)
    else
      READ(1, *)
    endif
  END DO

END SUBROUTINE

SUBROUTINE COARSE_GRAIN (frame)
  !Create a coarse grain function to determine (with time) if an atom is
  !or not inside a specific window
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, frame
  REAL*8                     :: z, z_pbc, lower, upper
  REAL*8                     :: binspace(6) 

  !Initializing the coarse graining function
  coarse(frame,:,:) = 0

  !Defining bins with respect to number density analysis
  binspace = (/ 14., 14.84, 16.51, 22.4, 24.12, 25.03 /)

  !DO j = 1, natoms
  DO j = 1, natO

    DO i = 1, nbins

      lower = binspace(i)
      upper = binspace(i+1)
      z = z_pbc( pos(3,j,frame), box(3) )

      !If molecule is inside the window, set value to 1
      IF ( ( z >= lower ) .and. ( z < upper ) ) THEN
        coarse(frame,i,j) = 1
        EXIT
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE COARSE_GRAIN 

SUBROUTINE GET_MSD 
  !Calculate the Mean Square Displacement for each 
  !of the layers in z direction
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j,k,l
  INTEGER*8                  :: counter,tframe
  REAL*8                     :: deltaX, Dist2xy
  REAL*8                     :: deltaX_tmp
  INTEGER*8                  :: counter_tmp

  tframe = nframes - nequil

  DO i = 1,nbins
    
    DO j = 1,tframe, 50
      deltaX = 0.
      counter = 0

     ! IF ( sum( coarse(1:tframe-j,i,:) ) /= 0 ) THEN !TODO: move IF

      DO k = 1, 2000, 50!tframe-j

        !Atom-specific msd already taken in account  
        deltaX_tmp =.0d0
        counter_tmp = 0
!$omp parallel do reduction(+:deltaX_tmp,counter_tmp)
        DO l = 1,natO
       
          !Using absorbing boundaries, considering only particles
          !That remains in the window during k to k+j step
          IF ( ALL( coarse(k:j+k,i,l) == 1 ) ) THEN

            deltaX_tmp = deltaX_tmp + Dist2xy(pos(:,l,k),pos(:,l,j+k))   
            counter_tmp = counter_tmp + 1
          
          ENDIF
           
        ENDDO
!$omp end parallel do
        deltaX = deltaX + deltaX_tmp
        counter = counter + counter_tmp
      ENDDO

      IF ( counter /= 0 ) THEN
        msd(j,i) = msd(j,i) + deltaX/(DBLE(counter)) 
      ENDIF

      !ENDIF

    ENDDO
  ENDDO

END SUBROUTINE GET_MSD

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  !Write the msd for Ow separated in bins
  write(*, fmt = '(a2,3x,<nbins>(f8.4,3x))'), "# ",&
            &( ( dble(i) ) * 0.67, i=1,nbins) 
  DO i = 1,nframes

    !Printing the number of bond lengths
    write(*,fmt = '(f10.3,3X,<nbins>(f10.4,3X))'), dble(i)*dt, sqrt(msd(i,:)) / 2.7

  ENDDO

END SUBROUTINE PRINT_RESULTS

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
