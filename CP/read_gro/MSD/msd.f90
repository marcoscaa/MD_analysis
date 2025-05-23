!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: msd
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE    :: pos
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE   :: coarse
  INTEGER, DIMENSION(:,:), ALLOCATABLE     :: hist_lt
  REAL*8, DIMENSION(3)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, natO
  INTEGER,PARAMETER                        :: nframes = 1000, nbins = 10
  INTEGER,PARAMETER                        :: nbins_hist = 50
  INTEGER*8, DIMENSION(nbins_hist)         :: histN, histO
  REAL*8, PARAMETER                        :: dt = 2. !in ps

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

  CALL GET_MSD 
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
  ALLOCATE(coarse(nframes,nbins,natO))
  ALLOCATE(msd(nframes,nbins)); msd = 0.

  histN = 0

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
      z = z_pbc( vtmp(3), box(3) )
      pos(3,io,frame) = z
      index_ = int( DBLE(nbins_hist) * z / box(3) ) + 1
      histO(index_) = histO(index_) + 1 
    
    !Making a distribution of Z coordinate for N (lipid head)
    elseif ( adjustl(trim(atomtype)).eq.'N' ) then
      z = z_pbc( vtmp(3), box(3) )
      index_ = int( DBLE(nbins_hist) * z / box(3) ) + 1
      histN(index_) = histN(index_) + 1 

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
  coarse(frame,:,:) = 0
  lower=0.34

  !DO j = 1, natoms
  DO j = 1, natO

    DO i = 1, nbins

      lower = 1. + (float(i) - 1.) * 0.66 
      upper = 1. + float(i) * 0.66 
      z = pos(3,j,frame)

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
  INTEGER                    :: i,j,k,l,counter
  REAL*8                     :: deltaX, Dist2xy
  REAL*8                     :: deltaX_tmp
  INTEGER                    :: counter_tmp

  DO i = 1,nbins
    
    DO j = 1,nframes
      deltaX = 0.
      counter = 0

      IF ( sum( coarse(1:nframes-j,i,:) ) /= 0 ) THEN !TODO: move IF

      DO k = 1, nframes-j

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

      ENDIF

    ENDDO
  ENDDO

END SUBROUTINE GET_MSD

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN( unit = 2, file='N_O_distribution.dat' )

  !Write the msd for Ow separated in bins
  write(*, fmt = '(a2,3x,<nbins>(f8.4,3x))'), "# ",&
            &( 1. + ( float(i) - 0.5 ) * 0.66, i=1,nbins) 
  DO i = 1,nframes

    !Printing the number of bond lengths
    write(*,fmt = '(f10.3,3X,<nbins>(f10.4,3X))'), dble(i)*dt, sqrt(msd(i,:)) / .27

  ENDDO

  !Write the histogram of N distribution
  DO i = 1,nbins_hist

    write(2,fmt = '(f10.4,2(3X,f10.4))'), box(3)*( DBLE(i) -.5D0 ) &
                                  & / DBLE( nbins_hist ),&
                                  & DBLE(histN(i)) / DBLE(sum(histN)),&
                                  & DBLE(histO(i)) / DBLE(sum(histO) )

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
