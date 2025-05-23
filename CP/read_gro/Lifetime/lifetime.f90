!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: hist_lt(:,:)
  DOUBLE PRECISION                         :: box(3)
  INTEGER*8, ALLOCATABLE                   :: coarse(:,:,:)
  INTEGER,PARAMETER                        :: nframes = 1000, nbins = 10
  REAL*8,PARAMETER                         :: dt = 2. ! in ps
  INTEGER                                  :: natoms, natO

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

  CALL GET_LIFETIME
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
  ALLOCATE(hist_lt(nbins,nframes)); hist_lt = 0

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
  coarse(frame,:,:) = 0

  !DO j = 1, natoms
  DO j = 1, natO

    DO i = 1, nbins

      z = z_pbc( pos(3,j,frame), box(3) )
      lower = DBLE(i-1) * box(3) / dble(nbins)
      upper = DBLE(i) * box(3) / dble(nbins)

      !If molecule is inside the window, set value to 1
      IF ( ( z >= lower ) .and. ( z < upper ) ) THEN
        coarse(frame,i,j) = 1
        EXIT
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE COARSE_GRAIN 

SUBROUTINE GET_LIFETIME
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame1, frame2, i_at, i_bin
  INTEGER                    :: counter, counter_tmp 
  DOUBLE PRECISION           :: sum_hist

  DO i_bin = 1,nbins
  
    !frame1 is the duration of the interval 
    DO frame1 = 1,nframes

      counter = 0

      !Change initial condition, look for combinations with same interval
      DO frame2 = 1,nframes-frame1

        counter_tmp = 0

!$omp parallel do reduction(+:counter_tmp)
        DO i_at = 1,natO

          !If the particle did not leave the window during the interval frame1
          IF ( ALL( coarse(frame2:frame1+frame2,i_bin,i_at) == 1 ) ) THEN

            counter_tmp = counter_tmp + 1

          ENDIF
        ENDDO
!$omp end parallel do

        counter = counter + counter_tmp

      ENDDO

      !Probability of finding lifitime of duration frame1
      hist_lt(i_bin,frame1) = dble(counter) / dble( nframes-frame1 ) 

    ENDDO

    !Generating a normalized probability
    hist_lt(i_bin,:) = hist_lt(i_bin,:) / sum( hist_lt(i_bin,1:nframes-1) ) 

  ENDDO

END SUBROUTINE GET_LIFETIME

SUBROUTINE PRINT_RESULTS (mode)
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j
  CHARACTER(10),INTENT(IN)   :: mode

  !Write the lifetime for Ow separated in bins
  write(*, fmt = '(a2,3x,<nbins>(f8.4,3x))'), "# ",&
            &( ( dble(i)-0.5 ) * box(3) / dble(nbins), i=1,nbins) 

  DO i = 1,nframes-1

    write(*,fmt = '(f8.3,3X,<nbins>(f10.6,3X))'), float(i)*dt, hist_lt(:,i) 

  ENDDO

  OPEN( unit = 2, file = "lifetime.dat" )

  DO i = 1, nbins

    write(2,fmt = '(f8.4,3x,f8.4)') (dble(i)-.5) * box(3) / dble(nbins), &
         & sum( (/ ( hist_lt(i,j) * dble(j), j=1,nframes-1) /) )

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

DOUBLE PRECISION FUNCTION sum_hist( hist, len )
  !Idiotic sum of an array
  IMPLICIT NONE
  INTEGER                    :: i, len
  DOUBLE PRECISION           :: sum_out
  DOUBLE PRECISION           :: hist(1,len) 

  sum_out = 0.

  DO i=1,len

    sum_out = sum_out + hist(1,i)

  ENDDO

  sum_hist = sum_out
  
END FUNCTION sum_hist
