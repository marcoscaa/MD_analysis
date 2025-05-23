!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:,:), dipol_dist(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION                         :: box(3)
  INTEGER*8, ALLOCATABLE                   :: coarse(:,:,:)
  INTEGER,PARAMETER                        :: nframes = 1000, nbins = 10
  REAL*8,PARAMETER                         :: dt = 2. ! in ps
  INTEGER                                  :: natoms, nwater

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_FRAME (frame)
    CALL COARSE_GRAIN (frame)
    CALL DIPOLE_DISTRIBUTION (frame )
  END DO

  CALL ROT_AUTOCORR 
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
  nwater = 0

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)

  !Assuming constant number of particles during the simulation
  !Getting the number of atoms only once
  OPEN(unit = 1, file = file_name)
  READ(1,*)
  READ(1,*), natoms

  !!!!!!!!get the number of O atoms !!!!!!!!!!
  DO i=1,natoms
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype
    if (adjustl(trim(atomtype)).eq.'OH2') nwater = nwater + 1
  ENDDO
  READ(1,*), box(1), box(2), box(3)
  CLOSE(1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!Allocating all the arrays to be used
  ALLOCATE(pos(3,3,nwater)) ! coordinates, atomtype (O, H1 and H2), n_water
  ALLOCATE(coarse(nframes,nbins,nwater))
  ALLOCATE(dipol_dist(3,nwater,nframes)); dipol_dist = 0.d0
  ALLOCATE(autocorr(nbins,nframes)); autocorr = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  i = 1
  DO WHILE ( i <= natoms )
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype, junk1, vtmp(:)
    if ( adjustl(trim(moltype)).eq.'TIP3' ) then
      io = io + 1
      !Oxygem coordinates
      pos(:,1,io) = vtmp(:)
      !Hydrogen coordinates
      READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype,junk1,pos(:,2,io)
      READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype,junk1,pos(:,3,io)
      i = i + 2
    endif
    i = i + 1
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
  DO j = 1, nwater

    DO i = 1, nbins

      !Considering pbc in the z direction
      z = z_pbc( pos(3,1,j), box(3) )
      !Defining upper and lower bins
      lower = 1. + 0.66 * float(i-1)
      upper = 1. + 0.66 * float(i)
      !If molecule is inside the window, set value to 1
      IF ( ( z >= lower ) .and. ( z < upper ) ) THEN
        coarse(frame,i,j) = 1
        EXIT
      ENDIF

    ENDDO
  ENDDO

END SUBROUTINE COARSE_GRAIN 

SUBROUTINE DIPOLE_DISTRIBUTION (frame)
  !Generate the dipole distribution for each water in a specific frame
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame, i

  DO i = 1,nwater

    !Calculate the dipole as:
    !dipole = ( ( O - H1 ) + ( O - H2) ) / 2.D0
    !dipol_dist (:,i,frame) = ( 2*pos(:,1,i) - pos(:,2,i) - pos(:,3,i) ) / 2.D0 
    dipol_dist (:,i,frame) = ( pos(:,1,i) - pos(:,2,i) ) / 2.D0 

    !Normalization
    dipol_dist (:,i,frame) = dipol_dist (:,i,frame) / norm2( dipol_dist(:,i,frame) )

  ENDDO

END SUBROUTINE DIPOLE_DISTRIBUTION

SUBROUTINE ROT_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE parameters
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip, sum_dip_temp
  INTEGER                    :: frame1, frame2, i_at, i_bin, i_coord
  INTEGER                    :: counter, counter_tmp 

  DO i_bin = 1,nbins
  
    !frame1 is the duration of the interval 
    DO frame1 = 1,nframes

      counter = 0
      sum_dip = 0.D0

      !Change initial condition, look for combinations with same interval
      DO frame2 = 1,10!nframes-frame1

        counter_tmp = 0
        sum_dip_temp = 0.D0

!$omp parallel do reduction(+:counter_tmp, sum_dip_temp)
        DO i_at = 1,nwater

          !If the particle dont leave the window during the interval frame1
          IF ( ALL( coarse(frame2:frame1+frame2,i_bin,i_at) == 1 ) ) THEN

            counter_tmp = counter_tmp + 1

            DO i_coord = 1,3
              sum_dip_temp = sum_dip_temp + &
                           & dipol_dist(i_coord, i_at, frame2) * &
                           & dipol_dist(i_coord, i_at, frame2+frame1)
            END DO

          ENDIF
        ENDDO
!$omp end parallel do

        counter = counter + counter_tmp
        sum_dip = sum_dip + sum_dip_temp

      ENDDO

      IF ( counter /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(i_bin,frame1) = sum_dip / ( dble( counter )  )!dble( nframes-frame1 ) )
      END IF

    ENDDO

  ENDDO

END SUBROUTINE ROT_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j

  !Write the correlation function separated in bins
  write(*, fmt = '(a2,3x,<nbins>(f8.4,3x))'), "# ",&
            &( 1. + ( float(i-1) ) * 0.66, i=1,nbins) 

  DO i = 1,nframes-1

    write(*,fmt = '(f8.3,3X,<nbins>(f10.6,3X))'), float(i)*dt,autocorr(:,i)/autocorr(:,1) 

  ENDDO

  OPEN( unit = 2, file = "lifetime.dat" )

  DO i = 1, nbins

    write(2,fmt = '(f8.4,3x,f8.4)') (dble(i)-.5) * box(3) / dble(nbins), &
         & sum( (/ ( autocorr(i,j) * dble(j), j=1,nframes-1) /) )

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
