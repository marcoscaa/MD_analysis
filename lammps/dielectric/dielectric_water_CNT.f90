!Compute dielectric constant of water using dipole fluctuation

MODULE wannmod
  IMPLICIT NONE 
  INTEGER                            :: atypeO, atypeH, atypeC
  REAL*8,parameter                   :: d_cut=1.3d0 !OH bond cutoff
  REAL*8,parameter                   :: electron_charge = 1.602176e-19
  REAL*8,parameter                   :: angstrom_to_m = 1.e-10
  REAL*8,parameter                   :: kB = 1.3806e-23 ! J/K
  REAL*8,parameter                   :: eps0 = 8.854e-12 ! C^2 N^-1 m-2
  REAL*8                             :: dipole(3), dipole_variance(3,2)
  REAL*8                             :: volume, T
END MODULE wannmod

PROGRAM Dielectric 
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  !CALL REMOVE_EQUIL
  !CALL REMOVE_EQUIL_WANN
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    CALL COMPUTE_CNT_VOLUME(frame) 
    CALL COMPUTE_DIPOLE 
    CALL COMPUTE_DIPOLE_VARIANCE(frame)
    IF (MOD(frame,10)==0) CALL PRINT_RESULTS
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2);CLOSE(3)

END PROGRAM Dielectric

SUBROUTINE INITIALIZE
  USE parameters, only : natoms,nframes,nwann, pos,&
                  wannier, atype, nequil 
  USE wannmod
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, wan_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, wan_file)
  CALL getarg(3, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwann, nframes, T
  READ(1,*) atypeO, atypeH, atypeC
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(wannier(3,natoms)); wannier=0.d0
  ALLOCATE(atype(natoms))

  dipole_variance = 0.d0
  volume = 0.d0

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'dielectric.dat')
  OPEN(unit = 3,file = wan_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_DIPOLE
  USE parameters, only : pos, atype, natoms, wannier, nwann
  USE wannmod, only : atypeO, atypeH, dipole
  IMPLICIT NONE
  INTEGER :: iat, iat2
  REAL*8 :: d(4)

  dipole=0.d0

  DO iat=1,natoms
    if (atype(iat).eq.atypeO) then
      DO iat2=1,natoms
        if (atype(iat2).eq.atypeH) then
          CALL DISTANCE_VECTOR_IND(iat2,iat,d)
          if (d(4)<1.3d0) then
            dipole=dipole+d(1:3)
          end if
        end if
      end do
      dipole=dipole-8*wannier(:,iat)
    end if
  end do

END SUBROUTINE COMPUTE_DIPOLE

SUBROUTINE COMPUTE_DIPOLE_VARIANCE(frame)
  USE wannmod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: frame
  dipole_variance(:,1) = (dipole_variance(:,1)*dble(frame-1) + dipole*dipole)/dble(frame)
  dipole_variance(:,2) = (dipole_variance(:,2)*dble(frame-1) + dipole)/dble(frame)
END SUBROUTINE COMPUTE_DIPOLE_VARIANCE

SUBROUTINE COMPUTE_CNT_VOLUME(frame)
  !Assume the axis of the cylinder is the z axis
  USE parameters, only : pos,box,atype,natoms
  USE wannmod, only : volume, atypeC
  IMPLICIT NONE
  INTEGER :: iat, nC, frame
  REAL*8, parameter :: pi=3.1415
  REAL*8 :: radius, sumC(3), d(4), center(3)

  sumC=0.d0
  center=0.d0

  !First, find the CNT axis (center)
  nC=0
  DO iat=1,natoms
    if (atype(iat)==atypeC) then
      CALL Vector_Distance(pos(:,iat),center,d)
      sumC=sumC+d(1:3)
      nC=nC+1
    end if
  END DO

  center=sumC/dble(nC)

  !Second, compute the radius
  nC=0
  DO iat=1,natoms
    if (atype(iat)==atypeC) then
      CALL Vector_Distance(pos(:,iat),center,d)
      radius=radius+norm2(d(1:2))
      nC=nC+1
    end if
  END DO

  radius = radius / dble(nC) - 1.7d0 !Exclusion radius
  volume = (volume*dble(frame-1) + pi*radius**2*box(3))/dble(frame)

END SUBROUTINE COMPUTE_CNT_VOLUME

SUBROUTINE PRINT_RESULTS
  !Radial component is wrong
  USE wannmod
  IMPLICIT NONE
  REAL*8                     :: eps(3),eps_(3), V, dip(3,2)
  
  dip=dipole_variance
  V = volume
  !eps = 1.d0/(kb*T*V)* &
  !        (dip(:,1) - dip(:,2)**2)
  eps = dip(:,1)/(kb*T*V)
  eps = eps * electron_charge**2/angstrom_to_m
  eps = eps/eps0+1.d0

  WRITE(2,fmt="(9(E15.8,3X))"), dip(:,1), dip(:,2)**2, eps

END SUBROUTINE

