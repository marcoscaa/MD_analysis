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
  REAL*8                             :: volume, T, dipole_moment
END MODULE wannmod

PROGRAM Dielectric 
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    CALL COMPUTE_VOLUME 
    CALL COMPUTE_DIPOLE 
    CALL COMPUTE_DIPOLE_VARIANCE
    IF (MOD(frame,10)==0) CALL PRINT_RESULTS(frame)
  END DO

  CALL PRINT_RESULTS(frame)

  CLOSE(1);CLOSE(3)

END PROGRAM Dielectric

SUBROUTINE INITIALIZE
  USE parameters, only : natoms,nframes,nwann, pos,&
                  wannier, atype 
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
  OPEN(unit = 3,file = wan_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_DIPOLE
  USE parameters, only : pos, atype, natoms, wannier, nwann
  USE wannmod, only : atypeO, atypeH, dipole, dipole_moment
  IMPLICIT NONE
  INTEGER :: iat, iat2,c
  REAL*8 :: d(4), diptmp(3)

  dipole=0.d0
  dipole_moment=0.d0
  c=0

  DO iat=1,natoms
    if (atype(iat).eq.atypeO) then
      diptmp=0.d0
      c=c+1
      DO iat2=1,natoms
        if (atype(iat2).eq.atypeH) then
          CALL DISTANCE_VECTOR_IND(iat2,iat,d)
          if (d(4)<1.3d0) then
            diptmp=diptmp+d(1:3)
          end if
        end if
      end do
      diptmp=diptmp-8*wannier(:,iat)
      dipole_moment=dipole_moment+norm2(diptmp)
      dipole=dipole+diptmp
    end if
  end do

  dipole_moment=dipole_moment/float(c)/0.2081943 

END SUBROUTINE COMPUTE_DIPOLE

SUBROUTINE COMPUTE_DIPOLE_VARIANCE
  USE wannmod
  IMPLICIT NONE
  dipole_variance(:,1) = dipole_variance(:,1) + dipole*dipole
  dipole_variance(:,2) = dipole_variance(:,2) + dipole
END SUBROUTINE COMPUTE_DIPOLE_VARIANCE

SUBROUTINE COMPUTE_VOLUME
  USE parameters, only : box
  USE wannmod, only : volume
  IMPLICIT NONE
  volume=volume+box(1)*box(2)*box(3)
END SUBROUTINE COMPUTE_VOLUME

SUBROUTINE PRINT_RESULTS(frame)
  USE wannmod
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  REAL*8                     :: eps, N, V, dip(3,2)
  REAL*8, PARAMETER          :: pi=3.1415
  
  N = dble(frame)
  dip=dipole_variance/N
  V = volume/N
  !eps = 4.d0*pi/(3.d0*kb*T*V)*dip(1,1)
  eps = 1.d0/(3.d0*kb*T*V)*sum(dip(:,1)-dip(:,2)**2)
  eps = eps * electron_charge**2/angstrom_to_m
  eps = eps/eps0+1.d0

  WRITE(*,fmt="(4(E15.8,3X))"), eps, dip(1,1), sum(dip(:,2)**2), dipole_moment

END SUBROUTINE

