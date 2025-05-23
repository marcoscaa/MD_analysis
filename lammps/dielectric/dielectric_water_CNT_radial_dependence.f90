!Compute dielectric constant of water using dipole fluctuation

MODULE wannmod
  IMPLICIT NONE 
  INTEGER                            :: atypeO, atypeH, atypeC
  INTEGER,parameter                  :: nbins=100
  REAL*8,parameter                   :: d_cut=1.3d0 !OH bond cutoff
  REAL*8,parameter                   :: electron_charge = 1.602176e-19
  REAL*8,parameter                   :: angstrom_to_m = 1.e-10
  REAL*8,parameter                   :: kB = 1.3806e-23 ! J/K
  REAL*8,parameter                   :: eps0 = 8.854e-12 ! C^2 N^-1 m-2
  REAL*8,parameter                   :: maxR = 30.d0 ! Angstroms
  REAL*8                             :: dipole(nbins), int_dipole_r
  REAL*8                             :: dipole_int_dipole_r(nbins)
  REAL*8                             :: dV(nbins), R(nbins)
  REAL*8                             :: T, center(3), volume
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
    CALL FIND_CENTER_CNT
    CALL COMPUTE_GRID
    CALL COMPUTE_CNT_VOLUME(frame) 
    CALL COMPUTE_DIPOLE 
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2);CLOSE(3)

END PROGRAM Dielectric

SUBROUTINE INITIALIZE
  USE parameters, only : natoms,nframes,nwann, pos,&
                  wannier, atype, nequil, Pi, box 
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

  dipole = 0.0
  int_dipole_r = 0.0
  dipole_int_dipole_r = 0.0
  volume = 0.0

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'dielectric_radial_dependence.dat')
  OPEN(unit = 3,file = wan_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_GRID
  USE parameters, only : box, Pi
  USE wannmod, only : maxR, nbins, R, dV
  IMPLICIT NONE
  INTEGER :: i

  DO i=1,nbins
    R(i) = maxR*(float(i)-0.5)/float(nbins)
    dV(i) = Pi*(i**2-(i-1)**2)*(maxR/float(nbins))**2*box(3) 
  END DO

END SUBROUTINE COMPUTE_GRID

SUBROUTINE COMPUTE_DIPOLE
  USE parameters, only : pos, atype, natoms, wannier, nwann
  USE wannmod, only : atypeO, atypeH, dipole, int_dipole_r,&
                      dipole_int_dipole_r, nbins, R, dV
  IMPLICIT NONE
  INTEGER :: iat, iat2, ind, indr
  INTEGER :: Index_Histogram
  REAL*8 :: d(4), dR
  REAL*8 :: dipole_(nbins), int_dipole_r_

  dipole_ = 0.0
  int_dipole_r_ = 0.0

  DO iat=1,natoms
    if (atype(iat).eq.atypeO) then
      ind = Index_Histogram(pos(:,iat))
      DO iat2=1,natoms
        if (atype(iat2).eq.atypeH) then
          CALL DISTANCE_VECTOR_IND(iat2,iat,d)
          if (d(4)<1.3d0) then
            dipole_(ind)=dipole_(ind)+d(3)
          end if
        end if
      end do
      dipole_(ind)=dipole_(ind)-8*wannier(3,iat)
    end if
  end do

  DO indr=1,nbins
    int_dipole_r_ = int_dipole_r_ + &
        dipole_(indr) 
  END DO

  dipole_ = dipole_ / dV

  dipole = dipole + dipole_
  int_dipole_r = int_dipole_r + int_dipole_r_
  dipole_int_dipole_r = dipole_int_dipole_r + &
      dipole_ * int_dipole_r_

END SUBROUTINE COMPUTE_DIPOLE

SUBROUTINE COMPUTE_CNT_VOLUME(frame)
  !Assume the axis of the cylinder is the z axis
  USE parameters, only : pos,box,atype,natoms
  USE wannmod, only : volume, atypeC, center
  IMPLICIT NONE
  INTEGER :: iat, nC, frame
  REAL*8, parameter :: pi=3.1415
  REAL*8 :: radius, sumC(3), d(4)

  !Compute the radius
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
  USE parameters, only: nframes, box, Pi
  USE wannmod
  IMPLICIT NONE
  REAL*8                     :: eps(nbins)
  INTEGER                    :: iV
  
  dipole = dipole / float(nframes)
  int_dipole_r = int_dipole_r / float(nframes)
  dipole_int_dipole_r = dipole_int_dipole_r / float(nframes)
  eps = dipole_int_dipole_r !- dipole * int_dipole_r
  !eps = eps * 2 * Pi * box(3) / (kb*T)  
  eps = eps / (kb*T)  
  eps = eps * electron_charge**2/angstrom_to_m
  eps = eps/eps0+1.d0

  WRITE(2,*) "# R(A) eps_z // Volume = ", volume 

  DO iV=1,nbins
    WRITE(2,fmt="(F12.8,3X,F15.8)") R(iV), eps(iV) 
  END DO

END SUBROUTINE PRINT_RESULTS

SUBROUTINE FIND_CENTER_CNT
  USE parameters, only : pos, box, natoms, atype
  USE wannmod, only : atypeC, center
  IMPLICIT NONE
  REAL*8 :: sum_coord(3),d(3)
  INTEGER :: i, nc

  sum_coord=0.d0
  center=0.0
  nc=0
  DO i = 1,natoms
    IF (atype(i)==atypeC) THEN
      CALL DISTANCE_VECTOR(pos(:,i),center,d)
      sum_coord = sum_coord + d(1:3)
      nc = nc+1
    END IF
  END DO

  center = sum_coord/float(nc)

END SUBROUTINE FIND_CENTER_CNT

INTEGER FUNCTION Index_Histogram(pos)
  USE wannmod, ONLY : maxR,nbins,center
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(3)
  REAL*8                     :: frac, d(4)

  CALL DISTANCE_VECTOR(pos,center,d)
  frac=norm2(d(1:2))/maxR*float(nbins)
  IF ((int(frac)+1>nbins).or.(int(frac)<0)) THEN
    PRINT *, 'Value ', norm2(d(1:2)), ' outside histogram range'
    STOP
  END IF
  Index_Histogram = int(frac) + 1

END FUNCTION Index_Histogram
