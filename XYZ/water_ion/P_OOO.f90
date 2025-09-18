!(O-O*-O) Angle distribution. O* is water ion O atom

MODULE histogram
  IMPLICIT NONE 
  INTEGER,PARAMETER                  :: maxcoord=16 !buffer
  INTEGER, ALLOCATABLE               :: hist_angle(:) 
  REAL*8                             :: d_cut 
  LOGICAL, ALLOCATABLE               :: is_water_ion(:)
END MODULE histogram

PROGRAM angle_dist 
  USE parameters, ONLY : nframes
  USE, intrinsic :: iso_fortran_env, Only : iostat_end
  IMPLICIT NONE
  INTEGER :: frame, iostat

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  frame=1
  DO 
    CALL READ_EXTXYZ_IO (iostat)
    IF (iostat == iostat_end ) THEN
      EXIT
    END IF
    CALL FIND_WATER_ION
    CALL MAKE_HISTOGRAM
    frame=frame+1
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM angle_dist

SUBROUTINE INITIALIZE
  USE parameters, only : natoms,nframes,nequil,nhist, pos, atype
  USE histogram
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  PRINT*, "Give the O-O cutoff radius (in A units):"
  READ *, d_cut
  nhist=1000 !HARDCODED

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist_angle(nhist));hist_angle=0
  ALLOCATE(is_water_ion(natoms))

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : natoms, nhist
  USE histogram, ONLY : hist_angle, maxcoord, is_water_ion
  IMPLICIT NONE
  INTEGER                     :: iat, n1, n2
  INTEGER                     :: nbonds,ind 
  REAL*8                      :: bond_dir(3,maxcoord), angle

  DO iat = 1, natoms
    IF (is_water_ion(iat)) THEN
      CALL CONNECT_NEAREST_NEIGHBOHRS(iat,bond_dir,nbonds)
      DO n1 = 1,nbonds
        DO n2 = n1+1,nbonds
          angle = sum(bond_dir(:,n1)*bond_dir(:,n2))
          ind = int(nhist*(angle+1)/2.) + 1
          hist_angle(ind) = hist_angle(ind) + 1
        END DO
      END DO
    END IF
  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist 
  USE histogram, ONLY : hist_angle
  IMPLICIT NONE
  INTEGER                    :: ic
  REAL*8                     :: bin
  
  OPEN(unit = 2,file = "hist_angle.dat")
  write(2,*) '# ', 'cos(Angle)', 'Count'

  bin = 2./real(nhist)

  DO ic = 1, nhist
    WRITE(2,*) bin*(real(ic)-0.5) - 1.0, hist_angle(ic) 
  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE CONNECT_NEAREST_NEIGHBOHRS( iat,bonds,nbonds )
  !Compute the coordination number of atom with index iat
  USE parameters, ONLY : natoms, atype
  USE histogram, ONLY : d_cut, maxcoord
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat 
  INTEGER                    :: iat2, nbonds
  REAL*8                     :: d(3),normd,bonds(3,maxcoord)

  nbonds = 0
  bonds=0

  DO iat2 = 1, natoms
    IF ( (iat2.ne.iat).and.(trim(atype(iat2)).eq."O") ) THEN
      CALL DIST_UNIT_VECTOR(iat,iat2,d,normd) 
      IF ( normd < d_cut ) THEN 
        nbonds = nbonds + 1
        bonds(:,nbonds) = d  
      END IF
    END IF 
  END DO
 
END SUBROUTINE CONNECT_NEAREST_NEIGHBOHRS

SUBROUTINE FIND_WATER_ION
  USE parameters, ONLY : atype, natoms
  USE histogram, ONLY : is_water_ion
  IMPLICIT NONE
  INTEGER :: iat, cn_ow(natoms), ind_Ow(natoms)
  INTEGER :: closest_atom

  cn_ow=0
  ind_Ow=0
  is_water_ion=.false.
  DO iat=1,natoms
    if (trim(atype(iat))=="H") then 
      ind_Ow(iat) = closest_atom(iat,"O    ")
      cn_ow(ind_Ow(iat)) = cn_ow(ind_Ow(iat)) + 1
    END IF
  END DO
      
  DO iat=1,natoms
    if (trim(atype(iat))=="O") then !1: O index
      if(cn_Ow(iat)/=2) then
        is_water_ion(iat)=.True.
      end if
    end if
  END DO

END SUBROUTINE FIND_WATER_ION

INTEGER FUNCTION closest_atom(icenter,atyp)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: icenter
  CHARACTER(5), INTENT(IN) :: atyp
  INTEGER :: iat, iat_min
  REAL*8 :: d, min_d, Dist

  min_d = 100.
  iat_min = 0

  DO iat=1,natoms

    IF (trim(atype(iat))==trim(atyp)) THEN
      d = Dist(iat,icenter)
      IF (d < min_d) THEN
        min_d=d
        iat_min=iat
      END IF
    END IF

  END DO

  IF (iat_min==0) THEN
    PRINT *, "Could not find closest atom to atom", icenter
    PRINT *, "STOP!!!!"
    STOP
  END IF

  closest_atom = iat_min

END FUNCTION closest_atom
