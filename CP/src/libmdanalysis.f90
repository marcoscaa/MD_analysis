
MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: vel(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: force(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: layers(:)
  DOUBLE PRECISION                         :: box(3,3)
  DOUBLE PRECISION                         :: dt !Read in the index file
  DOUBLE PRECISION, PARAMETER              :: cutoff_OH = 1.6d0 
  DOUBLE PRECISION, PARAMETER              :: pi = 3.1415926539d0 
  INTEGER                                  :: nlayers 
  INTEGER, PARAMETER                       :: iprint = 10
  INTEGER, ALLOCATABLE                     :: ind_atom(:)
  INTEGER, ALLOCATABLE                     :: n_water(:)
  INTEGER, ALLOCATABLE                     :: OH_bond(:,:)
  INTEGER, ALLOCATABLE                     :: OH_index(:)
  INTEGER*8, ALLOCATABLE                   :: coarse(:,:)
  INTEGER                                  :: natoms, ntype
  INTEGER                                  :: nwater 
  INTEGER                                  :: nOxygen 
  INTEGER                                  :: nframes
  INTEGER                                  :: nequil
  LOGICAL, PARAMETER                       :: inversion_symmetry=.false. 

END MODULE parameters

!!!!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_INDEX (index_file)
  USE parameters, ONLY : natoms, nframes, nequil, box, dt, &
                         ind_atom, nlayers, layers, nwater, &
                         nOxygen
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file
  LOGICAL                    :: is_water_oxygen, is_oxygen

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 3 lines of index file 
  READ(1, *) natoms, nframes, nequil, nlayers, dt
  READ(1, *) box(1,1), box(1,2), box(1,3), &
             box(2,1), box(2,2), box(2,3), &
             box(3,1), box(3,2), box(3,3) 

  IF(nlayers > 1) THEN
    ALLOCATE(layers(nlayers+1))
    READ(1, *) layers 
  END IF

  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)") ind_atom(i)  
  END DO
 
  !Count number of water molecules and oxygen atoms
  nwater = 0
  nOxygen = 0
  DO i = 1, natoms
    IF (is_water_oxygen(i)) nwater = nwater+1
    IF (is_oxygen(i)) nOxygen = nOxygen+1
  END DO 

  CLOSE(1)

  nframes = nframes - nequil

END SUBROUTINE READ_INDEX

SUBROUTINE READ_POS
  USE parameters, ONLY : pos, natoms
  IMPLICIT NONE
  INTEGER                    :: i

  READ(1,*)

  DO i = 1,natoms
    READ(1, *) pos(i,1), pos(i,2), pos(i,3)
  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE READ_CEL
  USE parameters, ONLY : box
  IMPLICIT NONE
  INTEGER                    :: i

  READ(5,*)

  DO i = 1,3
    READ(5, *) box(i,1), box(i,2), box(i,3)
  END DO

  !Bohr to angstrom
  box = box * 0.529177

END SUBROUTINE

SUBROUTINE READ_POS_XYZ
  USE parameters, ONLY : pos, natoms, box
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: atom_name

  box = 0
  READ(1,*)
  READ(1,*) box(1,1), box(2,2), box(3,3)

  DO i = 1,natoms
    READ(1, *) atom_name, pos(i,1), pos(i,2), pos(i,3)
  END DO

END SUBROUTINE

SUBROUTINE READ_POS_VEL
  USE parameters, ONLY : pos, vel, natoms
  IMPLICIT NONE
  INTEGER                    :: i

  READ(1,*)
  READ(2,*)

  DO i = 1,natoms
    READ(1, *) pos(i,1), pos(i,2), pos(i,3)
    READ(2, *) vel(i,1), vel(i,2), vel(i,3)
  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nequil
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE REMOVE_EQUIL_VEL
  USE parameters, ONLY : natoms, nequil
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
    READ(2,*)
  END DO

END SUBROUTINE REMOVE_EQUIL_VEL

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, ind_atom
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( ind_atom(i) == 2 ) THEN

      dist_OH = Dist(ind_O,i)

      !Select only H's within a threshold
      IF ( dist_OH < cutoff_OH+1 ) THEN
 
        !Very basic sorting of all OH distances within the threshold
        IF ( dist_OH < mindist(1) ) THEN 

          ind_H(2) = ind_H(1)
          ind_H(1) = i
          mindist(2) = mindist(1)
          mindist(1) = dist_OH

        ELSEIF ( dist_OH < mindist(2) ) THEN
 
          ind_H(2) = i
          mindist(2) = dist_OH

        END IF
 
      END IF
 
    END IF

  END DO

  !For consitency: higher index comes first
  IF ( ind_H(1) < ind_H(2) ) THEN
    i = ind_H(2)
    ind_H(2) = ind_H(1)
    ind_H(1) = i
  END IF

END SUBROUTINE get_H

SUBROUTINE ASSIGN_OH_BOND
  USE parameters, ONLY : natoms, OH_bond
  IMPLICIT NONE
  INTEGER                    :: iat, ind_H(2)
  LOGICAL                    :: is_water_oxygen
  
  DO iat = 1, natoms

    IF ( is_water_oxygen(iat) ) THEN
  
      CALL get_H(iat,ind_H)
      OH_bond(iat,:) = ind_H(:)

    END IF

  END DO
          
END SUBROUTINE ASSIGN_OH_BOND

SUBROUTINE COARSE_GRAIN_POS (frame)
  USE parameters, ONLY : natoms, coarse, box, &
                         pos, nlayers, ind_atom, nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ind_z, iwat
  INTEGER                              :: layer_index
  LOGICAL                              :: is_water_oxygen
  
  IF (nlayers < 2) RETURN

  coarse(frame,:) = 0
  iwat = 0

  DO iat = 1, natoms

    !Selecting only OW
    IF ( is_water_oxygen(iat) ) THEN

      iwat = iwat + 1
      ind_z = layer_index(pos(iat,3))
      coarse(frame,iwat) = ind_z

    END IF

  END DO

END SUBROUTINE COARSE_GRAIN_POS

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jump between layers in the coarse grain function
  USE parameters, ONLY : coarse, nwater, nframes, nlayers
  IMPLICIT NONE
  INTEGER                    :: iw, frame, layer, tout, np
  INTEGER, PARAMETER         :: forgiv=20
  
  IF (nlayers < 2) RETURN

  np = size(coarse,2)

  DO iw = 1, np
  
    frame = 1
    tout  = 0
    layer = coarse(frame,iw)

    DO WHILE ( frame <= nframes ) 

      IF ( coarse(frame,iw) /= layer ) THEN
        tout = tout + 1
      ELSEIF ( tout > 0 ) THEN 
        coarse(frame-tout:frame,iw) = layer
        tout = 0
      END IF
  
      !If the particle really changed layer
      IF ( tout > forgiv ) THEN

        layer = coarse(frame,iw)
        tout = 0

      END IF

      frame = frame + 1

    END DO

  END DO

END SUBROUTINE SMOOTH_COARSE_GRAIN

SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO
  USE parameters, ONLY : natoms, pos, box, ind_atom
  IMPLICIT NONE
  INTEGER                              :: iat, ipol
  DOUBLE PRECISION                     :: mass
  DOUBLE PRECISION                     :: mtot, m, cm(3)
 
  cm=0.d0
  mtot=0.d0
 
  DO iat = 1,natoms

    m = mass(ind_atom(iat))
    
    DO ipol = 1,3
      pos(iat,ipol) = pos(iat,ipol) & 
                    - nint(pos(iat,ipol)/box(ipol,ipol))*box(ipol,ipol)
    END DO  

    cm = cm + m*pos(iat,:)
    mtot = mtot + m

  END DO

  cm = cm / mtot
  DO iat = 1,natoms
    pos(iat,:) = pos(iat,:) - cm
  END DO

END SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO

SUBROUTINE SET_CENTER_OF_MASS_VELOCITY_TO_ZERO
  USE parameters, ONLY : natoms, vel, ind_atom
  IMPLICIT NONE
  INTEGER                              :: iat
  DOUBLE PRECISION                     :: mass
  DOUBLE PRECISION                     :: Ptot(3), m, mtot
 
  Ptot=0.d0
  mtot=0.d0
 
  DO iat = 1,natoms

    m = mass(ind_atom(iat))
    Ptot = Ptot + m*vel(iat,:)
    mtot = mtot + m

  END DO

  DO iat = 1,natoms
    vel(iat,:) = vel(iat,:) - Ptot / mtot
  END DO

END SUBROUTINE SET_CENTER_OF_MASS_VELOCITY_TO_ZERO

SUBROUTINE IDENTIFY_OH_GROUPS 
  !For each H in the system, assign the index of the O atom it is bound to.
  !Return the indexes in the array OH_index
  USE parameters, ONLY : natoms, OH_index
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2
  DOUBLE PRECISION           :: Dist, d, d_OH
  DOUBLE PRECISION           :: maxdist
  LOGICAL                    :: is_hydrogen, is_oxygen
  INTEGER                    :: nH(natoms), ind_bkp(natoms), ind_out
  
  OH_index = 0
  nH=0

  DO iat1 = 1, natoms

    IF ( is_hydrogen(iat1) ) THEN

      d = 100.d0

      DO iat2 = 1, natoms

        IF ( is_oxygen(iat2) ) THEN

          d_OH = Dist(iat1,iat2)

          IF ( d_OH < d ) THEN

            ind_bkp(iat1) = OH_index(iat1)
            OH_index(iat1) = iat2 
            d = d_OH

          END IF

        END IF

      END DO

      nH(OH_index(iat1)) = nH(OH_index(iat1)) + 1

      IF ( d > 2 ) THEN
        print *, 'Hydrogen ', iat1, ' is not bound' 
      END IF

    END IF

  END DO

  !Check if a oxygen has more than 2 hydrogens
  DO iat1 = 1,natoms

    IF ( nH(iat1) > 2 ) THEN
  
      DO iat2 = 1, natoms

        maxdist=0
        IF ( OH_index(iat2) .eq. iat1 ) THEN
          IF ( Dist(iat1,iat2) > maxdist ) THEN
            ind_out = iat2
            maxdist = Dist(iat1,iat2)
          END IF
        END IF

      END DO

      !Assign 2nd nearest neighbor to the H farthest from oxygen iat1
      OH_index(ind_out) = ind_bkp(ind_out)
      
      IF ( nH(OH_index(ind_out)) .eq. 2 ) print *, "No way for index ", ind_out

    END IF

  END DO

END SUBROUTINE IDENTIFY_OH_GROUPS

SUBROUTINE GET_IND_H ( iO, ind_H, n_h)
  !For an O atom with index iO, return the indexes of H atoms bound to it
  USE parameters, ONLY : natoms, OH_index
  INTEGER, INTENT(IN)        :: iO
  INTEGER                    :: iat, ind_H(2)
  INTEGER                    :: n_h

  ind_H = 0
  n_h = 0

  DO iat = 1, natoms

    IF ( OH_index(iat) == iO ) THEN 
      n_h = n_h + 1
      ind_H(n_h) = iat
    END IF

    IF (n_h .eq. 2) EXIT   

  END DO

END SUBROUTINE GET_IND_H

SUBROUTINE CROSSPROD( x, y, prod )
  !Cross product of 3D vectors x and y
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3)       :: prod
  DOUBLE PRECISION, INTENT(IN)         :: x(3),y(3)
  
  prod(1) = x(2)*y(3) - x(3)*y(2)
  prod(2) = x(3)*y(1) - x(1)*y(3)
  prod(3) = x(1)*y(2) - x(2)*y(1)

END SUBROUTINE CROSSPROD

SUBROUTINE NEAREST_TI( iO, d )
  !Find Ti nearest to O with index iO
  USE parameters, ONLY : pos, box, ind_atom, natoms
  INTEGER, INTENT(IN)        :: iO
  DOUBLE PRECISION           :: Dist, d(3), md
  DOUBLE PRECISION           :: mindist 
  INTEGER                    :: iTi, iat, ipol

  mindist=1000

  DO iat = 1,natoms

    IF ( ind_atom(iat) == 1 ) THEN
  
      md = Dist(iO, iat)

      IF ( md < mindist ) THEN
        iTi = iat
        mindist=md
      END IF

    END IF

  END DO

  DO ipol=1,3
    d(ipol) = pos(iTi,ipol) - pos(iO,ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol,ipol))*box(ipol,ipol)
  END DO

  d = d/mindist

END SUBROUTINE NEAREST_TI

SUBROUTINE OH_VECT(O, H, OH)
  ! Normalized OH vectors of a water molecule 
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H(2)
  DOUBLE PRECISION                     :: OH(2,3), unit_vector(3)
  INTEGER                              :: ih, ipol

  IF ( H(1) .eq. 0 ) THEN
    print *, 'ERROR: First index is zero ', O
    STOP
  END IF

  DO ih = 1,2

    IF ( H(ih) .ne. 0 ) THEN

      DO ipol=1,3

        OH(ih,ipol) = pos(H(ih),ipol) - pos(O,ipol)
        OH(ih,ipol) = OH(ih,ipol) - nint( OH(ih,ipol)/box(ipol,ipol) ) * box(ipol,ipol)

      END DO

      OH(ih,:) = OH(ih,:) / norm2(OH(ih,:))

    ELSE

      CALL NEAREST_TI(O, unit_vector) 

      OH(ih,:) = unit_vector 
 
    END IF

  END DO

END SUBROUTINE OH_VECT

SUBROUTINE DISTANCE_VECTOR( coord, wfc, d ) 
  USE parameters, ONLY : box
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: coord(3), wfc(3)
  DOUBLE PRECISION             :: d(4)
  INTEGER                      :: ipol
  
  DO ipol = 1,3
    d(ipol) = wfc(ipol) - coord(ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol,ipol))*box(ipol,ipol)
  END DO

  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR

SUBROUTINE MATINV3(A,B)
    IMPLICIT NONE
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real, intent(in) :: A(3,3)   !! Matrix
    real             :: B(3,3)   !! Inverse matrix
    real             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

END SUBROUTINE

!!!!!!!!!!!!!!! FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION layer_index(pos_z)
  !For layer resolved analysis - Get the index of Ow with z coordinate
  USE parameters, ONLY : box, nlayers, layers
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: pos_z
  DOUBLE PRECISION                     :: ind, z
  INTEGER                              :: i
 
  IF (nlayers < 2) THEN
    layer_index = 1
    RETURN
  END IF

  z = pos_z - nint( pos_z / box(3,3) ) * box(3,3)
  z = z + box(3,3)/2.d0
  i = 0 

  DO WHILE ( i < nlayers ) 
 
    i = i + 1
    IF ( ( z > layers(i) ) .and. ( z <= layers(i+1) ) ) THEN 
      layer_index = i
      RETURN
    END IF

  END DO

  layer_index = 0 
  PRINT *, "Could not find layer for particle with z ", z
  STOP

END FUNCTION layer_index

INTEGER FUNCTION layer_index_equal(pos_z)
  !For layer resolved analysis - Get the index of Ow with z coordinate
  USE parameters, ONLY : box, nlayers, layers
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN)         :: pos_z
  DOUBLE PRECISION                     :: ind, z
  INTEGER                              :: i
 
  IF (nlayers < 2) THEN
    layer_index_equal = 1
    RETURN
  END IF

  z = pos_z - nint( pos_z / box(3,3) ) * box(3,3)
  z = z + box(3,3)/2.d0

  layer_index_equal = int( nlayers * z / box(3,3) ) + 1 

END FUNCTION layer_index_equal

DOUBLE PRECISION FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  DOUBLE PRECISION                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i,i) ) * box(i,i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

DOUBLE PRECISION FUNCTION Dist_cart(ind1,ind2,icart)
  ! Distance between two points including pbc. Only 
  ! cartesian coordinate `icart` is computed
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: ind1,ind2,icart
  DOUBLE PRECISION                     :: d

  d = pos(ind1,icart) - pos(ind2,icart)
  d = d - nint( d/box(icart,icart) ) * box(icart,icart)

  Dist_cart = d

END FUNCTION Dist_cart

DOUBLE PRECISION FUNCTION dipole(O, H1, H2, ipol)
  !Water dipole along direction ipol 
  USE parameters,ONLY: box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H1, H2
  INTEGER, INTENT(IN)                  :: ipol 
  INTEGER                              :: ix
  DOUBLE PRECISION                     :: v1(3), v2(3), vout(3)

  DO ix=1,3
    v1(ix) = pos(O,ix) - pos(H1,ix)
    v1(ix) = v1(ix) - nint( v1(ix)/box(ix,ix) ) * box(ix,ix)

    v2(ix) = pos(O,ix) - pos(H2,ix)
    v2(ix) = v2(ix) - nint( v2(ix)/box(ix,ix) ) * box(ix,ix)
  END DO
  
  vout = v1 + v2

  dipole = -1.*vout(ipol) / norm2(vout) 

END FUNCTION dipole

DOUBLE PRECISION FUNCTION Angle(ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms, pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  INTEGER                    :: ipol
  DOUBLE PRECISION           :: v1(3), v2(3)

  DO ipol = 1,3

    v1(ipol) = pos(ind2,ipol) - pos(ind1,ipol)
    v1(ipol) = v1(ipol) - nint( v1(ipol) / box(ipol,ipol) ) * box(ipol,ipol)
    v2(ipol) = pos(ind3,ipol) - pos(ind1,ipol)
    v2(ipol) = v2(ipol) - nint( v2(ipol) / box(ipol,ipol) ) * box(ipol,ipol)

  END DO

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )

END FUNCTION Angle

DOUBLE PRECISION FUNCTION second_leg (x)
  !Second order legendre polynomial
  IMPLICIT NONE
  DOUBLE PRECISION                :: x 

  second_leg = ( 3d0 * x*x - 1d0 ) / 2d0

END FUNCTION second_leg

LOGICAL FUNCTION in_layer( iat, ti, tf, sep )
  !Test if atom is inside a layer from time ti to tf
  USE parameters, ONLY : coarse, nlayers
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: iat, ti, tf
  INTEGER                    :: layer, ts
  INTEGER, OPTIONAL          :: sep

  ts = 1
  IF (present(sep)) ts=sep

  IF ( nlayers > 1 ) THEN

    layer = coarse(ti,iat)
    in_layer = ALL( coarse(ti:tf:ts,iat) == layer )

  ELSE 

    in_layer = .true.

  END IF

END FUNCTION in_layer

LOGICAL FUNCTION is_oxygen(ind)
  USE parameters, ONLY : ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  IF ( ( ind_atom(ind) == 3 ) .or. ( ind_atom(ind) == 4 ) ) THEN
    is_oxygen = .true.
  ELSE
    is_oxygen = .false.
  END IF

END FUNCTION is_oxygen

LOGICAL FUNCTION is_water_oxygen(ind)
  USE parameters, ONLY : ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  IF ( ind_atom(ind) == 4 ) THEN
    is_water_oxygen = .true.
  ELSE
    is_water_oxygen = .false.
  END IF

END FUNCTION is_water_oxygen

LOGICAL FUNCTION is_hydrogen(ind)
  USE parameters, ONLY : ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  IF ( ind_atom(ind) == 2 ) THEN
    is_hydrogen = .true.
  ELSE
    is_hydrogen = .false.
  END IF

END FUNCTION is_hydrogen

LOGICAL FUNCTION is_titanium(ind)
  USE parameters, ONLY : ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  IF ( ind_atom(ind) == 1 ) THEN
    is_titanium = .true.
  ELSE
    is_titanium = .false.
  END IF

END FUNCTION is_titanium

LOGICAL FUNCTION is_TiO2_oxygen(ind)
  USE parameters, ONLY : ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  
  IF ( ind_atom(ind) == 3 ) THEN
    is_TiO2_oxygen = .true.
  ELSE
    is_TiO2_oxygen = .false.
  END IF

END FUNCTION is_TiO2_oxygen

LOGICAL FUNCTION is_hbonded( ind1, ind2 )
  !Determine if oxygen ind1 H-bonds to ind2 based on 
  !Chandler's definition of H-bond
  USE parameters, ONLY : OH_bond
  IMPLICIT NONE 
  INTEGER, INTENT(IN)                  :: ind1, ind2
  INTEGER                              :: i, index_H1(2), index_H2(2)
  DOUBLE PRECISION                     :: Angle, OHO_ang, Dist
  DOUBLE PRECISION, PARAMETER          :: cutoff_OO = 3.5d0
  
  index_H1(:) = OH_bond(ind1,:)
  index_H2(:) = OH_bond(ind2,:)

  IF ( Dist(ind1,ind2) < cutoff_OO ) THEN
    DO i = 1,2
      !H-bond donation by ind1
      OHO_ang = Angle( ind1, index_H1(i), ind2 )
      IF ( OHO_ang > 0.8660d0 ) THEN
        is_hbonded = .true.
        RETURN
      END IF
      IF ( ALL( index_H2 /= 0 ) ) THEN
        !H-bond donation to ind1
        OHO_ang = Angle( ind2, index_H2(i), ind1 )
        IF ( OHO_ang > 0.8660d0 ) THEN
          is_hbonded = .true.
          RETURN
        END IF
      END IF
    END DO
  END IF

  is_hbonded = .false. 

END FUNCTION is_hbonded

LOGICAL FUNCTION is_nn( ind1, ind2 )
  IMPLICIT NONE 
  INTEGER, INTENT(IN)                  :: ind1, ind2
  DOUBLE PRECISION                     :: Dist
  DOUBLE PRECISION, PARAMETER          :: cutoff_OO = 3.5d0
  
  IF ( Dist(ind1,ind2) < cutoff_OO ) THEN
    is_nn = .true.
  ELSE
    is_nn = .false. 
  END IF

END FUNCTION is_nn

DOUBLE PRECISION FUNCTION mass(ind_atom)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                   :: ind_atom
  
  SELECT CASE (ind_atom)
    CASE (1)
      mass = 47.867
    CASE (2)
      mass = 2.016
    CASE (3)
      mass = 15.999
    CASE (4) 
      mass = 15.999
  END SELECT

END FUNCTION mass

INTEGER FUNCTION ind_layer( ind_O )
  !Layer index for prefactor. ind_O is the index of oxygen in pos matrix
  USE parameters, ONLY : pos, box, layers, nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O
  INTEGER                    :: il
  DOUBLE PRECISION           :: z

  z = pos(ind_O,3)
  z = z - nint(z/box(3,3))*box(3,3) + box(3,3)/2.

  IF ( ( z < layers(1) ) .or. ( z > layers(nlayers+1) ) ) THEN
    ind_layer = 0
    RETURN
  END IF

  DO il = 1,nlayers

    if ( z <= layers(il+1) ) THEN
      ind_layer = il
      RETURN
    endif

  END DO

END FUNCTION ind_layer

REAL FUNCTION det_3x3(f)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: f(3,3)
  
  det_3x3 =  f(1,1)*f(2,2)*f(3,3) + f(2,1)*f(3,2)*f(1,3) + f(3,1)*f(1,2)*f(2,3) &
       - f(1,3)*f(2,2)*f(3,1) - f(2,3)*f(3,2)*f(1,1) - f(3,3)*f(1,2)*f(2,1)
  return

END FUNCTION det_3x3
