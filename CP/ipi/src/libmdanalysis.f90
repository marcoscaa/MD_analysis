
MODULE parameters
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE            :: pos(:,:)
  REAL*8, ALLOCATABLE            :: vel(:,:)
  REAL*8, ALLOCATABLE            :: force(:,:)
  REAL*8, ALLOCATABLE            :: layers(:)
  REAL*8                         :: box(3,3)
  REAL*8                         :: dt !Read in the index file
  REAL*8, PARAMETER              :: cutoff_OH = 1.6d0 
  REAL*8, PARAMETER              :: pi = 3.1415926539d0 
  INTEGER                                  :: nlayers 
  INTEGER, PARAMETER                       :: iprint = 10
  CHARACTER(5), ALLOCATABLE                :: atype(:)
  INTEGER, ALLOCATABLE                     :: n_water(:)
  INTEGER, ALLOCATABLE                     :: OH_bond(:,:)
  INTEGER, ALLOCATABLE                     :: OH_index(:)
  INTEGER*8, ALLOCATABLE                   :: coarse(:,:)
  INTEGER                                  :: natoms
  INTEGER                                  :: nwater 
  INTEGER                                  :: nOxygen 
  INTEGER                                  :: nframes
  INTEGER                                  :: nequil
  LOGICAL, PARAMETER                       :: inversion_symmetry=.false. 

END MODULE parameters

!!!!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_POS_CEL
  USE parameters, ONLY : pos, natoms, box, atype
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: trash

  box = 0
  READ(1,*)
  READ(1,*) trash, trash, box(1,1), box(2,2), box(3,3)

  DO i = 1,natoms
    READ(1, *) atype(i), pos(1,i), pos(2,i), pos(3,i)
    atype(i)=TRIM(atype(i))
  END DO

  box = box * 0.529177
  !pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : nframes, natoms, nequil
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+2)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : natoms, cutoff_OH, atype
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  REAL*8                     :: Dist, dist_OH, mindist(2)

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( atype(i) == 'H' ) THEN

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
  USE parameters, ONLY : natoms, OH_bond, atype
  IMPLICIT NONE
  INTEGER                    :: iat, ind_H(2)
  
  DO iat = 1, natoms

    IF ( atype(iat) == 'Ow' ) THEN
  
      CALL get_H(iat,ind_H)
      OH_bond(iat,:) = ind_H(:)

    END IF

  END DO
          
END SUBROUTINE ASSIGN_OH_BOND

SUBROUTINE COARSE_GRAIN_POS (frame)
  USE parameters, ONLY : natoms, coarse, box, &
                         pos, nlayers, nlayers, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ind_z, iwat
  INTEGER                              :: layer_index
  
  IF (nlayers < 2) RETURN

  coarse(frame,:) = 0
  iwat = 0

  DO iat = 1, natoms

    !Selecting only OW
    IF ( atype(iat)=='Ow' ) THEN

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
  USE parameters, ONLY : natoms, pos, box, atype
  IMPLICIT NONE
  INTEGER                              :: iat, ipol
  REAL*8                     :: mass
  REAL*8                     :: mtot, m, cm(3)
 
  cm=0.d0
  mtot=0.d0
 
  DO iat = 1,natoms

    m = mass(atype(iat))
    
    DO ipol = 1,3
      pos(ipol,iat) = pos(ipol,iat) & 
                    - nint(pos(ipol,iat)/box(ipol,ipol))*box(ipol,ipol)
    END DO  

    cm = cm + m*pos(:,iat)
    mtot = mtot + m

  END DO

  cm = cm / mtot
  DO iat = 1,natoms
    pos(:,iat) = pos(:,iat) - cm
  END DO

END SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO

SUBROUTINE SET_CENTER_OF_MASS_VELOCITY_TO_ZERO
  USE parameters, ONLY : natoms, vel, atype
  IMPLICIT NONE
  INTEGER                              :: iat
  REAL*8                     :: mass
  REAL*8                     :: Ptot(3), m, mtot
 
  Ptot=0.d0
  mtot=0.d0
 
  DO iat = 1,natoms

    m = mass(atype(iat))
    Ptot = Ptot + m*vel(:,iat)
    mtot = mtot + m

  END DO

  DO iat = 1,natoms
    vel(:,iat) = vel(:,iat) - Ptot / mtot
  END DO

END SUBROUTINE SET_CENTER_OF_MASS_VELOCITY_TO_ZERO

SUBROUTINE IDENTIFY_OH_GROUPS 
  !For each H in the system, assign the index of the O atom it is bound to.
  !Return the indexes in the array OH_index
  USE parameters, ONLY : natoms, OH_index, atype
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2
  REAL*8           :: Dist, d, d_OH
  REAL*8           :: maxdist
  INTEGER                    :: nH(natoms), ind_bkp(natoms), ind_out
  
  OH_index = 0
  nH=0

  DO iat1 = 1, natoms

    IF ( atype(iat1)=='H' ) THEN

      d = 100.d0

      DO iat2 = 1, natoms

        IF ( ( atype(iat2)=='O' ) .or. ( atype(iat2)=='Ow' ) ) THEN

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
  REAL*8, DIMENSION(3)       :: prod
  REAL*8, INTENT(IN)         :: x(3),y(3)
  
  prod(1) = x(2)*y(3) - x(3)*y(2)
  prod(2) = x(3)*y(1) - x(1)*y(3)
  prod(3) = x(1)*y(2) - x(2)*y(1)

END SUBROUTINE CROSSPROD

SUBROUTINE NEAREST_TI( iO, d )
  !Find Ti nearest to O with index iO
  USE parameters, ONLY : pos, box, atype, natoms
  INTEGER, INTENT(IN)        :: iO
  REAL*8           :: Dist, d(3), md
  REAL*8           :: mindist 
  INTEGER                    :: iTi, iat, ipol

  mindist=1000

  DO iat = 1,natoms

    IF ( atype(iat) == 'Ti' ) THEN
  
      md = Dist(iO, iat)

      IF ( md < mindist ) THEN
        iTi = iat
        mindist=md
      END IF

    END IF

  END DO

  DO ipol=1,3
    d(ipol) = pos(ipol,iTi) - pos(ipol,iO)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol,ipol))*box(ipol,ipol)
  END DO

  d = d/mindist

END SUBROUTINE NEAREST_TI

SUBROUTINE OH_VECT(O, H, OH)
  ! Normalized OH vectors of a water molecule 
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H(2)
  REAL*8                     :: OH(2,3), unit_vector(3)
  INTEGER                              :: ih, ipol

  IF ( H(1) .eq. 0 ) THEN
    print *, 'ERROR: First index is zero ', O
    STOP
  END IF

  DO ih = 1,2

    IF ( H(ih) .ne. 0 ) THEN

      DO ipol=1,3

        OH(ih,ipol) = pos(ipol,H(ih)) - pos(ipol,O)
        OH(ih,ipol) = OH(ipol,iH) - nint( OH(ipol,iH)/box(ipol,ipol) ) * box(ipol,ipol)

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
  REAL*8, INTENT(IN) :: coord(3), wfc(3)
  REAL*8             :: d(4)
  INTEGER                      :: ipol
  
  DO ipol = 1,3
    d(ipol) = wfc(ipol) - coord(ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol,ipol))*box(ipol,ipol)
  END DO

  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR

!!!!!!!!!!!!!!! FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION layer_index(pos_z)
  !For layer resolved analysis - Get the index of Ow with z coordinate
  USE parameters, ONLY : box, nlayers, layers
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos_z
  REAL*8                     :: ind, z
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
  REAL*8, INTENT(IN)         :: pos_z
  REAL*8                     :: ind, z
  INTEGER                              :: i
 
  IF (nlayers < 2) THEN
    layer_index_equal = 1
    RETURN
  END IF

  z = pos_z - nint( pos_z / box(3,3) ) * box(3,3)
  z = z + box(3,3)/2.d0

  layer_index_equal = int( nlayers * z / box(3,3) ) + 1 

END FUNCTION layer_index_equal

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(i,ind1) - pos(i,ind2)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i,i) ) * box(i,i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Dist_cart(ind1,ind2,icart)
  ! Distance between two points including pbc. Only 
  ! cartesian coordinate `icart` is computed
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: ind1,ind2,icart
  REAL*8                     :: d

  d = pos(icart,ind1) - pos(icart,ind2)
  d = d - nint( d/box(icart,icart) ) * box(icart,icart)

  Dist_cart = d

END FUNCTION Dist_cart

REAL*8 FUNCTION dipole(O, H1, H2, ipol)
  !Water dipole along direction ipol 
  USE parameters,ONLY: box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H1, H2
  INTEGER, INTENT(IN)                  :: ipol 
  INTEGER                              :: ix
  REAL*8                     :: v1(3), v2(3), vout(3)

  DO ix=1,3
    v1(ix) = pos(ix,O) - pos(ix,H1)
    v1(ix) = v1(ix) - nint( v1(ix)/box(ix,ix) ) * box(ix,ix)

    v2(ix) = pos(ix,O) - pos(ix,H2)
    v2(ix) = v2(ix) - nint( v2(ix)/box(ix,ix) ) * box(ix,ix)
  END DO
  
  vout = v1 + v2

  dipole = vout(ipol) / norm2(vout) 

END FUNCTION dipole

REAL*8 FUNCTION Angle(ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms, pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  INTEGER                    :: ipol
  REAL*8           :: v1(3), v2(3)

  DO ipol = 1,3

    v1(ipol) = pos(ipol,ind2) - pos(ipol,ind1)
    v1(ipol) = v1(ipol) - nint( v1(ipol) / box(ipol,ipol) ) * box(ipol,ipol)
    v2(ipol) = pos(ipol,ind3) - pos(ipol,ind1)
    v2(ipol) = v2(ipol) - nint( v2(ipol) / box(ipol,ipol) ) * box(ipol,ipol)

  END DO

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )

END FUNCTION Angle

REAL*8 FUNCTION second_leg (x)
  !Second order legendre polynomial
  IMPLICIT NONE
  REAL*8                :: x 

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

LOGICAL FUNCTION is_hbonded( ind1, ind2 )
  !Determine if oxygen ind1 H-bonds to ind2 based on 
  !Chandler's definition of H-bond
  USE parameters, ONLY : OH_bond
  IMPLICIT NONE 
  INTEGER, INTENT(IN)                  :: ind1, ind2
  INTEGER                              :: i, index_H1(2), index_H2(2)
  REAL*8                     :: Angle, OHO_ang, Dist
  REAL*8, PARAMETER          :: cutoff_OO = 3.5d0
  
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
  REAL*8                     :: Dist
  REAL*8, PARAMETER          :: cutoff_OO = 3.5d0
  
  IF ( Dist(ind1,ind2) < cutoff_OO ) THEN
    is_nn = .true.
  ELSE
    is_nn = .false. 
  END IF

END FUNCTION is_nn

REAL*8 FUNCTION mass(atype)
  IMPLICIT NONE
  CHARACTER(5),INTENT(IN)             :: atype
  
  SELECT CASE (atype)
    CASE ('Ti')
      mass = 47.867
    CASE ('H')
      mass = 2.016
    CASE ('O')
      mass = 15.999
    CASE ('Ow') 
      mass = 15.999
  END SELECT

END FUNCTION mass

INTEGER FUNCTION ind_layer( ind_O )
  !Layer index for prefactor. ind_O is the index of oxygen in pos matrix
  USE parameters, ONLY : pos, box, layers, nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O
  INTEGER                    :: il
  REAL*8           :: z

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

