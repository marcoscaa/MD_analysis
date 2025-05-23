
MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER         :: sp = selected_real_kind(6, 37)
  INTEGER                    :: natoms, nwater, nsets
  INTEGER                    :: nframes, nequil, nhist
  INTEGER                    :: nlayers, nskip, ntype
  INTEGER, ALLOCATABLE       :: moltype(:), atype(:) 
  INTEGER, ALLOCATABLE       :: gofr(:) 
  INTEGER, ALLOCATABLE       :: coarse(:,:) 
  REAL*8                     :: box(3,3), boxinv(3,3)
  REAL*8, ALLOCATABLE        :: boxt(:,:,:)
  REAL*8, ALLOCATABLE        :: pos(:,:), vel(:,:), force(:,:)
  REAL*8, ALLOCATABLE        :: post(:,:,:), velt(:,:,:)
  REAL*8, ALLOCATABLE        :: wannier(:,:,:)
  REAL*8, ALLOCATABLE        :: dipole(:,:), polar(:,:,:) 
  REAL*8, ALLOCATABLE        :: polart(:,:,:,:), dipolet(:,:,:)
  REAL*8, ALLOCATABLE        :: layers(:) 
  REAL*8                     :: dt
  REAL*8                     :: mean_pol
  REAL*8, PARAMETER          :: cutoff_OH = 1.5d0
  LOGICAL,ALLOCATABLE        :: within_cutoff(:,:,:)
  LOGICAL,ALLOCATABLE        :: at_surface(:,:)
END MODULE parameters

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil
    READ(1,*)
    READ(2,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE REMOVE_EQUIL_WANNIER
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil
    READ(1,*)
    READ(2,*)
    READ(4,*) !Wannier or Polarizability file
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL_WANNIER

SUBROUTINE READ_RAW_POS_BOX 
  !Read DPMD raw file
  USE parameters, ONLY : pos, natoms, box, boxinv
  IMPLICIT NONE
  INTEGER                    :: iat, i, j
 
  !Read Box
  READ(2,*) ( ( box(i, j), i=1,3 ), j=1,3 )
  call matinv3(box,boxinv) 

  !READ POS
  READ(1,*) ( pos(1,iat), pos(2,iat), pos(3,iat), iat=1,natoms )

END SUBROUTINE READ_RAW_POS_BOX

SUBROUTINE READ_RAW_FORCE 
  !Read DPMD raw file
  USE parameters, ONLY : force, natoms 
  IMPLICIT NONE
  INTEGER                    :: iat
 
  !READ POS
  READ(5,*) ( force(1,iat), force(2,iat), force(3,iat), iat=1,natoms )

END SUBROUTINE READ_RAW_FORCE

SUBROUTINE READ_BIN 
  USE parameters, ONLY : natoms, boxt, post, nframes
  use iso_fortran_env
  IMPLICIT NONE

  ALLOCATE(post(3,natoms,nframes))
  ALLOCATE(boxt(3,3,nframes))

  OPEN(1, file='coord.bin', action="read", form='unformatted', access='stream')
  OPEN(2, file='box.bin'  , action="read", form='unformatted', access='stream')
  READ(1) post
  READ(2) boxt

  CLOSE(1);CLOSE(2)

END SUBROUTINE READ_BIN

SUBROUTINE READ_RAW_WANNIER 
  !Read DPMD raw file
  USE parameters, ONLY : wannier, nwater
  IMPLICIT NONE
  INTEGER                    :: iw, i, j
 
  !Read Wannier centers
  READ(4,*) ( ( ( wannier(i, j, iw), i=1,3 ), j=1,4 ), iw=1,nwater )

END SUBROUTINE READ_RAW_WANNIER

SUBROUTINE READ_RAW_DIPOLE 
  !Read DPMD raw file
  USE parameters, ONLY : dipole, nwater
  IMPLICIT NONE
  INTEGER                    :: iw, i
 
  !Read Wannier centers
  READ(4,*) ( ( dipole(i, iw), i=1,3 ), iw=1,nwater )

END SUBROUTINE READ_RAW_DIPOLE

SUBROUTINE READ_RAW_DIPOLE_COMPONENTS (n) 
  !Read DPMD raw file
  USE parameters, ONLY : dipole, nwater
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: n
  INTEGER                    :: iw, i
 
  !Read Wannier centers
  READ(4,*) ( ( dipole(i, iw), i=1,n ), iw=1,nwater )

END SUBROUTINE READ_RAW_DIPOLE_COMPONENTS

SUBROUTINE READ_RAW_POLARIZABILITY 
  !Read DPMD raw file
  USE parameters, ONLY : polar, nwater
  IMPLICIT NONE
  INTEGER                    :: iw, i, j
 
  !Read local polarizabilities
  READ(5,*) ( ( ( polar(i, j, iw), i=1,3 ), j=1,3 ), iw=1,nwater )

END SUBROUTINE READ_RAW_POLARIZABILITY

SUBROUTINE READ_RAW_POLARIZABILITY_COMPONENTS (n1,n2) 
  !Read DPMD raw file
  USE parameters, ONLY : polar, nwater
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: n1,n2
  INTEGER                    :: iw, i, j
 
  !Read local polarizabilities
  READ(5,*) ( ( ( polar(i, j, iw), i=1,n1 ), j=1,n2 ), iw=1,nwater )

END SUBROUTINE READ_RAW_POLARIZABILITY_COMPONENTS

!SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO
!  USE parameters, ONLY : natoms, pos, box, atype
!  IMPLICIT NONE
!  INTEGER                              :: iat, ipol
!  REAL*8                               :: mass
!  REAL*8                               :: mtot, m, cm(3)
!  REAL*8                               :: boxinv(3,3)
! 
!  cm=0.d0
!  mtot=0.d0
!  boxinv=matinv3(box)
! 
!  DO iat = 1,natoms
!
!    m = mass(atype(iat))
!    
!    pos(:,iat) = matmul( boxinv, pos(:,iat) )
!    pos(:,iat) = pos(:,iat) - nint(pos(:,iat))
!    pos(:,iat) = matmul( box, pos(:,iat) )
!
!    cm = cm + m*pos(:,iat)
!    mtot = mtot + m
!
!  END DO
!
!  cm = cm / mtot
!  DO iat = 1,natoms
!    pos(:,iat) = pos(:,iat) - cm
!  END DO
!
!END SUBROUTINE SET_CENTER_OF_MASS_TO_ZERO

SUBROUTINE COARSE_GRAIN_POS (frame)
  USE parameters, ONLY : natoms, coarse, box, &
                         pos, nlayers
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
      ind_z = layer_index(pos(3,iat))
      coarse(frame,iwat) = ind_z

    END IF

  END DO

END SUBROUTINE COARSE_GRAIN_POS

SUBROUTINE SMOOTH_COARSE_GRAIN
  !Avoid fast jump between layers in the coarse grain function
  USE parameters, ONLY : coarse, nframes, nlayers
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

SUBROUTINE APPLY_PBC(x_in,x_out)
  USE parameters, ONLY : box, boxinv
  IMPLICIT NONE
  REAL*8, INTENT(IN)        :: x_in(3)
  REAL*8                    :: x_out(3)
  INTEGER                             :: ipol

  x_out = matmul( boxinv, x_in )
  x_out = x_out - nint(x_out)
  x_out = matmul( box, x_out)

END SUBROUTINE APPLY_PBC 

SUBROUTINE DISTANCE_VECTOR( coord, wfc, d ) 
  USE parameters, ONLY : box, boxinv
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: coord(3), wfc(3)
  REAL*8             :: d(4)

  d(1:3) = wfc - coord
  d(1:3) = matmul( boxinv, d(1:3) )
  d(1:3) = d(1:3) - nint(d(1:3))
  d(1:3) = matmul( box, d(1:3) ) 

  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR

SUBROUTINE get_H2(ind_O, ind_H)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  REAL*8                     :: Dist, dist_OH, mindist(2)
  REAL*8, PARAMETER          :: cutoff_OH = 1.6

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( atype(i) == 1 ) THEN

      dist_OH = Dist(ind_O,i)

      !Select only H's within a threshold
      IF ( dist_OH < cutoff_OH ) THEN
 
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

  IF (ANY(ind_H == 0)) THEN
    PRINT *, 'Could not find Hydrogen atoms'
    PRINT *, ind_H, ind_O
    STOP
  END IF

END SUBROUTINE get_H2

SUBROUTINE MAKE_NEIGHBOR_LIST(frame,rcut)
  USE PARAMETERS, ONLY : natoms, nwater, within_cutoff
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  REAL*8, INTENT(IN) :: rcut
  INTEGER                    :: iat, jat, iw, jw
  REAL*8                     :: d, Dist
  LOGICAL                    :: is_water_oxygen
  
  iw = 0
  jw=0

  DO iat = 1,natoms 
 
    IF (is_water_oxygen(iat)) THEN
      
      jw = iw
      iw = iw + 1

      DO jat = iat,natoms

        IF (is_water_oxygen(jat)) THEN

          jw = jw + 1

          IF (iat==jat) THEN 
            within_cutoff(iw,jw,frame) = .false.
          ELSE
            IF ( Dist(iat,jat) < rcut ) THEN
              within_cutoff(iw,jw,frame) = .true.
            ELSE
              within_cutoff(iw,jw,frame) = .false.
            ENDIF
            within_cutoff(jw,iw,frame) = within_cutoff(iw,jw,frame)
          ENDIF
              
        END IF

      END DO

    END IF
  
  END DO

END SUBROUTINE MAKE_NEIGHBOR_LIST

REAL*8 FUNCTION mass(ind_atom)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                   :: ind_atom
  
  SELECT CASE (ind_atom)
    CASE (0)
      mass = 16.000
    CASE (1)
      mass = 2.018
  END SELECT

END FUNCTION mass

LOGICAL FUNCTION is_water_oxygen(iat)
  USE parameters, only: atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: iat

  is_water_oxygen=.false.
  IF (atype(iat).eq.0) is_water_oxygen=.true.

END FUNCTION is_water_oxygen

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

  z = pos_z !- nint( pos_z / box(3) ) * box(3)
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

DOUBLE PRECISION FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box, boxinv
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                              :: i, ind1, ind2

  xyz = pos(:,ind1) - pos(:,ind2)
  xyz = matmul( boxinv, xyz )
  xyz = xyz - nint(xyz)
  xyz = matmul( box, xyz)

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Distxyz(pos1,pos2,box,boxinv)
  ! Distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos1(3),pos2(3),box(3,3),boxinv(3,3)
  REAL*8                     :: xyz(3)

  xyz = pos1 - pos2
  xyz = matmul( boxinv, xyz )
  xyz = xyz - nint(xyz)
  xyz = matmul( box, xyz)

  Distxyz = SQRT( SUM(xyz*xyz) )

END FUNCTION Distxyz

subroutine matinv3(A,B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real*8, intent(in) :: A(3,3)   !! Matrix
    real*8             :: B(3,3)   !! Inverse matrix
    real*8             :: detinv

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
  end subroutine

REAL*8 FUNCTION Angle(C, a1, a2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  IMPLICIT NONE
  REAL*8, INTENT(IN), DIMENSION(3)       :: C, a1, a2
  REAL*8 :: d1(4), d2(4)

  CALL DISTANCE_VECTOR(C,a1,d1)
  CALL DISTANCE_VECTOR(C,a2,d2)
  Angle = dot_product(d1(1:3),d2(1:3)) / ( d1(4) * d2(4) )

END FUNCTION Angle
