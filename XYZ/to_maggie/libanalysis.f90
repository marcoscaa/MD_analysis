
MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, nwater, corr_time, nattype
  INTEGER                            :: nframes, nequil, nhist, nwann
  INTEGER                            :: nlayers, nskip, stride, nH
  INTEGER                            :: maxneighbor, nTi5c, natoms_buff
  CHARACTER(5)                       :: atypeO,atypeH, indexC, ind_rdf(2)
  CHARACTER(5), ALLOCATABLE          :: moltype(:), atype(:), atype_wann(:) 
  INTEGER, ALLOCATABLE               :: gofr(:) 
  INTEGER, ALLOCATABLE               :: coarse(:,:) 
  INTEGER, ALLOCATABLE               :: OH_index(:)
  INTEGER, ALLOCATABLE               :: ind_atom(:)
  INTEGER, ALLOCATABLE               :: neighborlist(:,:)
  INTEGER, ALLOCATABLE               :: index_O2c(:)
  INTEGER, ALLOCATABLE               :: molid(:), index_bonded(:,:)
  REAL*8                             :: box(3,3), boxinv(3,3), zoffset
  REAL*8, ALLOCATABLE                :: pos(:,:), vel(:,:), wannier(:,:)
  REAL*8, ALLOCATABLE                :: layers(:) 
  REAL*8,ALLOCATABLE                 :: rcut(:,:)
  REAL*8                             :: dt
  REAL*8                             :: maxr, number_density, fixed_coord(3)
  REAL*8, PARAMETER                  :: cutoff_OH = 1.3d0, cutoff_NH=1.5d0
  REAL*8, PARAMETER                  :: cutoff_TiO = 2.6d0
  REAL*8, PARAMETER                  :: cutoffOwTi=2.6
  REAL*8, PARAMETER                  :: Pi=3.1415
END MODULE parameters

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+2)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE READ_XYZ
  !Read xyz file. 
  USE parameters, ONLY : pos, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  CHARACTER(5)               :: indatom
 
  READ(1,*)
  READ(1,*) !box

  DO iat = 1, natoms
    READ(1,*) atype(iat), pos(1,iat), pos(2,iat), pos(3,iat)
  END DO

END SUBROUTINE READ_XYZ

SUBROUTINE READ_EXTXYZ
  !Read extxyz file
  USE parameters, ONLY : pos, natoms, atype, box, boxinv
  IMPLICIT NONE
  INTEGER                    :: iat, il 
  REAL*8                     :: a
  CHARACTER(4)               :: typ
  CHARACTER(len=500)         :: second_line
 
  READ(1,*)
  READ(1,fmt="(A500)") second_line

  !Extract lattice
  DO iat=1,500
    IF (second_line(iat:iat+7)=='Lattice=') THEN
      DO il = 1,490-iat
        IF (second_line(il+iat+8:il+iat+8)=='"') THEN
          READ(second_line(iat+9:iat+7+il),*) box(1,1), box(1,2), box(1,3), &
                                              box(2,1), box(2,2), box(2,3), &
                                              box(3,1), box(3,2), box(3,3)
          GO TO 100
        END IF
      END DO
    END IF
  END DO

  PRINT *, "STOP: Your XYZ file does not have Lattice information"
  STOP

100  DO iat = 1, natoms
    READ(1,*) atype(iat), pos(1,iat), pos(2,iat), pos(3,iat)
  END DO

  CALL M33INV(box,boxinv)

END SUBROUTINE READ_EXTXYZ

SUBROUTINE READ_EXTXYZ_POS_VEL
  !Read extxyz file
  USE parameters, ONLY : pos, vel, natoms, atype, box, boxinv
  IMPLICIT NONE
  INTEGER                    :: iat, il 
  REAL*8                     :: a
  CHARACTER(4)               :: typ
  CHARACTER(len=500)         :: second_line
 
  READ(1,*)
  READ(1,fmt="(A500)") second_line

  !Extract lattice
  DO iat=1,500
    IF (second_line(iat:iat+7)=='Lattice=') THEN
      DO il = 1,490-iat
        IF (second_line(il+iat+8:il+iat+8)=='"') THEN
          READ(second_line(iat+9:iat+7+il),*) box(1,1), box(1,2), box(1,3), &
                                              box(2,1), box(2,2), box(2,3), &
                                              box(3,1), box(3,2), box(3,3)
          GO TO 100
        END IF
      END DO
    END IF
  END DO

100  DO iat = 1, natoms
    READ(1,*) atype(iat), pos(1,iat), pos(2,iat), pos(3,iat), &
                          vel(1,iat), vel(2,iat), vel(3,iat)
  END DO

  CALL M33INV(box,boxinv)

END SUBROUTINE READ_EXTXYZ_POS_VEL

SUBROUTINE READ_EXTXYZ_IO (iostat)
  !Read extxyz file
  USE parameters, ONLY : pos, natoms, atype, box, boxinv
  USE, intrinsic :: iso_fortran_env, Only : iostat_end
  IMPLICIT NONE
  INTEGER                    :: iat, il, iostat 
  REAL*8                     :: a
  CHARACTER(4)               :: typ
  CHARACTER(len=500)         :: second_line
 
  READ(1,*, iostat = iostat) natoms

  IF (iostat == iostat_end) THEN
    RETURN
  END IF

  READ(1,fmt="(A500)") second_line

  IF (allocated(pos)) THEN
    deallocate(pos,atype)
  END IF

  allocate(pos(3,natoms))
  allocate(atype(natoms))

  !Extract lattice
  DO iat=1,500
    IF (second_line(iat:iat+7)=='Lattice=') THEN
      DO il = 1,490-iat
        IF (second_line(il+iat+8:il+iat+8)=='"') THEN
          READ(second_line(iat+9:iat+7+il),*) box(1,1), box(1,2), box(1,3), &
                                              box(2,1), box(2,2), box(2,3), &
                                              box(3,1), box(3,2), box(3,3)
          GO TO 100
        END IF
      END DO
    END IF
  END DO

  PRINT *, "STOP: Your XYZ file does not have Lattice information"
  STOP

100  DO iat = 1, natoms
    READ(1,*,iostat=iostat) atype(iat), pos(1,iat), pos(2,iat), pos(3,iat)
  END DO

  CALL M33INV(box,boxinv)

END SUBROUTINE READ_EXTXYZ_IO

SUBROUTINE APPLY_PBC(x_in,x_out)
  USE parameters, ONLY : box,boxinv
  IMPLICIT NONE
  REAL*8, INTENT(IN)                  :: x_in(3)
  REAL*8                              :: x_out(3)
  INTEGER                             :: ipol

  x_out = MATMUL(boxinv,x_in)
  DO ipol=1,3
    x_out(ipol) = x_out(ipol) - nint(x_out(ipol))
  END DO
  x_out = MATMUL(box,x_out)

END SUBROUTINE APPLY_PBC 

SUBROUTINE DISTANCE_VECTOR( coord1, coord2, d ) 
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: coord1(3), coord2(3)
  REAL*8             :: d(4), d_(3)
  
  d_ = coord1 - coord2
  CALL APPLY_PBC(d_,d(1:3))
  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR

SUBROUTINE MINDISTANCE_VECTOR_IND( ind1, atype2, dmin, ind_min )
  USE parameters, only: natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ind1
  CHARACTER(*), INTENT(IN)     :: atype2
  INTEGER, INTENT(INOUT)       :: ind_min
  REAL*8, INTENT(INOUT)        :: dmin(3)
  INTEGER                      :: iat
  REAL*8                       :: d(4), rmin

  rmin = 100.

  DO iat=1,natoms
    IF( trim(atype(iat))==trim(atype2)) then
        CALL DISTANCE_VECTOR_IND( ind1, iat, d) 
        if (d(4)<rmin) then
            ind_min = iat
            dmin = d(1:3)
            rmin = d(4)
        end if
    END IF
  END DO

END SUBROUTINE MINDISTANCE_VECTOR_IND

SUBROUTINE DISTANCE_VECTOR_IND( ind1, ind2, d ) 
  USE parameters, ONLY : pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ind1,ind2
  REAL*8                       :: d(4)
  
  CALL DISTANCE_VECTOR(pos(:,ind1),pos(:,ind2),d)

END SUBROUTINE DISTANCE_VECTOR_IND

SUBROUTINE DIST_UNIT_VECTOR(ind1,ind2,xyz,normvec)
  ! Unit vector from ind1 to ind2 
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: ind1, ind2
  REAL*8                      :: xyz(3), d(4), normvec

  CALL DISTANCE_VECTOR_IND(ind1,ind2,d)
  normvec=d(4)
  xyz = d(1:3) / d(4) 

END SUBROUTINE DIST_UNIT_VECTOR

SUBROUTINE COARSE_GRAIN_POS (frame)
  USE parameters, ONLY : natoms,pos, coarse, nlayers, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ind, ind_z
  INTEGER                              :: layer_index, assign_o_to_h
  
  IF (nlayers < 2) RETURN

  coarse(frame,:) = 0

  DO iat = 1, natoms

    !H should have the same index as the heavy atoms it is attached to.
    IF ( trim(atype(iat))=="H" ) THEN
      ind = assign_O_to_h(iat)
    ELSE
      ind = iat
    END IF

    ind_z = layer_index(pos(3,ind),3)
    coarse(frame,iat) = ind_z

  END DO

END SUBROUTINE COARSE_GRAIN_POS

INTEGER FUNCTION CoordNumb(iat,neighbor_type,cutoff)
  USE parameters, ONLY: natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat 
  REAL*8, INTENT(IN)         :: cutoff
  CHARACTER(5), INTENT(IN)   :: neighbor_type
  INTEGER                    :: i, cn
  REAL*8                     :: Dist

  cn=0
  DO i=1,natoms
    IF (i.ne.iat) THEN
      IF (trim(atype(i))==trim(neighbor_type)) THEN
        if (Dist(iat,i)<cutoff) cn=cn+1
      end if
    END IF
  END DO

  CoordNumb=cn

END FUNCTION CoordNumb

SUBROUTINE IDENTIFY_OH_GROUPS 
  !For each H in the system, assign the index of the O atom it is bound to.
  !Return the indexes in the array OH_index
  USE parameters, ONLY : natoms, OH_index, box, atype
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2
  DOUBLE PRECISION           :: Dist, d, d_OH
  DOUBLE PRECISION           :: maxdist!, z
  INTEGER                    :: nH(natoms), ind_bkp(natoms), ind_out
  
  OH_index = 0
  nH=0

  DO iat1 = 1, natoms

    IF ( trim(atype(iat1))=="H" ) THEN

      d = 100.d0

      DO iat2 = 1, natoms

        IF ( trim(atype(iat2))=="O" ) THEN

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

END SUBROUTINE IDENTIFY_OH_GROUPS 

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

INTEGER FUNCTION CordNumb(ind,ind_type,cutoff)
  ! Number of nearest atoms with type atype 
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  CHARACTER(5), INTENT(IN)   :: ind_type
  INTEGER                    :: iat, cn
  REAL*8, INTENT(IN)         :: cutoff !angstrom
  REAL*8                     :: Dist

  cn=0

  DO iat = 1,natoms
    IF (TRIM(atype(iat))==TRIM(ind_type)) THEN
      IF (Dist(ind,iat)<=cutoff) cn=cn+1
    ENDIF
  ENDDO
 
  CordNumb=cn 

END FUNCTION CordNumb

SUBROUTINE CROSSPROD( x, y, prod )
  !Cross product of 3D vectors x and y
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: prod
  REAL*8, INTENT(IN)         :: x(3),y(3)
  
  prod(1) = x(2)*y(3) - x(3)*y(2)
  prod(2) = x(3)*y(1) - x(1)*y(3)
  prod(3) = x(1)*y(2) - x(2)*y(1)

END SUBROUTINE CROSSPROD

SUBROUTINE Bubble_Sort(a,ind,n,ntot)
  !Sort array a in ascending order. Ind is an index array
  !a will contain the smallest n entries sorted from 1 to n
  INTEGER, INTENT(in) :: n, ntot
  REAL*8, DIMENSION(ntot) :: a
  INTEGER, DIMENSION(ntot) :: ind
  REAL*8 :: temp
  INTEGER :: i, j, itemp
  LOGICAL :: swapped
 
  !DO j = SIZE(a)-1, 1, -1
  DO j = 1, n 
    swapped = .FALSE.
    !DO i = 1, j
    DO i = ntot-1, j, -1
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        itemp = ind(i)
        ind(i) = ind(i+1)
        ind(i+1) = itemp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort

! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
RECURSIVE SUBROUTINE quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
END SUBROUTINE quicksort

INTEGER FUNCTION layer_index(pos_z,zdir)
  !For layer resolved analysis - Get the index of Ow with z coordinate
  USE parameters, ONLY : box, nlayers, layers, zoffset
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos_z
  REAL*8                     :: ind, z
  INTEGER                    :: i, zdir
 
  IF (nlayers < 2) THEN
    layer_index = 1
    RETURN
  END IF

  z = pos_z + zoffset
  z = z - nint( z / box(zdir,zdir) ) * box(zdir,zdir)
  i = 0 

  DO WHILE ( i < nlayers ) 
 
    i = i + 1
    IF ( ( z > layers(i) ) .and. ( z <= layers(i+1) ) ) THEN 
      layer_index = i
      RETURN
    END IF

  END DO

  layer_index = 0 

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

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind1, ind2
  REAL*8                     :: d(4)

  CALL DISTANCE_VECTOR_IND(ind1,ind2,d)
  Dist = d(4) 
  RETURN

END FUNCTION Dist

REAL*8 FUNCTION DistXY(A,B)
  ! XY distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: A(3), B(3)
  REAL*8                     :: d(4)

  CALL DISTANCE_VECTOR(A,B,d)
  DistXY = sqrt(d(1)**2+d(2)**2) 
  RETURN

END FUNCTION DistXY

REAL*8 FUNCTION DistXYZ(A,B)
  ! XYZ distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: A(3), B(3)
  REAL*8                     :: d(4)

  CALL DISTANCE_VECTOR(A,B,d)
  DistXYZ = d(4) 
  RETURN

END FUNCTION DistXYZ

REAL*8 FUNCTION Dist_Center(center,coord)
  ! Distance from a point coord to the center
  USE parameters, ONLY :  box
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: center(3),coord(3)
  REAL*8                     :: d(4)

  CALL DISTANCE_VECTOR( coord, center, d)
  Dist_Center = d(4)

END FUNCTION Dist_Center

SUBROUTINE M33INV (A, AINV)
  !Code from David G. Simpson
  !https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt
  IMPLICIT NONE
  
  REAL*8, DIMENSION(3,3), INTENT(IN)  :: A
  REAL*8, DIMENSION(3,3), INTENT(OUT) :: AINV
  
  REAL*8, PARAMETER :: EPS = 1.0D-10
  REAL*8 :: DET
  REAL*8, DIMENSION(3,3) :: COFACTOR
  
  
  DET =   A(1,1)*A(2,2)*A(3,3)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)
  
  IF (ABS(DET) .LE. EPS) THEN
     AINV = 0.0D0
     PRINT *, "STOP: Could not invert Matrix"
     STOP
  END IF
  
  COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
  COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
  COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
  COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
  COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))
  
  AINV = TRANSPOSE(COFACTOR) / DET
  
  RETURN

END SUBROUTINE M33INV

REAL*8 FUNCTION mindist(ind, type_mindist)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: ind
  CHARACTER(5), INTENT(IN) :: type_mindist
  INTEGER                  :: iat
  REAL*8                   :: dmin, Dist, d

  dmin=100.

  DO iat=1,natoms
    IF ( (TRIM(atype(iat))==TRIM(type_mindist) ) .and. (iat.ne.ind) ) THEN
      d=Dist(ind, iat) 
      IF ( d < dmin ) THEN
        dmin=d 
      ENDIF
    ENDIF
  ENDDO

  mindist=dmin 

END FUNCTION mindist

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

REAL*8 FUNCTION Angle_reference(C, a1, reference)
  !Angle between C-a1 and reference
  !reference should be a unit vector
  !Return value in cos(Angle)
  IMPLICIT NONE
  REAL*8, INTENT(IN), DIMENSION(3)       :: C, a1, reference
  REAL*8 :: d1(4)

  CALL DISTANCE_VECTOR(C,a1,d1)
  Angle_reference = dot_product(d1(1:3),reference) / d1(4)

END FUNCTION Angle_reference

REAL*8 FUNCTION Anglei(iC, ia1, ia2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters, only : pos
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iC, ia1, ia2
  REAL*8 :: Angle

  Anglei = Angle(pos(:,iC),pos(:,ia1),pos(:,ia2)) 

END FUNCTION Anglei

REAL*8 FUNCTION Angle_Bissector(C, a1, reference)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters, only : pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)       :: C, a1(2)
  REAL*8 :: reference(3)
  REAL*8 :: d1(4), d2(4)

  CALL DISTANCE_VECTOR_IND(C,a1(1),d1)
  CALL DISTANCE_VECTOR_IND(C,a1(2),d2)
  d1 = (d1+d2)/2.
  d1(4) = sqrt(sum(d1(1:3)*d1(1:3)))
  !CALL DISTANCE_VECTOR(pos(:,C),reference,d2)
  !Angle_Bissector = dot_product(d1(1:3),d2(1:3)) / ( d1(4) * d2(4) )
  Angle_Bissector = dot_product(d1(1:3),reference) / d1(4) 

END FUNCTION Angle_Bissector

INTEGER FUNCTION get_number_atype(atyp)
    USE parameters, only : atype,natoms
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN) :: atyp
    INTEGER                 :: iat,c

    c=0
    DO iat=1,natoms
        if (trim(atype(iat))==trim(atyp)) then
            c=c+1
        end if
    END DO 

    get_number_atype=c

END FUNCTION get_number_atype 

FUNCTION M33DET (A) RESULT (DET)
    !Code by David Simpson
    !https://caps.gsfc.nasa.gov/simpson/software/m33det_f90.txt
    IMPLICIT NONE
    
    REAL*8, DIMENSION(3,3), INTENT(IN)  :: A
    REAL*8 :: DET
    
    DET =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)
    
    RETURN

END FUNCTION M33DET

INTEGER FUNCTION assign_o_to_h (indH)
  !For H with index indH, assign the index of the O atom it is bound to.
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: indH
  INTEGER                    :: indO, O_index
  DOUBLE PRECISION           :: Dist, d, d_OH
  
  O_index = 0

  d = 100.d0

  DO indO = 1, natoms

    IF ( trim(atype(indO))=="O" ) THEN

      d_OH = Dist(indH,indO)

      IF ( d_OH < d ) THEN

        O_index = indO 
        d = d_OH

      END IF

    END IF

  END DO

  IF ( O_index == 0 ) THEN
    print *, 'Hydrogen ', indH, ' is not bound' 
    O_index = indH
  END IF

  assign_o_to_h = O_index

END FUNCTION assign_o_to_h 
